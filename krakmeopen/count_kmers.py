#!/usr/bin/env python3

from collections import Counter
from stringmeup import taxonomy
from pathlib import Path
import pandas as pd
import argparse
import logging
import os.path
import pickle
import gzip

# Logging
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d [%H:%M:%S]')
log = logging.getLogger(os.path.basename(__file__))


# Argparse
parser = argparse.ArgumentParser(
    prog='count_kmers.py',
    usage='python count_kmers.py --input_file FILE [--tax_id INT | --tax_id_file FILE] --nodes FILE --names FILE --output_folder DIRECTORY',
    description='Tallies the kmers from the reads that are classified to the supplied tax_ids. Saves the data in a pickle. Input is Kraken2 read-by-read classification files (can be gzipped).')

# Inputs and outputs
parser.add_argument(
    '--input_file',
    metavar='FILE',
    type=str,
    required=True,
    help='The kraken2 read-by-read output file that this script should read from.')
parser.add_argument(
    '--output_folder',
    metavar='DIRECTORY',
    type=str,
    required=True,
    help='The output folder to place the pickled file in.')

# Supply relevant taxonomic ID on command line, or one or multiple taxonomic IDs
# through a text file.
tax_id_group = parser.add_mutually_exclusive_group(required=True)
tax_id_group.add_argument(
    '--tax_id',
    metavar='Taxonomic ID',
    type=int,
    help='A taxonomic ID for a clade that you wish to calculate quality metrics for.')
tax_id_group.add_argument(
    '--tax_id_file',
    metavar='FILE',
    type=str,
    help='Supply multiple taxonomic IDs at once. A textfile with one taxonomic ID per line. Calculate quality metrics for the clades rooted at the taxonomic IDs in the file.')

# The taxonomy
parser.add_argument(
    '--names',
    metavar='FILE',
    required=True,
    help='taxonomy names dump file (names.dmp)')
parser.add_argument(
    '--nodes',
    metavar='FILE',
    required=True,
    help='taxonomy nodes dump file (nodes.dmp)')

# Parse the arguments
args = parser.parse_args()


## FUNCTIONS
def read_file(filename):
    """
    Wrapper to read either gzipped or ordinary text file input.
    """
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')

def get_tax_ids(args):
    """
    Returns a set containing ints of taxonomic IDs. Reads either from a file
    with potentially multiple taxonomic IDs, or a single taxonomic ID from
    the command line.
    """
    tax_id_set = set()

    # Loop over each line in the supplied file and add the taxonomic IDs to
    # the list
    if args.tax_id_file:
        log.info('Reading taxonomic IDs from {}...'.format(args.tax_id_file))
        with read_file(args.tax_id_file) as f:
            for l in f:
                l = int(l.strip())
                tax_id_set.add(l)
            log.info('Found {} taxonomic IDs.'.format(len(tax_id_set)))

    # If no file, taxonomic ID must come from the command line
    else:
        log.info('Using a single taxonomic ID from the command line: {}.'.format(args.tax_id))
        tax_id_set.add(args.tax_id)

    return tax_id_set

def get_tax_id_maps(tax_id_set, taxonomy_tree):
    """
    Create a dictionary mapping of all children (and root nodes) in all clades
    rooted at the tax_ids in the tax_id_set. The keys are the tax_ids in the
    clades, the values are the clade root node's tax_id.

    Create another dictionary where the keys are the clade tax_ids from
    tax_id_set, and the values are sets containing the tax_ids for the members
    of the specific clade.
    """

    log.info('Gathering all nodes in all clades rooted at supplied taxIDs...')

    # Main dictionary mappings
    children2clade_roots_map = {}
    clade_roots2children_map = {}

    for tax_id in tax_id_set:
        # Get the tax_ids that are in the clade rooted at tax_id
        clade_tax_ids = taxonomy_tree.get_clade([tax_id])[tax_id]

        # Create the children2root dictionary mapping for the current clade
        clade_children2root_map = {
            clade_tax_id: tax_id for clade_tax_id in clade_tax_ids}

        # Update the main dictionary mappings to contain the tax_ids from
        # current clade
        children2clade_roots_map.update(clade_children2root_map)

        # Update the clade roots to members dictionary mapping
        clade_roots2children_map[tax_id] = clade_tax_ids

    # Done
    log.info('Done.')
    return clade_roots2children_map, children2clade_roots_map

def process_kmer_string(kmer_info_string, paired_input):
    """
    Process a kmer info string (last column of a Kraken 2 output file), so that
    we get a count of total number of kmers and number of ambiguous kmers.
    """
    kmer_info_string = kmer_info_string.split()

    # Kraken2 classifications file for paired data contain the "|:|" delimiter
    if paired_input:
        kmer_info_string.remove('|:|')

    # Messy list comprehension. Converts all "taxa":"num_kmer" string pairs
    # into integer tuples like (taxa, num_kmers), and saves them in a list.
    kmer_classifications = [
        (int(x[0]), int(x[1])) for x in (
            kmer_info.split(':') for kmer_info in kmer_info_string)
        if x[0] != 'A']

    # Further processes the (taxa, num_kmers) tuples into a dict where each
    # tax_id stores the total sum of kmer hits to that tax_id.
    taxa_kmer_dict = {}
    for kmer_info in kmer_classifications:
        if kmer_info[0] not in taxa_kmer_dict:
            taxa_kmer_dict[kmer_info[0]] = kmer_info[1]
        else:
            taxa_kmer_dict[kmer_info[0]] += kmer_info[1]

    return taxa_kmer_dict

def collect_kmers_per_clade(input_file, supplied_tax_ids_set, taxonomy_tree, clade_tax_id_map):
    """
    Parse an input file. For all clades rooted at the tax_ids in
    supplied_tax_ids_set, find all reads that are classified there and count
    the kmers and where they hit. Will create a dict of Counters that have
    the clade root tax_ids as keys. The Counters will in turn contain keys for
    all tax_ids that kmers have hit, and values that represent how many kmers
    hit the specific tax_id.
    returns {clade_root_tax_id: Counter{tax_id: NUM_KMERS}}
    """
    report_frequency = 2500000

    # The main datastructure that we will populate
    # Using the Counter class from collections. When updating it, it will
    # automatically add together the values of the matching keys.
    taxa_kmer_hit_dict = {tax_id: Counter() for tax_id in supplied_tax_ids_set}

    # Open the classifications file and parse it, one file at a time
    with read_file(input_file) as f:

        # Setup
        first_read = True
        mhg_present = False   # Is minimizer hit groups present as a column?
        paired_input = False  # Are the reads paired?
        i = 0

        # Parse the reads in the file, one by one
        for r in f:

            # Skip if the read is unclassified
            if r.startswith('U'):
                i += 1
                if i % report_frequency == 0:
                    log.info('Processed {} reads...'.format(i))
                continue

            # Enter here only once per input file
            if first_read:
                first_read = False
                r_split = r.strip().split('\t')

                # Naive check if there's a column in the input file for minimizer_hit_groups
                # See https://github.com/danisven/StringMeUp
                if len(r_split) == 6:
                    mhg_present = True

                # Paired reads?
                if "|" in r_split[3]:
                    paired_input = True

            # Split the read classification line
            r_split = r.strip().split('\t')

            # The tax_id that the read is classified to
            tax_id_classification = int(r_split[2])

            # Is the tax_id in any of the clades rooted at the supplied tax_ids?
            if tax_id_classification in clade_tax_id_map:

                # The clade root tax_id
                tax_id_root = clade_tax_id_map[tax_id_classification]

                # The kmer string
                kmer_string = r_split[-1]

                # The kmer dict (from the kmer string of the input line)
                kmer_dict = process_kmer_string(kmer_string, paired_input)

                # Add the read's kmer information to the main datastructure
                taxa_kmer_hit_dict[tax_id_root].update(kmer_dict)

            i += 1
            if i % report_frequency == 0:
                log.info('Processed {} reads...'.format(i))

    log.info('Finished processing file. Read {} reads in total.'.format(i))
    return taxa_kmer_hit_dict


### CODE
log.info('Going to tally the kmers from the reads that are classified to your supplied taxonomic IDs.')

# Get the tax IDs that we should use
supplied_tax_ids_set = get_tax_ids(args)

# Create a taxonomy tree
tt = taxonomy.TaxonomyTree(
    names_filename=args.names,
    nodes_filename=args.nodes)

# Get a mapping of all clade tax_ids to their respective clade root
clade2tax_id_map, tax_id2clade_map = get_tax_id_maps(supplied_tax_ids_set, tt)

# Logging
log.info('Processing file {}...'.format(
    os.path.split(args.input_file)[1]))


# Get the kmer information from the input file
sample_kmer_collection = collect_kmers_per_clade(
    args.input_file, supplied_tax_ids_set, tt, tax_id2clade_map)

# Output file name
filename = Path(args.input_file)
base_filename_prefix = filename.stem.split('.')[0]
output_file_name = base_filename_prefix + '_kmers.pickle'
out_file_path = os.path.join(args.output_folder, output_file_name)

# Pickle the sample_kmer_collection dictionary
log.info('Saving to {}...'.format(out_file_path))
pickle.dump(sample_kmer_collection, open(out_file_path, 'wb'))
log.info('Done. Exiting.')

#!/usr/bin/env python3

from collections import Counter
from stringmeup import taxonomy
import pandas as pd
import argparse
import logging
import os.path
import pickle
import gzip
import glob
import sys


# Logging
logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d [%H:%M:%S]')
log = logging.getLogger(os.path.basename(__file__))


# Argparse
parser = argparse.ArgumentParser(
    prog='classification_quality_metrics.py',
    usage='python classification_quality_metrics.py --input_pattern "FILENAME_PATTERN" [--tax_id INT | --tax_id_file FILE] --output_prefix PREFIX --nodes FILE --names FILE --supporting_dir DIRECTORY',
    description='Calculates a number of classification quality metrics at the genus level for Kraken2 classifications, across multiple samples. Input is Kraken2 read-by-read classification files (can be gzipped).')

# Inputs and outputs
parser.add_argument(
    '--supporting_dir',
    metavar='DIRECTORY',
    type=str,
    required=True,
    help='Directory where supporting files are located.')
parser.add_argument(
    '--input_pattern',
    metavar='PATTERN',
    type=str,
    required=True,
    help='The filename pattern of the pickled kmer collection files created with count_kmers.py. Use quotes! For example "Ki*_kmers.pickle"')
parser.add_argument(
    '--output_prefix',
    metavar='PREFIX',
    type=str,
    required=True,
    help='The output file prefix.')

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


### SETUP

# Supporting files
metadata_file = os.path.join(args.supporting_dir, 'master_db_191220_minimizers_sequenceLengths_taxonomy.tsv.gz')
sample_key_file = os.path.join(args.supporting_dir, 'new_sample_name_key.csv')

# Find the report files
pattern = args.input_pattern
kmer_collection_files = [f for f in glob.glob(pattern)]

# Logging and input check
if len(kmer_collection_files) == 0:
    log.info(
        'Couldn\'t find any input files. Check your supplied pattern ({}). Exiting.'.format(pattern))
    sys.exit()
else:
    log.info('Found {} files matching the supplied pattern.'.format(
        len(kmer_collection_files)))


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

def load_kmer_collection_pickle(pickle_file):
    """
    Read a pickle file created with count_kmers.py.
    """
    kmer_collection = pickle.load(open(pickle_file, 'rb'))
    return kmer_collection

def get_clade_kmer_metrics(clade_id, taxonomy_tree, clade2members_dict, kmer_collection):
    """
    Function to calculate a number of metrics for the supplied clade_id. Uses
    only the number of kmers that hit throughout the taxonomy, not reads.
    """

    # A dictionary that will be populated with the different metrics
    metrics_dict = {
        'confidence_original': None,
        'confidence_classified': None,
        'other_kmers_lineage_ratio': None,
        'other_kmers_root_ratio': None,
        'other_kmers_classified_ratio': None,
        'other_kmers_distance': None,
        'other_kmers_distance_lineage_excluded': None}

    # Get the relevant kmers
    clade_kmer_collection = kmer_collection[clade_id]

    # If the tax_id has no reads classified to it, it will have no kmers to
    # operate with. Return.
    if len(clade_kmer_collection) == 0:
        return metrics_dict

    # Get the tax_ids for the members of the clade rooted at clade_id
    clade_members = clade2members_dict[clade_id]

    # The lineage of the clade, make a set of it
    #lineage_full = tt.get_lineage([clade_id])[clade_id]
    #lineage = set([x for x in lineage_full])  # make a copy of the lineage and turn into a set
    lineage = set(tt.get_lineage([clade_id])[clade_id])
    lineage.remove(clade_id) # remove the clade_id from the lineage (only want tax_ids above the clade tax_id)

    # Dictionary comprehension and set intersection to extract the
    # <tax_id>: <num_kmers> for members of the clade
    clade_kmers = {
        key: clade_kmer_collection[key]
        for key in clade_members & clade_kmer_collection.keys()}

    # Dictionary comprehension and set complement to extract the
    # <tax_id>: <num_kmers> for non-members of the clade
    other_kmers = {
        key: clade_kmer_collection[key]
        for key in clade_kmer_collection.keys() - clade_members}

    # Dictionary comprehension and set intersection to extract the
    # <tax_id>: <num_kmers> for tax_ids in the lineage
    lineage_kmers = {
        key: clade_kmer_collection[key]
        for key in lineage & clade_kmer_collection.keys()}

    # Sum the number of kmers for each class of kmer
    num_clade_kmers = sum(clade_kmers.values())
    num_other_kmers = sum(other_kmers.values())
    num_lineage_kmers = sum(lineage_kmers.values())
    total_kmers = num_clade_kmers + num_other_kmers

    # The ratio of kmers hitting the root over the number of kmers not hitting the clade
    if 1 in clade_kmer_collection:
        num_root_kmers = clade_kmer_collection[1]
        other_kmers_root_ratio = num_root_kmers / num_other_kmers
    else:
        other_kmers_root_ratio = 0.0

    # The number of kmers that have been classified to a tax_id, but not in the clade
    if 0 in other_kmers:
        num_other_kmers_classified = num_other_kmers - other_kmers[0]
    else:
        num_other_kmers_classified = num_other_kmers

    # The two different confidence scores
    confidence_original = num_clade_kmers / total_kmers
    confidence_classified = num_clade_kmers / (num_other_kmers_classified + num_clade_kmers)

    # Kmers hitting the lineage / kmers not hitting the clade members
    other_kmers_lineage_ratio = num_lineage_kmers / num_other_kmers if num_other_kmers > 0 else None

    # For the kmers not hitting in the clade, the ratio of kmers that hit a
    # tax_id over those that are not classified
    other_kmers_classified_ratio = num_other_kmers_classified / num_other_kmers if num_other_kmers > 0 else None

    # Getting distance metrics
    distance_dict = {tax_id: None for tax_id in other_kmers if tax_id != 0}  # We can't calculate a distance to 'Unclassified' (tax_id = 0)
    for tax_id in distance_dict:

        # The tax_id is an ancestor of the clade root
        if tax_id in lineage:

            # The distance between the clade root and the ancestor
            distance = tt.get_distance(tax_id, clade_id)

        # Not an ancestor, must compute two distances and add them together
        else:

            # Get the lineage of the tax_id
            tax_id_lineage = tt.get_lineage([tax_id])[tax_id]
            tax_id_lineage.reverse()  # Flip the lineage so that it goes from leaf to root

            # Loop to find the lowest common ancestor (lca) of the clade id and
            # the tax_id that we are currently getting the distance to
            lca = None
            for ancestor in tax_id_lineage:
                if ancestor in lineage:
                    lca = ancestor
                    break

            # Find the distance between the lca and the clade id
            clade_lca_distance = tt.get_distance(ancestor, clade_id)

            # Find the distance between the lca and the tax_id currently
            # being investigated
            tax_id_lca_distance = tt.get_distance(ancestor, tax_id)

            # The distance between the clade id and the tax_id is the sum of
            # the two distances
            distance = clade_lca_distance + tax_id_lca_distance

        # Save the distance between the clade id and the tax_id
        distance_dict[tax_id] = distance

    # Multiply the distances with the number of kmers for each tax_id
    weighted_distance_dict = {
        tax_id: distance_dict[tax_id] * other_kmers[tax_id]
        for tax_id in distance_dict.keys()}

    # Average distance to kmers not hitting within the clade
    other_kmers_distance = sum(
        weighted_distance_dict.values()) / num_other_kmers_classified if num_other_kmers_classified > 0 else None

    # Average distance to kmers not hitting within the clade,
    # nor in the lineage of the clade
    weighted_distance_dict_excl_lineage = {
        key: weighted_distance_dict[key]
        for key in distance_dict.keys() - lineage}
    other_kmers_distance_lineage_excluded = sum(
        weighted_distance_dict_excl_lineage.values()) / (
            num_other_kmers_classified - num_lineage_kmers) if num_other_kmers_classified - num_lineage_kmers > 0 else None  # If no other_kmers have hit outside of lineage, we want to avoid dividing by 0

    # Save the values of all variables
    metrics_dict = {
        'confidence_original': confidence_original,
        'confidence_classified': confidence_classified,
        'other_kmers_lineage_ratio': other_kmers_lineage_ratio,
        'other_kmers_root_ratio': other_kmers_root_ratio,
        'other_kmers_classified_ratio': other_kmers_classified_ratio,
        'other_kmers_distance': other_kmers_distance,
        'other_kmers_distance_lineage_excluded': other_kmers_distance_lineage_excluded}

    # Done
    return metrics_dict

def get_sample_name(filename):
    """
    Returns
        1) a sample name made from components of the input file name: YYYY_WW
        2) the year
    """
    year = filename.split('-')[1]
    week = filename.split('-')[2]
    sample_name = '_'.join([year, week])
    return sample_name, year


### CODE
# Get the tax IDs that we should use
supplied_tax_ids_set = get_tax_ids(args)

# Create a taxonomy tree
tt = taxonomy.TaxonomyTree(
    names_filename=args.names,
    nodes_filename=args.nodes)

## Metadata for the supplied tax_ids
# Read metadata file
log.info('Loading metadata information...')
meta_df = pd.read_csv(metadata_file, compression='gzip', sep='\t')

# Create an old_sample_names to new_sample_names mapping
new_sample_name = {}
with read_file(sample_key_file) as f:
    i = 1
    for line in f:

        # Skip header
        if i == 1:
            i += 1
            continue

        original_name = line.split(';')[0].strip()
        new_name = line.split(';')[1].strip()

        new_sample_name[original_name] = new_name

        i += 1

log.info('Done.')

# Get a mapping of all clade tax_ids to their respective clade root
clade2tax_id_map, tax_id2clade_map = get_tax_id_maps(supplied_tax_ids_set, tt)

# Preparations
samples_kmer_collection_dict = {}
samples_in_years = {}
i = 1

# Load the pickled kmer collections
log.info('Loading input files...')
for file_name in kmer_collection_files:
    # Logging
    log.info('Loading file {} ({}/{}).'.format(
        os.path.split(file_name)[1], i, len(kmer_collection_files)))

    # The name and year of the sample
    sample_name, year = get_sample_name(file_name)

    # Load the pickled kmer collection
    sample_kmer_collection = load_kmer_collection_pickle(file_name)

    # Add the data on a per-sample basis to the samples_kmer_collection_dict
    samples_kmer_collection_dict[sample_name] = sample_kmer_collection

    # Keep track of which year the sample belongs to
    if year in samples_in_years:
        samples_in_years[year].add(sample_name)
    else:
        samples_in_years[year] = set([sample_name])

    # Increment
    i += 1

# Done
log.info('Done loading the input files.')

## Calculate metrics and output results
# Weekly metrics
log.info('Calculating sample-wise metrics...')
sample_metrics_dict = {}
i = 1
for sample in samples_kmer_collection_dict.keys():
    # Logging
    log.info('Processing sample {} ({}/{})...'.format(
        sample, i, len(samples_kmer_collection_dict.keys())))

    # Calculate the sample-wise metrics and store them in the metrics_dict, keys are clade_ids
    metrics_dict = {clade_id: None for clade_id in supplied_tax_ids_set}

    for clade_id in supplied_tax_ids_set:

        # Calculate several metrics for the kmers of this clade
        clade_metrics_dict = get_clade_kmer_metrics(
            clade_id, tt, clade2tax_id_map, samples_kmer_collection_dict[sample])

        # Add the metrics to the main dict
        metrics_dict[clade_id] = clade_metrics_dict

    # Make a dataframe from the metrics_dict
    df = pd.DataFrame.from_dict(metrics_dict, orient='index')
    df = df.reset_index()
    df = df.rename(columns={'index': 'tax_id'})

    # Save in sample_metrics_dict using an updated sample name, if there is one
    if sample in new_sample_name:
        sample_name = new_sample_name[sample]
    else:
        sample_name = sample
    sample_metrics_dict[sample_name] = df

    # Increment
    i += 1

log.info('Done.')

a_sample = [sample for sample in sample_metrics_dict.keys()][0]
metric_names = [metric for metric in sample_metrics_dict[a_sample].columns if metric != 'tax_id']

# Iterate over the different metrics, concatenate the specific metric from all
# samples, and output to file
log.info('Saving weekly metrics to files...')
i = 1
for metric in metric_names:
    # Logging
    log.info('Processing metric {} ({}/{})...'.format(
        metric, i, len(metric_names)))

    # Dict of Series
    metric_dict = {}

    # For each sample, extract the metric column of interest and save it as a
    # series in metric_dict
    for sample in sample_metrics_dict.keys():
        df = sample_metrics_dict[sample]
        df = df.set_index('tax_id')
        metric_series = df[metric]
        metric_dict[sample] = metric_series

    # Concatenate the metric_dict into a DF
    df = pd.concat(metric_dict, axis=1)

    # Sort the sample columns
    df = df[sorted(df.columns)]

    # Add metadata columns
    df = df.reset_index()
    df = meta_df.merge(df, on='tax_id')

    # Save to file
    file_name = '_'.join([args.output_prefix, 'weekly', metric]) + '.tsv'
    df.to_csv(file_name, sep='\t', index=False)

    # Increment
    i += 1

log.info('Done.')

log.info('Calculating year-wise metrics...')
year_metrics_dict = {}
i = 1
for year in samples_in_years.keys():
    # Logging
    log.info('Processing year {} ({}/{})...'.format(
        year, i, len(samples_in_years.keys())))

    # The samples that belong to this year
    year_samples = samples_in_years[year]

    # Sum up all kmers for all tax_ids from the samples that belong to the
    # current year
    year_kmers = {tax_id: Counter() for tax_id in supplied_tax_ids_set}
    for sample in year_samples:
        for tax_id in supplied_tax_ids_set:
            year_kmers[tax_id].update(samples_kmer_collection_dict[sample][tax_id])

    # Calculate the year-wise metrics and store them in the metrics_dict, keys
    # are clade_ids
    metrics_dict = {clade_id: None for clade_id in supplied_tax_ids_set}

    for clade_id in supplied_tax_ids_set:

        # Calculate several metrics for the kmers of this clade
        clade_metrics_dict = get_clade_kmer_metrics(
            clade_id, tt, clade2tax_id_map, year_kmers)

        # Add the metrics to the main dict
        metrics_dict[clade_id] = clade_metrics_dict

    # Make a dataframe from the metrics_dict
    df = pd.DataFrame.from_dict(metrics_dict, orient='index')
    df = df.reset_index()
    df = df.rename(columns={'index': 'tax_id'})

    # Save in year_metrics_dict
    year_metrics_dict[year] = df

    # Increment
    i += 1

log.info('Done.')

# Iterate over the different metrics, concatenate the specific metric from all
# years, and output to file
log.info('Saving year-wise metrics to file...')
i = 1
for metric in metric_names:
    # Logging
    log.info('Processing metric {} ({}/{})...'.format(
        metric, i, len(metric_names)))

    # Dict of Series
    metric_dict = {}

    # For each year, extract the metric column of interest and save it as a
    # series in metric_dict
    for year in samples_in_years.keys():
        df = year_metrics_dict[year]
        df = df.set_index('tax_id')
        metric_series = df[metric]
        metric_dict[year] = metric_series

    # Concatenate the metric_dict into a DF
    df = pd.concat(metric_dict, axis=1)

    # Sort the year columns
    df = df[sorted(df.columns)]

    # Add metadata columns
    df = df.reset_index()
    df = meta_df.merge(df, on='tax_id')

    # Save to file
    file_name = '_'.join([args.output_prefix, 'yearly', metric]) + '.tsv'
    df.to_csv(file_name, sep='\t', index=False)

    # Increment
    i += 1

log.info('Done.')

log.info('Calculating a grand total metrics for all samples...')
# Sum up all kmers for all tax_ids from all samples
total_kmers = {tax_id: Counter() for tax_id in supplied_tax_ids_set}
for sample in samples_kmer_collection_dict.keys():
    for tax_id in supplied_tax_ids_set:
        total_kmers[tax_id].update(samples_kmer_collection_dict[sample][tax_id])

# Calculate the grand total metrics and store them in the metrics_dict, keys
# are clade_ids
metrics_dict = {clade_id: None for clade_id in supplied_tax_ids_set}

for clade_id in supplied_tax_ids_set:

    # Calculate several metrics for the kmers of this clade
    clade_metrics_dict = get_clade_kmer_metrics(
        clade_id, tt, clade2tax_id_map, total_kmers)

    # Add the metrics to the main dict
    metrics_dict[clade_id] = clade_metrics_dict

# Make a dataframe from the metrics_dict
df = pd.DataFrame.from_dict(metrics_dict, orient='index')
df = df.reset_index()
df = df.rename(columns={'index': 'tax_id'})
log.info('Done.')

log.info('Saving the grand total metrics to file...')

# Add metadata columns
df = df.reset_index()
df = meta_df.merge(df, on='tax_id')

# Save to file
file_name = '_'.join([args.output_prefix, 'grand_total']) + '.tsv'
df.to_csv(file_name, sep='\t', index=False)

log.info('Done.')

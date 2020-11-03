#!/usr/bin/env python3

import logging
import pathlib
import pickle
from collections import Counter
from utility_functions import read_file
from stringmeup.taxonomy import TaxonomyTree

# Set up logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
format = '%(asctime)s: %(name)s: %(levelname)s:    %(message)s'
date_format = '%Y-%m-%d [%H:%M:%S]'
formatter = logging.Formatter(fmt=format, datefmt=date_format)
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)


class KmerCounterFaultyArgumentsError(Exception):
    """Raised when a wrong or missing argument is given at instantiation."""
    pass

class KmerCounterNotImplementedError(Exception):
    """Raised when trying to use functionality that is not implemented yet."""
    pass


class KmerCounter:

    def __init__(self, tax_id_set, names=None, nodes=None, taxonomy_tree=None):
        self.taxonomy_tree = self.get_taxonomy(taxonomy_tree, names, nodes)
        self.tax_id_set = tax_id_set
        self.nfiles_processed = 0

        # Dictionary mappings of clade root taxIDs to children,
        # and clade children to clade root taxIDs
        self.clade_roots2children_map, \
        self.children2clade_roots_map = self.get_tax_id_mappings()

        # Main data structure to hold the kmer counts for each taxon.
        # Using the Counter class from collections. When updating it, it will
        # automatically add together the values of the matching keys.
        self.kmer_tallies = {
            tax_id: Counter() for tax_id in self.tax_id_set
        }

    def tally_kmers(self, classifications_path, output_pickle_path=None):
        """
        Function to count kmers and return a final result.

        Takes a kraken2 classifications filename (pathlib.Path).

        Optionally takes a filename (pathlib.Path) to pickle the kmer counts
        from the given classifications file to.
        """

        # Check inputs
        files = [path for path in [classifications_path, output_pickle_path]
                 if path is not None]
        for putative_path in files:
            if not isinstance(putative_path, pathlib.Path):
                msg = ('Paths to files must be pathlib.Path objects.')
                raise KmerCounterFaultyArgumentsError(msg)

        logger.info('Processing file {}...'.format(classifications_path.name))

        # When adding kmer counts from multiple files
        if self.nfiles_processed > 0:
            raise KmerCounterNotImplementedError('You tried to add multiple results.')
            logger.info('Adding kmer counts from an additional file.')
            logger.debug('Creating a copy of main self.kmer_tallies ' + \
                         'to compare with after counting this file.')

        # Count the kmers
        self.tally_per_clade(classifications_path)

        if output_pickle_path:
            raise KmerCounterNotImplementedError('You tried to save output.')
            logger.info('Pickling result to {}'.format(output_pickle))

        self.nfiles_processed += 1

    def get_tally(self):
        """
        Getter for self.kmer_tallies.
        """

        if self.nfiles_processed > 0:
            return self.kmer_tallies
        else:
            logger.warning('Have not counted kmers yet. Run tally_kmers first.')
            return None

    def tally_per_clade(self, input_file):
        """
        Parse the input file. For all clades rooted at the tax_ids in
        self.tax_id_set, find all reads that are classified there and count
        the kmers and where they hit.

        Saves info in the main datastructure (self.kmer_tallies), which is a
        dict of Counters that have the clade root tax_ids as keys. The Counters
        will in turn contain keys for all tax_ids that kmers have hit, and values
        that represent how many that kmers hit the specific tax_id.
        """

        # Log message every Nth parsed line of the input file
        report_frequency = 2500000

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

                    # Naive check if there's a column in the input file for
                    # minimizer_hit_groups. See https://github.com/danisven/StringMeUp
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
                if tax_id_classification in self.children2clade_roots_map:

                    # The clade root tax_id
                    tax_id_root = self.children2clade_roots_map[tax_id_classification]

                    # The kmer string
                    kmer_string = r_split[-1]

                    # The kmer dict (from the kmer string of the input line)
                    kmer_dict = self.process_kmer_string(kmer_string, paired_input)

                    # Add the read's kmer information to the main datastructure
                    self.kmer_tallies[tax_id_root].update(kmer_dict)

                i += 1
                if i % report_frequency == 0:
                    logger.info('Processed {} reads...'.format(i))

        logger.info('Finished processing file. Read {} reads in total.'.format(i))

    def process_kmer_string(self, kmer_info_string, paired_input):
        """
        Process a kmer info string (last column of a Kraken 2 output file),
        so that we get a count of total number of kmers and number of
        ambiguous kmers.
        """

        kmer_info_string = kmer_info_string.split()

        # Kraken2 classifications file for paired data contain the "|:|" delimiter
        if paired_input:
            kmer_info_string.remove('|:|')

        # Convert all "taxa":"num_kmer" string pairs into integer tuples
        # like (taxa, num_kmers), and save them in a list.
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

    def get_tax_id_mappings(self):
        """
        Create a dictionary mapping of all children (and root nodes) in all clades
        rooted at the tax_ids in the tax_id_set. The keys are the tax_ids in the
        clades, the values are the clade root node's tax_id.

        Create another dictionary where the keys are the clade tax_ids from
        tax_id_set, and the values are sets containing the tax_ids for the members
        of the specific clade.
        """

        logger.info('Gathering all nodes in all clades rooted at supplied taxIDs...')

        # Main dictionary mappings
        children2clade_roots_map = {}
        clade_roots2children_map = {}

        for tax_id in self.tax_id_set:
            # Get the tax_ids that are in the clade rooted at tax_id
            clade_tax_ids = self.taxonomy_tree.get_clade([tax_id])[tax_id]

            # Create the children2root dictionary mapping for the current clade
            clade_children2root_map = {
                clade_tax_id: tax_id for clade_tax_id in clade_tax_ids}

            # Update the main dictionary mappings to contain the tax_ids from
            # current clade
            children2clade_roots_map.update(clade_children2root_map)

            # Update the clade roots to members dictionary mapping
            clade_roots2children_map[tax_id] = clade_tax_ids

        logger.info('Done.')

        return clade_roots2children_map, children2clade_roots_map


    def get_taxonomy(self, taxonomy_tree_input, names_input, nodes_input):
        """
        Use supplied TaxonomyTree instance, or create a new one from
        the supplied names and nodes.
        """

        taxonomy_tree = None

        if taxonomy_tree_input is None:
            if (names_input is None) or (nodes_input is None):
                msg = ("If not providing an instance of TaxonomyTree, " +
                       "you need to supply both names and nodes files.")
                raise KmerCounterFaultyArgumentsError(msg)

            else:
                logger.info("Building taxonomy tree.")
                taxonomy_tree = TaxonomyTree(
                    names_filename=names_input,
                    nodes_filename=nodes_input,
                )

        elif isinstance(taxonomy_tree_input, TaxonomyTree):
            taxonomy_tree = taxonomy_tree_input

        else:
            msg = ("The class of the supplied taxonomy_tree argument was not TaxonomyTree, " +
                   "instead it was {}.".format(type(taxonomy_tree_input)))
            raise KmerCounterFaultyArgumentsError(msg)

        return taxonomy_tree

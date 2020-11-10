#!/usr/bin/env python3

import logging
import pathlib
import pickle
from collections import Counter
from krakmeopen.utilities import read_file, vet_input_path, vet_output_path, vet_kraken2_file
from stringmeup.taxonomy import TaxonomyTree

logger = logging.getLogger(__name__)


class KmerCounterException(Exception):
    """Raised when an error occurs in KmerCounter."""
    pass


class KmerCounter:
    """
    Tallies kmers from Kraken2 classifications. Tallies all kmers from all
    reads that are classified to any tax_id within the clades rooted at the
    supplied tax IDs.

    Creates a datastructure like so:
        {
            clade_root_id: collections.Counter{
                tax_id: n_kmers
            }
        }
    where clade_root_id are the tax IDs supplied to the __init__ function in
    the tax_id_set. And the n_kmers (values) within the Counter objects are
    the number of kmers that hit the tax_id (keys). N.B. the number of kmers
    are tallied from all reads.

    Usage:
    kmer_counter = KmerCounter(
        tax_id_set,
        clade_roots2children_map,
        children2clade_roots_map)      # Initialize
    kmer_counter.tally_kmers()         # Start counting
    result = kmer_counter.get_tally()  # Get the result
    """

    def __init__(self, tax_id_set, clade_roots2children_map, children2clade_roots_map):

        self.tax_id_set = tax_id_set
        self.has_tally = False

        # Dictionary mappings of clade root taxIDs to children,
        # and clade children to clade root taxIDs
        self.clade_roots2children_map = clade_roots2children_map
        self.children2clade_roots_map = children2clade_roots_map

        # Main data structure to hold the kmer counts for each taxon.
        # Using the Counter class from collections. When updating it, it will
        # automatically add together the values of the matching keys.
        self.kmer_tally = {
            tax_id: Counter() for tax_id in self.tax_id_set}

    def tally_kmers(self, classifications=None, read_pickle=None, read_pickle_list=None, output_pickle=None):
        """
        Function to count kmers and return a final result.

        Input one of
            (1) a kraken2 classifications filename
            (2) a pickle of kmer tallies created with this function (see below)
            (3) a file containing a list of pickle filenames (one per line)

        Optionally takes a filename to pickle the kmer counts
        from the given classifications file to.
        """

        if self.has_tally:
            logger.warning('Will overwrite current kmer tally.')

        # Depending on type of input, get kmer tally
        if classifications:
            # From kraken 2 classifications
            classifications = vet_input_path(classifications)
            vet_kraken2_file(classifications)
            logger.info('Tallying kmers from {}...'.format(
                classifications.name))
            self.tally_per_clade(classifications)  # Count kmers

        elif read_pickle:
            # From a kmer_tally saved in a pickle
            logger.info('Reading already tallied kmers from {}...'.format(
                read_pickle))
            read_pickle = vet_input_path(read_pickle)
            self.kmer_tally = pickle.load(open(read_pickle, 'rb'))

        elif read_pickle_list:
            # From multiple pickles, add together to one kmer_tally
            read_pickle_list = vet_input_path(read_pickle_list)
            logger.info('Will combine kmer counts from pickles in {}.'.format(
                read_pickle_list))
            self.kmer_tally = self.combine_pickles(read_pickle_list)

        else:
            msg = 'You need to input one of [\'classifications\', ' + \
                  '\'read_pickle\', \'read_pickle_list\'].'
            raise KmerCounterException(msg)

        self.has_tally = True

        if output_pickle:
            output_pickle = vet_output_path(output_pickle)
            logger.info('Saving the kmer tally datastructure in {}...'.format(
                output_pickle))
            pickle.dump(self.kmer_tally, open(output_pickle, 'wb'))

    def combine_pickles(self, input_file_list):
        """Add kmer counts from multiple pickles."""

        pickles_list = []
        with read_file(input_file_list) as f:
            for line in f:
                pickles_list.append(line.strip())

        kmer_tally = {tax_id: Counter() for tax_id in self.tax_id_set}

        for i, pickle_file in enumerate(pickles_list):
            pickle_data = self.load_pickle(pickle_file)

            for tax_id, pickled_tally in pickle_data.items():
                kmer_tally[tax_id].update(pickled_tally)

            logger.info('Loaded kmer tally from {} ({}/{})...'.format(
                pickle_file, i+1, len(pickles_list)))

        return kmer_tally

    def load_pickle(self, pickle_file):
        """Load a pickle and return it. Check that the tax_ids are the same."""

        def taxa_error():
            msg = 'Pickle files appear to have been created with ' + \
                  'different tax_ids as clade definitions. Make ' + \
                  'sure to always use the same tax_ids when ' + \
                  'tallying kmers and calculating metrics.'
            raise KmerCounterException(msg)

        pickle_data = pickle.load(open(pickle_file, 'rb'))

        # The taxa in the pickle should be the same as in self.tax_id_set
        same_taxa = pickle_data.keys() == self.tax_id_set
        if not same_taxa:
            taxa_error()

        return pickle_data

    def get_tally(self):
        """
        Getter for self.kmer_tally.
        """

        if not self.has_tally:
            msg = 'Have not counted kmers yet. Run tally_kmers first.'
            raise KmerCounterException(msg)

        return self.kmer_tally

    def tally_per_clade(self, input_file):
        """
        Parse a classifications file. For all clades rooted at the tax_ids in
        self.tax_id_set, find all reads that are classified there and count
        the kmers and where they hit.

        Saves info in the main datastructure (self.kmer_tally), which is a
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
                    self.kmer_tally[tax_id_root].update(kmer_dict)

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

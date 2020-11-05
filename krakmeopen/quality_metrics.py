#!/usr/bin/env python3

import logging
import pathlib
import pandas as pd
from tally_kmers import KmerCounter
from utility_functions import read_file
from stringmeup.taxonomy import TaxonomyTree

# Set up logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
format = '%(asctime)s : %(name)-12s: %(levelname)-6s    %(message)s'
date_format = '%Y-%m-%d [%H:%M:%S]'
formatter = logging.Formatter(fmt=format, datefmt=date_format)
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)
logger.propagate = False


class ClassificationMetricsFaultyArguments(Exception):
    """Raised when a wrong or missing argument is given at instantiation."""
    pass


class ClassificationMetricsNotImplementedError(Exception):
    """Raised when trying to use functionality that is not implemented yet."""
    pass


class MetricsTabulator:
    """
    Calculates a number of metrics based on the kmer classifications from the
    reads classified to the clades rooted at the supplied tax IDs.

    The following metrics are calculated:
    1  nkmers_total                            Total number of kmers
    2  nkmers_classified                       Total number of classified kmers
    3  nkmers_unclassified                     Total number of unclassified kmers
    4  nkmers_clade                            Total number of kmers classified to any tax ID within the clade
    5  nkmers_lineage                          Total number of kmers classified to any tax ID directly above the clade root tax ID
    6  confidence_origin                       The confidence score for the clade, calculated as described by Kraken2
    7  confidence_classified                   An alternative confidence score where the unclassified kmers are removed from the denominator
    8  other_kmers_lineage_ratio               Ratio of nkmers_lineage / (nkmers_total - nkmers_clade)
    9  other_kmers_root_ratio                  Ratio of "kmers classified to root" / (nkmers_total - nkmers_clade)
    10 other_kmers_classified_ratio            Ratio of (nkmers_total - nkmers_clade - nkmers_unclassified) / (nkmers_total - nkmers_clade)
    11 other_kmers_distance                    Average distance between the clade root tax ID and the tax IDs which kmers are classified to
    12 other_kmers_distance_lineage_excluded   Like other_kmers_distance but kmers classified to tax IDs above the clade are excluded

    Usage:
    metrics_tabulator = MetricsTabulator(
        classifications_file,
        tax_id_file | tax_id,
        names,
        nodes)  # Initialize
    metrics_df = metrics_tabulator.tabulate_metrics()  # Get metrics

    metrics_df: a pandas.DataFrame object. Tax IDs are rows and metrics are
                columns.
    """

    def __init__(self,
                 classifications_file,
                 names,
                 nodes,
                 tax_id_file=None,
                 tax_id=None):

        self.tax_id_set = self.get_tax_ids(tax_id_file, tax_id)
        self.classifications = self.vet_kraken2_input(classifications_file)
        self.taxonomy_tree = self.get_taxonomy(names, nodes)

        # Dictionary mappings of clade root taxIDs to children,
        # and clade children to clade root taxIDs
        self.clade_roots2children_map, \
        self.children2clade_roots_map = self.clade_mappings()

    def get_taxonomy(self, names, nodes):
        """Build a taxonomy tree from names and nodes"""

        taxonomy_tree = TaxonomyTree(
            names_filename=names,
            nodes_filename=nodes)

        return taxonomy_tree

    def make_data_frame(self, a_dict):
        """Summarize a dict into a pandas.DataFrame."""

        df = pd.DataFrame.from_dict(a_dict, orient='index')
        df = df.reset_index()
        df = df.rename(columns={'index': 'tax_id'})

        return df

    def get_kmer_tally(self):
        """Returns a pandas.DataFrame object representing the self.kmer_tally"""

        kmer_tally_df = self.make_data_frame(self.kmer_tally)

        return kmer_tally_df

    def tabulate_metrics(self, tally_only=False, pickle_output=None, pickle_input=None):
        """Main loop to calculate the metrics."""

        self.kmer_tally = self.tally_kmers(pickle_input=pickle_input)

        # Tally the kmers and save the result in a pickle (Not Implemented)
        if tally_only:
            msg = "Tallying kmers without calculating metrics is not implemented yet."
            raise ClassificationMetricsNotImplementedError(msg)

            if not pickle_output:
                msg = "Must provide an output filename for the pickle."
                raise ClassificationMetricsFaultyArguments(msg)

            pickle_output = self.vet_output_path(pickle_output)

            # pickle output
            logger.info('Saving the kmer tally datastructure in {}'.format(
                pickle_output))

            return

        logger.info('Tabulating classification metrics per clade...')

        # Set up stuff for progress reporting
        n_taxa = len(self.tax_id_set)
        tenth = max(round(n_taxa / 10), 1)  # Report every 10th percent
        report_points = set(range(tenth, n_taxa+1, tenth))

        # Calculate metrics for each clade and store them in metrics_dict
        metrics_dict = {clade_id: None for clade_id in self.tax_id_set}
        for i, clade_id in enumerate(self.tax_id_set):
            clade_metrics_dict = self.tabulate_clade_metrics(clade_id)

            # Save to datastructure
            metrics_dict[clade_id] = clade_metrics_dict

            # Progress reporting
            if (i + 1) in report_points:
                logger.info('Processed {}/{} clades ({}%)...'.format(
                    i+1, n_taxa, round(i+1/n_taxa*100)))

        logger.info('Finished tabulating metrics.')

        logger.info('Converting result to table...')
        metrics_df = self.make_data_frame(metrics_dict)
        logger.info('Done.')

        return metrics_df

    def tally_kmers(self, pickle_input=None):
        """
        Get the kmer tallies from a KmerCounter object, or read a pickle.
        """

        # If kmer tallies are being read from a pickle (Not Implemented)
        if pickle_input:
            logger.info('Reading already tallied kmers from {}'.format(pickle_input))
            result = self.read_tally_pickle(pickle_input)
            logger.info('Done.')

        else:
            logger.info('Tallying kmers...')
            kmer_counter = KmerCounter(
                self.tax_id_set,
                self.clade_roots2children_map,
                self.children2clade_roots_map)
            kmer_counter.tally_kmers(self.classifications)
            result = kmer_counter.get_tally()
            logger.info('Done tallying.')

        return result

    def clade_mappings(self):
        """
        Create a dictionary mapping of all children (and root nodes) in all clades
        rooted at the tax_ids in the tax_id_set. The keys are the tax_ids in the
        clades, the values are the clade root node's tax_id.

        Create another dictionary where the keys are the clade tax_ids from
        tax_id_set, and the values are sets containing the tax_ids for the members
        of the specific clade.
        """

        logger.info('Creating maps of clades rooted at supplied tax IDs...')

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

    def tabulate_clade_metrics(self, clade_id):
        """
        Function to calculate a number of metrics for the supplied clade_id.
        Uses only the number of kmers that hit throughout the taxonomy,
        not reads.

        kmer_tally: Result returned from tally_kmers.KmerCounter.get_tally()
                    A dict <clade_id: collections.Counter>, where the Counter
                    objects have tax_ids as keys and kmer counts as values.
        """

        # A dictionary that will be populated with the different metrics
        metrics_dict = {
            'nkmers_total': None,
            'nkmers_classified': None,
            'nkmers_unclassified': None,
            'nkmers_clade': None,
            'nkmers_lineage': None,
            'confidence_original': None,
            'confidence_classified': None,
            'other_kmers_lineage_ratio': None,
            'other_kmers_root_ratio': None,
            'other_kmers_classified_ratio': None,
            'other_kmers_distance': None,
            'other_kmers_distance_lineage_excluded': None}

        # Get the relevant kmers
        clade_kmers_tally = self.kmer_tally[clade_id]

        # If the tax_id has no reads classified to it, it will have no kmers to
        # operate with. Return.
        if len(clade_kmers_tally) == 0:
            return metrics_dict

        # Get the tax_ids for the members of the clade rooted at clade_id
        clade_members = self.clade_roots2children_map[clade_id]

        # The lineage of the clade, make a set of it
        lineage = set(self.taxonomy_tree.get_lineage([clade_id])[clade_id])
        lineage.remove(clade_id) # Only want tax_ids above the clade tax_id

        # Dictionary comprehension and set intersection to extract the
        # <tax_id>: <num_kmers> for members of the clade
        clade_kmers = {
            key: clade_kmers_tally[key]
            for key in clade_members & clade_kmers_tally.keys()}

        # Dictionary comprehension and set complement to extract the
        # <tax_id>: <num_kmers> for non-members of the clade
        other_kmers = {
            key: clade_kmers_tally[key]
            for key in clade_kmers_tally.keys() - clade_members}

        # Dictionary comprehension and set intersection to extract the
        # <tax_id>: <num_kmers> for tax_ids in the lineage
        lineage_kmers = {
            key: clade_kmers_tally[key]
            for key in lineage & clade_kmers_tally.keys()}

        # Sum the number of kmers for each class of kmer
        num_clade_kmers = sum(clade_kmers.values())
        num_other_kmers = sum(other_kmers.values())
        num_lineage_kmers = sum(lineage_kmers.values())
        total_kmers = num_clade_kmers + num_other_kmers

        # The ratio of kmers hitting the root node over the number of kmers
        # not hitting within the clade
        if 1 in clade_kmers_tally:
            num_root_kmers = clade_kmers_tally[1]
            other_kmers_root_ratio = num_root_kmers / num_other_kmers
        else:
            other_kmers_root_ratio = 0.0

        # The number of kmers that have been classified to a tax_id, but not
        # within the clade
        if 0 in other_kmers:
            total_kmers_unclassified = other_kmers[0]
            num_other_kmers_classified = num_other_kmers - total_kmers_unclassified
        else:
            total_kmers_unclassified = 0
            num_other_kmers_classified = num_other_kmers

        # Total classified kmers
        total_kmers_classified = num_other_kmers_classified + num_clade_kmers

        # The two different confidence scores
        confidence_original = num_clade_kmers / total_kmers
        confidence_classified = num_clade_kmers / total_kmers_classified

        # Kmers hitting the lineage / kmers not hitting the clade members
        other_kmers_lineage_ratio = num_lineage_kmers / num_other_kmers if num_other_kmers > 0 else None

        # For the kmers not hitting in the clade, the ratio of kmers that hit a
        # tax_id (classified) over all other kmers not witin the clade (including non-classified)
        other_kmers_classified_ratio = num_other_kmers_classified / num_other_kmers if num_other_kmers > 0 else None

        # Getting distance metrics
        distance_dict = {tax_id: None for tax_id in other_kmers if tax_id != 0}  # We can't calculate a distance to 'Unclassified' (tax_id = 0)
        for tax_id in distance_dict:

            # The tax_id is an ancestor of the clade root
            if tax_id in lineage:

                # The distance between the clade root and the ancestor
                distance = self.taxonomy_tree.get_distance(tax_id, clade_id)

            # Not an ancestor, must compute two distances and add them together
            else:

                # Get the lineage of the tax_id
                tax_id_lineage = self.taxonomy_tree.get_lineage([tax_id])[tax_id]
                tax_id_lineage.reverse()  # Flip the lineage so that it goes from leaf to root

                # Loop to find the lowest common ancestor (lca) of the clade id and
                # the tax_id that we are currently getting the distance to
                lca = None
                for ancestor in tax_id_lineage:
                    if ancestor in lineage:
                        lca = ancestor
                        break

                # Find the distance between the lca and the clade id
                clade_lca_distance = self.taxonomy_tree.get_distance(ancestor, clade_id)

                # Find the distance between the lca and the tax_id currently
                # being investigated
                tax_id_lca_distance = self.taxonomy_tree.get_distance(ancestor, tax_id)

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
                num_other_kmers_classified - num_lineage_kmers) \
                if (num_other_kmers_classified - num_lineage_kmers) > 0 \
                else None  # If no other_kmers have hit outside of lineage, we want to avoid dividing by 0

        # Save the values of all variables
        metrics_dict = {
            'nkmers_total': total_kmers,
            'nkmers_classified': total_kmers_classified,
            'nkmers_unclassified': total_kmers_unclassified,
            'nkmers_clade': num_clade_kmers,
            'nkmers_lineage': num_lineage_kmers,
            'confidence_original': confidence_original,
            'confidence_classified': confidence_classified,
            'other_kmers_lineage_ratio': other_kmers_lineage_ratio,
            'other_kmers_root_ratio': other_kmers_root_ratio,
            'other_kmers_classified_ratio': other_kmers_classified_ratio,
            'other_kmers_distance': other_kmers_distance,
            'other_kmers_distance_lineage_excluded': other_kmers_distance_lineage_excluded}

        # Done
        return metrics_dict

    def get_tax_ids(self, tax_id_file, tax_id):
        """
        Returns a set containing ints of taxonomic IDs. Reads either from a file
        with potentially multiple taxonomic IDs, or a single taxonomic ID from
        the command line.
        """
        tax_id_set = set()

        # Loop over each line in the supplied file and add the taxonomic IDs to
        # the list
        if tax_id_file:
            logger.info('Reading taxonomic IDs from {}...'.format(tax_id_file))
            with read_file(tax_id_file) as f:
                for l in f:
                    l = int(l.strip())
                    tax_id_set.add(l)
                logger.info('Found {} taxonomic IDs.'.format(len(tax_id_set)))

        # If no file, taxonomic ID must come from the command line
        else:
            logger.info('Using a single taxonomic ID from the command line: {}.'.format(tax_id))
            tax_id_set.add(tax_id)

        return tax_id_set

    def vet_kraken2_input(self, filename):
        """
        Check if the given input file exists, and if it appears to be of
        correct format (a kraken2 classifications file).

        Returns a pathlib.Path object if everything is OK, raises
        KmerCounterFaultyArguments exception if not.
        """

        putative_path = pathlib.Path(filename)

        if putative_path.exists():
            if putative_path.is_file():
                self.conforms_to_kraken2_format(putative_path)
            else:
                msg = ('The given kraken2 classifications file is not ' + \
                'infact a file. You input {}.'.format(putative_path))
                raise ClassificationMetricsFaultyArguments(msg)

        else:
            msg = ('Could not find the specified kraken2 classifications ' + \
            'file. You input {}.'.format(putative_path))
            raise ClassificationMetricsFaultyArguments(msg)

        return putative_path

    def conforms_to_kraken2_format(self, kraken2_file):
        """
        Raises an error if the file doesn't conform to certain properties
        of a Kraken2 classifications file.
        """

        logger.debug('Validating input classifications file...')
        with read_file(kraken2_file) as f:
            line = f.readline()
            line_elements = line.strip().split('\t')

            # Number of columns as expected?
            # 5 columns in a file from Kraken2
            # 6 columns when using the fork of Kraken2, see
            # https://github.com/danisven/StringMeUp
            num_cols = len(line_elements) in [5, 6]

            # Line must start with C or U (as in Classified/unclassified)
            line_start = line_elements[0] in ['U', 'C']

            # Is the data paired or not?
            paired_input = '|' in line_elements[3]

            if paired_input:
                kmer_string_ok = len(line_elements[-1].split('|:|')) == 2
                read_len_ok = len(line_elements[3].split('|')) == 2
            else:
                # Last column should contain colons between kmer/taxon pairs
                kmer_string_ok = ":" in line_elements[-1]

                # Read length column should be convertable to int
                try:
                    int(line_elements[3])
                except ValueError:
                    read_len_ok = False
                else:
                    read_len_ok = True

            if num_cols and line_start and read_len_ok and kmer_string_ok:
                logger.debug('Validation OK.')
                return
            else:
                msg = ("The given file ({}) doesn't appear to be a" + \
                       "kraken2 classifications file.".format(kraken2_file))
                logger.debug('First line of input: {}'.format(line))
                logger.debug('num_cols: {}'.format(num_cols))
                logger.debug('line_start: {}'.format(line_start))
                logger.debug('read_len_col: {}'.format(read_len_ok))
                logger.debug('kmer_string: {}'.format(kmer_string_ok))
                raise ClassificationMetricsFaultyArguments(msg)

    def read_tally_pickle(self, pickle_input):
        """
        Read a pickle file created with count_kmers.py.
        """
        msg = "Loading a kmer tally pickle is not yet implemented."
        raise ClassificationMetricsNotImplementedError(msg)

        kmer_tally = pickle.load(pickle_input, 'rb')
        return kmer_tally

    def vet_output_path(self, filename):
        """
        Check if the given output file name is OK. Return the path as a
        pathlib.Path object if it is, else will raise
        KmerCounterFaultyArguments exception.
        """

        putative_path = pathlib.Path(filename)

        if putative_path.exists():
            if putative_path.is_dir():
                msg = ("The given output filename ({}) is a directory. " +
                        "You need to specify a file.".format(putative_path))
                raise ClassificationMetricsFaultyArguments(msg)

            else:
                msg = ("The given output filename ({}) already exists. " +
                       "I won't overwrite it. Remove it or specify another " +
                       "output file, then retry.".format(putative_path))
                raise ClassificationMetricsFaultyArguments(msg)

        else:
            if not putative_path.parent.is_dir():
                msg = ("The given output filename ({}) suggests writing to " +
                      "a directory that doesn't exist.".format(putative_path))
                raise ClassificationMetricsFaultyArguments(msg)
            else:
                return putative_path


if __name__ == '__main__':
    names = '/home/daniel/work/development/StringMeUp/data/names.dmp'
    nodes = '/home/daniel/work/development/StringMeUp/data/nodes.dmp'
    tax_id = 9606
    tax_id_file = '/home/daniel/work/development/scratch/KrakMeOpen/tax_id_file.txt'
    # classifications_file = '/home/daniel/work/development/StringMeUp/data/Ki-2014-27-74_sample.kraken2'
    classifications_file = '/home/daniel/work/development/scratch/KrakMeOpen/Ki-2012-33-500_sample100k.kraken2'
    metrics_tabulator = MetricsTabulator(
        classifications_file=classifications_file,
        tax_id_file=tax_id_file,
        names=names,
        nodes=nodes)
    metrics_df = metrics_tabulator.tabulate_metrics()
    print(metrics_df)

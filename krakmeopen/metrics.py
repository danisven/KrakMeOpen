#!/usr/bin/env python3

import pickle
import logging
import pathlib
import pandas as pd
from collections import Counter
from krakmeopen.kmers import KmerCounter
from krakmeopen.utilities import read_file
from stringmeup.taxonomy import TaxonomyTree

logger = logging.getLogger(__name__)


class MetricsTabulatorException(Exception):
    """Raised when an error occurs in MetricsTabulator."""
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

    # Tally kmers and calculate metrics for a kraken 2 classifications file
    metrics_tabulator = MetricsTabulator(
        input_classifications,
        tax_id_file | tax_id,
        names,
        nodes)  # Initialize
    metrics_df = metrics_tabulator.tabulate_metrics()  # Get metrics

    metrics_df: a pandas.DataFrame object. Tax IDs are rows and metrics are
                columns.
    """

    def __init__(self,
                 names,
                 nodes,
                 tax_id_file=None,
                 tax_id=None):

        self.tax_id_set = self.get_tax_ids(tax_id_file, tax_id)
        self.taxonomy_tree = self.get_taxonomy(names, nodes)
        self.kmer_tally = None

        # Dictionary mappings of clade root taxIDs to children,
        # and clade children to clade root taxIDs
        self.clade_roots2children_map, \
        self.children2clade_roots_map = self.clade_mappings()

    def get_taxonomy(self, names, nodes):
        """Build a taxonomy tree from names and nodes"""

        logger.info('Creating taxonomy...')
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

        if self.kmer_tally:
            kmer_tally_df = self.make_data_frame(self.kmer_tally)
            return kmer_tally_df

        else:
            msg = 'There is no kmer tally available. Run tabulate_metrics() first.'
            raise MetricsTabulatorException(msg)

    def tabulate_metrics(self, input_classifications=None, input_pickle=None, input_file_list=None, tally_only=False):
        """
        Main loop to calculate the metrics from a kmer tally.

        Kmer tally can be created from
            (1) self.classifications (kraken 2 output, set at instantiation)
            (2) input_pickle (pickle file, created with krakmeopen --output_pickle)
            (3) input_file_list (file with one pickle file name per line)

        Give tally_only=FILE to output only a pickle file and no metrics.
        """

        inputs = [input_classifications, input_pickle, input_file_list]

        if input_pickle and input_file_list:
            msg = 'Cannot input a pickle file and a list of files at the same time.'
            raise MetricsTabulatorException(msg)

        elif all(f is None for f in inputs):
            msg = 'Must supply an input for MetricsTabulator.tabulate_metrics()'
            raise MetricsTabulatorException(msg)

        kmer_counter = KmerCounter(
            self.tax_id_set,
            self.clade_roots2children_map,
            self.children2clade_roots_map)
        kmer_counter.tally_kmers(
            classifications=input_classifications,
            read_pickle=input_pickle,
            read_pickle_list=input_file_list,
            output_pickle=tally_only)
        self.kmer_tally = kmer_counter.get_tally()

        # If only tallying kmers (no metrics) and saving as pickle
        if tally_only:
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
                    i+1, n_taxa, round((i+1)/n_taxa*100)))

        logger.info('Converting result to table...')
        metrics_df = self.make_data_frame(metrics_dict)

        return metrics_df

    def clade_mappings(self):
        """
        Create a dictionary mapping of all children (and root nodes) in all clades
        rooted at the tax_ids in the tax_id_set. The keys are the tax_ids in the
        clades, the values are the clade root node's tax_id.

        Create another dictionary where the keys are the clade tax_ids from
        tax_id_set, and the values are sets containing the tax_ids for the members
        of the specific clade.
        """

        logger.info('Creating clade maps...')

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

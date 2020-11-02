#!/usr/bin/env python3

import logging
import pathlib
import pickle
from utility_functions import read_file
from stringmeup.taxonomy import TaxonomyTree

class KmerCounterFaultyArguments(Exception):
    """Raised when a wrong or missing argument is given at instantiation."""
    pass


class KmerCounter:

    def __init__(self,
                 kraken2_classifications,
                 tax_id_list,
                 logger=None,
                 names=None,
                 nodes=None,
                 taxonomy_tree=None,
                 output_pickle_fn=None):

        self.logger = logger or logging.get_logger(__name__)
        self.output_pickle = self.vet_output(output_pickle_fn) if output_pickle_fn else False
        self.taxonomy_tree = self.get_taxonomy(taxonomy_tree, names, nodes)


    def vet_output(self, filename):
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
                raise KmerCounterFaultyArguments(msg)

            else:
                msg = ("The given output filename ({}) already exists. " +
                       "I won't overwrite it. Remove it or specify another " +
                       "output file, then retry.".format(putative_path))
                raise KmerCounterFaultyArguments(msg)

        else:
            if not putative_path.parent.is_dir():
                msg = ("The given output filename ({}) suggests writing to " +
                      "a directory that doesn't exist.".format(putative_path))
                raise KmerCounterFaultyArguments(msg)
            else:
                self.logger("Will pickle the kmer counts and save to {}.".format(putative_path))

        return filename

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
                raise KmerCounterFaultyArguments(msg)

            else:
                self.logger.info("Building taxonomy tree.")
                taxonomy_tree = TaxonomyTree(
                    names_filename=names_input,
                    nodes_filename=nodes_input,
                )

        elif isinstance(taxonomy_tree_input, TaxonomyTree):
            taxonomy_tree = taxonomy_tree_input

        else:
            msg = ("The class of the supplied taxonomy_tree argument was not TaxonomyTree, " +
                   "instead it was {}.".format(type(taxonomy_tree_input)))
            raise KmerCounterFaultyArguments(msg)

        return taxonomy_tree

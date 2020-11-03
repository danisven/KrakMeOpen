#!/usr/bin/env python3

import logging
from tally_kmers import KmerCounter

# Set up logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
format = '%(asctime)s: %(name)s: %(levelname)s:    %(message)s'
date_format = '%Y-%m-%d [%H:%M:%S]'
formatter = logging.Formatter(fmt=format, datefmt=date_format)
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)

class ClassificationMetrics:
    def __init__(self, names, nodes):
        pass
        # Most of what is in classification_quality_metrics.py will go here.

    def collect_tallies():
        """
        Get the kmer tallies... Call on KmerCounter.
        """
        # Getting the kmer tallies
        kmer_counter = KmerCounter(tax_ids_set, taxonomy_tree)
        kmer_counter.tally_kmers(classifications_file)
        result = kmer_counter.get_tally()


# Check if putative output pickle filename is OK to write to.
# output_pickle: pathlib.Path object.
output_pickle = self.vet_output(output_pickle_fn) \
    if output_pickle_fn else False

# Conversely, check if input is OK.
# classifications: pathlib.Path object.
classifications = self.vet_input(kraken2_classifications)

def vet_input(self, filename):
    """
    Check if the given input file exists, and if it appears to be of
    correct format (a kraken2 classifications file).

    Returns a pathlib.Path object if everything is OK, raises
    KmerCounterFaultyArguments exception if not.
    """

    def conforms_to_kraken2_format(kraken2_file):
        """
        Raises an error if the file doesn't conform to certain properties
        of a Kraken2 classifications file.
        """


    putative_path = pathlib.Path(filename)

    if putative_path.exists():
        if putative_path.is_file():
            conforms_to_kraken2_format(putative_path)
        else:
            msg = ('The given kraken2 classifications file is not ' + \
                   'infact a file. You input {}.'.format(putative_path))
            raise KmerCounterFaultyArguments(msg)

    else:
        msg = ('Could not find the specified kraken2 classifications ' + \
               'file. You input {}.'.format(putative_path))
        raise KmerCounterFaultyArguments(msg)

    return putative_path

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
            logger("Will pickle the kmer counts and save to {}.".format(putative_path))

    return putative_path

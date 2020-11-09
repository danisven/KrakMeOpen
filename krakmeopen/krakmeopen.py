#!/usr/bin/env python3

__version__ = "0.1.0"

import yaml
import pickle
import logging
import logging.config
import argparse
from krakmeopen import conf
from krakmeopen.log_formatter import CustomFormatter
from krakmeopen.quality_metrics import MetricsTabulator
import importlib.resources as pkg_resources

# TODO: requires pyyaml

def get_arguments():
    """
    Wrapper function to get the command line arguments. Inserting this piece of code
    into its own function for conda compatibility.
    """

    parser = argparse.ArgumentParser(
        prog='KrakMeOpen',
        usage='krakmeopen [--input FILE | --input_pickle FILE | --input_file_list FILE] [--output FILE | --output_pickle FILE] --names FILE --nodes FILE [--tax_id INT | --tax_id_file FILE] --kmer_tally_table FILE',
        description='''
            A Kraken2 downstream analysis toolkit. More specifically, calculate
            a series of quality metrics for Kraken2 classifications.''',
        epilog='''
            The metrics are calculated on the clade-level. All kmers
            from all reads that are classified to any of the nodes in the
            clades rooted at the supplied tax IDs are aggregated, and metrics
            are calculated on those aggregations.

            Input is Kraken2 read-by-read classification files
            (can be gzipped).

            Output is a tab separated file containing the metrics.''')

    # Input arguments
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '--input',
        metavar='FILE',
        type=str,
        help='Kraken2 read-by-read classifications file.')
    input_group.add_argument(
        '--input_pickle',
        metavar='FILE',
        type=str,
        help='A pickle file containing kmer tallies, produced with --output_pickle')
    input_group.add_argument(
        '--input_file_list',
        metavar='FILE',
        type=str,
        help='''A file containing file paths to multiple pickles, one per line.
             Will calculate metrics on the sum of kmer counts from all pickles.''')

    # Output arguments
    output_group = parser.add_mutually_exclusive_group(required=True)
    output_group.add_argument(
        '--output',
        metavar='FILE',
        type=str,
        help='The file to write the quality metrics output to.')
    output_group.add_argument(
        '--output_pickle',
        metavar='FILE',
        type=str,
        help='''The pickle file to write kmer tallies to. Use this argument
             to supress calculation of quality metrics and only output kmer
             counts to a pickled file. Input the pickled file using
             --input_pickle.''')
    parser.add_argument(
        '--kmer_tally_table',
        metavar='FILE',
        required=False,
        help='File to output the complete kmer tally table for each tax ID to. Optional.')

    # The taxonomy
    parser.add_argument(
        '--names',
        metavar='FILE',
        required=True,
        help='NCBI style taxonomy names dump file (names.dmp). Required.')
    parser.add_argument(
        '--nodes',
        metavar='FILE',
        required=True,
        help='NCBI style taxonomy nodes dump file (nodes.dmp). Required.')

    # Supply relevant taxonomic ID on command line, or one or multiple taxonomic IDs
    # through a text file.
    tax_id_group = parser.add_mutually_exclusive_group(required=True)
    tax_id_group.add_argument(
        '--tax_id',
        metavar='INT',
        type=int,
        help='A taxonomic ID for a clade that you wish to calculate quality metrics for.')
    tax_id_group.add_argument(
        '--tax_id_file',
        metavar='FILE',
        type=str,
        help='''Supply multiple taxonomic IDs at once. A textfile with one
            taxonomic ID per line. Calculate quality metrics for the clades
            rooted at the taxonomic IDs in the file.''')

    return parser.parse_args()


def tally_only(metrics_tabulator, args, logger):
    """
    Tally kmers and pickle result. The output file can be loaded later and
    be used to calculate metrics at a later stage. Metrics can also be
    calculated on a combination of multiple pickles.
    """

    # Check output path and count kmers
    pickle_output = metrics_tabulator.vet_output_path(args.output_pickle)
    kmer_tally = metrics_tabulator.tally_kmers()

    logger.info('Saving the kmer tally datastructure in {} ...'.format(
        pickle_output))

    # Pickle
    pickle.dump(kmer_tally, open(pickle_output, 'wb'))

    logger.info('Saved.')
    logger.info('You can calculate quality metrics from this file at a ' + \
                'later stage by calling krakmeopen with the --input_pickle ' + \
                'or --input_file_list argument.')


def setup_logging(log_level=logging.DEBUG):
    conf_content = pkg_resources.read_text(conf, 'log_config.yaml')
    log_config = yaml.safe_load(conf_content)
    logging.config.dictConfig(log_config)
    logger = logging.getLogger(__name__)
    return logger


# Main func
def krakmeopen():
    args = get_arguments()
    logger = setup_logging()

    metrics_tabulator = MetricsTabulator(
        names = args.names,
        nodes = args.nodes,
        input_classifications = args.input,
        tax_id = args.tax_id,
        tax_id_file = args.tax_id_file)

    if args.output_pickle:
        logger.info('Will not calculate quality metrics.')
        tally_only(metrics_tabulator, args, logger)

    else:

        # Tabulate metrics and output to file
        metrics_df = metrics_tabulator.tabulate_metrics(
            input_pickle=args.input_pickle,
            input_file_list=args.input_file_list)
        logger.info('Saving quality metrics to {}.'.format(args.output))
        metrics_df.to_csv(args.output, sep='\t', index=False)

    # If user wants a kmer tally table to be output to file (pd.DataFrame)
    if args.kmer_tally_table:
        logger.info('Saving human raedable kmer tallies in {}'.format(
            args.kmer_tally_table))
        kmer_tally_df = metrics_tabulator.get_kmer_tally()
        kmer_tally_df.to_csv(args.kmer_tally_table, sep='\t', index=False)

    logger.info('Finished. Exiting KrakMeOpen.')

if __name__ == '__main__':
    krakmeopen()

#!/usr/bin/env python3

__version__ = "0.1.0"

import argparse
from krakmeopen.quality_metrics import MetricsTabulator


def get_arguments():
    """
    Wrapper function to get the command line arguments. Inserting this piece of code
    into its own function for conda compatibility.
    """

    parser = argparse.ArgumentParser(
        prog='KrakMeOpen',
        usage='krakmeopen --input FILE --output FILE --names FILE --nodes FILE [--tax_id INT | --tax_id_file FILE]',
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

    # Input Kraken2 classifications file
    parser.add_argument(
        '--input',
        metavar='FILE',
        type=str,
        required=True,
        help='Kraken2 read-by-read classifications file.')

    # Output file
    parser.add_argument(
        '--output',
        metavar='FILE',
        type=str,
        required=True,
        help='The file to write the output to.')

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

    # Output kmer tally to file
    parser.add_argument(
        '--output_kmer_tally',
        metavar='FILE',
        required=False,
        help='File to output the complete kmer tally for each tax ID to.')

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
        help='''
            Supply multiple taxonomic IDs at once. A textfile with one
            taxonomic ID per line. Calculate quality metrics for the clades
            rooted at the taxonomic IDs in the file.''')

    return parser.parse_args()


# Main func
def krakmeopen():
    args = get_arguments()

    metrics_tabulator = MetricsTabulator(
        classifications_file = args.input,
        names = args.names,
        nodes = args.nodes,
        tax_id = args.tax_id,
        tax_id_file = args.tax_id_file)

    # Tabulate metrics and output to file
    metrics_df = metrics_tabulator.tabulate_metrics()
    metrics_df.to_csv(args.output, sep='\t', index=False)

    # If user wants the full kmer tally to be output to file
    if args.output_kmer_tally:
        kmer_tally_df = metrics_tabulator.get_kmer_tally()
        kmer_tally_df.to_csv(args.output_kmer_tally, sep='\t', index=False)


if __name__ == '__main__':
    krakmeopen()

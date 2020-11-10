#!/usr/bin/env python3

import pathlib
import gzip


class OutputFileException(Exception):
    """Raised when something is wrong with the output file."""
    pass


class InputFileException(Exception):
    """Raised when something is wrong with the input file."""
    pass


def read_file(file_path):
    """
    Wrapper to read either gzipped or ordinary text file input.
    """
    if not isinstance(file_path, pathlib.Path):
        file_path = pathlib.Path(file_path)

    if file_path.suffix == '.gz':
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

def vet_output_path(filename):
    """
    Check if the given output file name is OK. Return the path as a
    pathlib.Path object if it is, else will raise OutputFileException.
    """

    putative_path = pathlib.Path(filename)

    if putative_path.exists():
        if putative_path.is_dir():
            msg = ('A given output filename ({}) is '.format(putative_path) + \
                   'a directory. You need to specify a file.')
            raise OutputFileException(msg)

        else:
            msg = ('A given output filename ({}) '.format(putative_path) + \
                   'already exists. I won\'t overwrite it. Remove it or ' + \
                   'specify another output file, then retry.')
            raise OutputFileException(msg)

    else:
        if not putative_path.parent.is_dir():
            msg = ('A given output filename ({}) '.format(putative_path) + \
                   'suggests writing to a directory that doesn\'t exist.')
            raise OutputFileException(msg)
        else:
            return putative_path

def vet_input_path(filename):
    """
    Check if the given input file exists.

    Returns a pathlib.Path object if everything is OK, raises
    InputFileException if not.
    """

    putative_path = pathlib.Path(filename)

    if putative_path.exists():
        if not putative_path.is_file():
            msg = ('A given input file is not infact a file. ' + \
                   'You input {}.'.format(putative_path))
            raise InputFileException(msg)

    else:
        msg = ('Could not find a specified input file. You input {}.'.format(
            putative_path))
        raise InputFileException(msg)

    return putative_path

def vet_kraken2_file(kraken2_file):
    """
    Raises an error if the file doesn't conform to certain properties
    of a Kraken2 classifications file.
    """

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
            return

        else:
            msg = ('The given file ({}) doesn\'t appear to be a ' + \
                   'kraken2 classifications file.'.format(kraken2_file))
            logger.debug('First line of input: {}'.format(line))
            logger.debug('num_cols: {}'.format(num_cols))
            logger.debug('line_start: {}'.format(line_start))
            logger.debug('read_len_col: {}'.format(read_len_ok))
            logger.debug('kmer_string: {}'.format(kmer_string_ok))

            raise InputFileException(msg)

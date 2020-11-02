#! /usr/bin/env python3

def read_file(filename):
    """
    Wrapper to read either gzipped or ordinary text file input.
    """
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')

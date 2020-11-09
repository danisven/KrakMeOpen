#!/usr/bin/env python3

import pathlib
import gzip

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

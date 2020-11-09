#!/usr/bin/env python3

import logging

class CustomFormatter(logging.Formatter):

    def __init__(self, fmt=None, datefmt=None, info_format=None):
        super().__init__(fmt=fmt, datefmt=datefmt, style='%')
        self.info_format = info_format
        # self.debug_format = debug_format  # Uncomment if needed

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._style._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.INFO:
            self._style._fmt = self.info_format

        # Uncomment to add debug format
        # elif record.levelno == logging.DEBUG:
        #     self._style._fmt = self.debug_format

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._style._fmt = format_orig

        return result

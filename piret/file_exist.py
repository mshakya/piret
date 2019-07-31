#! /usr/bin/env python

import os.path


def file_exist(file_path):
    """
    Check if file exists.

    returns True if file exists, if it doesn't it sends out False
    """
    if os.path.isfile(file_path):
        return True
    else:
        return False

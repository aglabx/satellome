#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 05.06.2011
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import os


def iter_filepath_folder(folder_path, recursive=True):
    """
    Iterate through all file paths in a given folder.

    :param folder_path: The path to the folder.
    :param recursive: If True, it will search subfolders as well. Default is False.
    :return: Generator yielding file paths.
    """
    if recursive:
        for dirpath, _, filenames in os.walk(folder_path):
            for filename in filenames:
                yield os.path.join(dirpath, filename)
    else:
        for filename in os.listdir(folder_path):
            full_path = os.path.join(folder_path, filename)
            if os.path.isfile(full_path):
                yield full_path

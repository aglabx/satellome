#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 06.09.2011
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com
"""
satelome settings loader.
"""
import os
import platform
import satelome

import yaml

satelome_path = satelome.__path__._path[0]

SETTINGS_FILENAME = os.path.join(satelome_path, "settings.yaml")
NGRAM_LENGTH = 23
NGRAM_N = 100000000


def load_settings():
    """Load settings from yaml file.
    @return settings
    """
    file_name = os.path.abspath(__file__)
    settings_file = os.path.join(os.path.split(file_name)[0], SETTINGS_FILENAME)
    with open(settings_file) as fh:
        settings = yaml.load(fh, Loader=yaml.FullLoader)
    myos = platform.system()
    if myos == "Windows":
        settings["satelome"]["os"] = "WIN"
    elif myos == "Darwin":
        settings["satelome"]["os"] = "OSX"
    else:
        settings["satelome"]["os"] = "NIX"
    return settings


def save_settings(settings):
    """Save settings to yaml file.
    @param settings: satelome settings
    """
    file_name = os.path.abspath(__file__)
    settings_file = os.path.join(os.path.split(file_name)[0], SETTINGS_FILENAME)
    with open(settings_file, "w") as fh:
        yaml.dump(fh, settings)

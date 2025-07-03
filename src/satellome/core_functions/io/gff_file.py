#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 10.11.2022
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com


from satellome.core_functions.models.gff3_model import Gff3FileIO


def sc_gff3_reader(gff3_file, only_fields=None):
    """Iter over gff3 file."""
    reader = Gff3FileIO()
    for gff3_obj in reader.read_online(gff3_file, only_fields=only_fields):
        yield gff3_obj

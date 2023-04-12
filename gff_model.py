#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 28.10.2014
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com 
"""
"""


from PyExp import AbstractModel
import csv


class Gff3Model(AbstractModel):
    """ Class for gff3 data wrapping.
    """

    dumpable_attributes = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]
    
    int_attributes = [
        "start",
        "end",
    ]

    @property
    def length(self):
        return abs(self.end - self.start)

    def get_coordinates(self):
        if self.end < self.start:
            return (self.end, self.start)
        return (self.start, self.end)

    def save_original(self, line):
        self.original = line

    def as_gff3(self):
        '''
        '''
        s = []
        for attr in self.dumpable_attributes:
            if not attr == "attributes":
                s.append(str(getattr(self, attr)))
        if not hasattr(self, "attributes"):
            self.attributes = {}     
        if self.attributes:
            attr_keys = [key for key in self.attributes.keys() if not isinstance(self.attributes[key], dict)]
            dict_attr_keys = [key for key in self.attributes.keys() if isinstance(self.attributes[key], dict)]
            attr_keys.sort()
            attr = []
            for k in attr_keys:
                attr.append("%s=%s" % (k, self.attributes[k]))  
            for k in dict_attr_keys:
                data = []
                for key in self.attributes[k]:
                    data.append(f"{key}:{self.attributes[k][key]}")
                attr.append("%s=%s" % (k, ",".join(data)))  
            self.raw_features = ";".join(attr)
        s.append(self.raw_features)
        s = "\t".join(s)
        return "%s\n" % s

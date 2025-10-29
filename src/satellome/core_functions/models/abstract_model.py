#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.06.2011
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import json


class AbstractModel(object):
    """Сlass for data wrapping.


    Private methods:

    - __str__(self) used dumpable_attributes
    - set_with_dict(self, dictionary)

    Initiation. Create attributes accordong to

    public properties:


    - dumpable_attributes
    - int_attributes
    - float_attributes
    - list_attributes
    - list_attributes_types
    - other_attributes

    """

    dumpable_attributes = []
    int_attributes = []
    float_attributes = []
    list_attributes = []
    list_attributes_types = {}
    other_attributes = {}

    def __init__(self):
        """Create attributes accordong to

        - dumpable_attributes
        - int_attributes
        - float_attributes
        - list_attributes
        - list_attributes_types
        - other_attributes
        """
        for attr in self.dumpable_attributes:
            setattr(self, attr, None)
        for attr in self.int_attributes:
            setattr(self, attr, 0)
        for attr in self.float_attributes:
            setattr(self, attr, 0.0)
        for attr in self.other_attributes:
            setattr(self, attr, self.other_attributes[attr])

    def __str__(self):
        """Get string representation with fields
        defined in dumpable_attributes."""
        self.preprocess_data()
        result = []
        for attr in self.dumpable_attributes:
            data = getattr(self, attr)
            if attr in self.list_attributes:
                if data is None:
                    data = []
                data = ",".join([str(x) for x in data])
            try:
                result.append(str(data).strip())
            except UnicodeEncodeError:
                result.append(str(data).strip())
        result = "%s\n" % "\t".join(result)
        return result

    def print_human_friendly(self):
        """Print human friendly string representation with fields
        defined in dumpable_attributes."""
        self.preprocess_data()
        result = []
        largest_name_length = max([len(x) for x in self.dumpable_attributes])
        print_string = "{0:%s} => {1}" % largest_name_length
        for attr in self.dumpable_attributes:
            data = getattr(self, attr)
            if attr in self.list_attributes:
                if data is None:
                    data = []
                data = ",".join([str(x) for x in data])
            try:
                print(print_string.format(attr, data.strip()))
            except UnicodeEncodeError:
                data = str(data).strip()
                print(print_string.format(attr, data.strip()))

    def get_as_string(self, dumpable_attributes):
        """Get string representation with fields
        defined in dumpable_attributes."""
        return self.__str__()

    def set_with_dict(self, dictionary):
        """Set object with dictionaty."""
        for key, value in dictionary.items():
            key, value = self.preprocess_pair(key, value)
            try:
                if value == "None" or value is None:
                    value = None
                elif key in self.int_attributes:
                    value = int(value)
                elif key in self.float_attributes:
                    value = float(value)
                elif key in self.list_attributes:
                    if not value:
                        value = []
                        continue
                    value = value.split(",")
                    value = [self.list_attributes_types[key](x) for x in value]
                setattr(self, key, value)
            except ValueError as e:
                print(self.dumpable_attributes)
                print(dictionary.items())
                raise ValueError(e)
            except TypeError as e:
                print(self.dumpable_attributes)
                print(dictionary.items())
                raise TypeError(e)

    def set_with_list(self, data):
        """Set object with list."""
        n = len(data)
        dumpable_attributes = self.dumpable_attributes
        if n != len(self.dumpable_attributes):
            if (
                hasattr(self, "alt_dumpable_attributes")
                and len(self.alt_dumpable_attributes) == n
            ):
                dumpable_attributes = self.alt_dumpable_attributes
            else:
                print(data)
                raise Exception("Wrong number of fields in data.")
        for i, value in enumerate(data):
            key = dumpable_attributes[i]
            if value == "None":
                value = None
            elif key in self.int_attributes:
                value = int(value)
            elif key in self.float_attributes:
                value = float(value)
            elif key in self.list_attributes:
                if value:
                    value = value.split(",")
                    value = [self.list_attributes_types[key](x) for x in value]
                else:
                    value = []
            setattr(self, key, value)

    def as_dict(self):
        """ """
        return self.get_as_dict()

    def get_as_dict(self):
        """Get dictionary representation with fields
        defined in dumpable_attributes"""
        self.preprocess_data()
        result = {}
        for attr in self.dumpable_attributes:
            result[attr] = getattr(self, attr)
        return result

    def get_as_json(self, preprocess_func=None):
        """Return JSON representation."""
        self.preprocess_data()
        d = self.get_as_dict()
        if preprocess_func:
            d = preprocess_func(d)
        return json.dumps(d)

    def preprocess_data(self):
        """Any data preprocessing before returning."""
        pass

    def preprocess_pair(self, key, value):
        """Any data preprocessing before initiation from dictionary."""
        return key, value

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        return setattr(self, key, value)


class AbstractMongoModel(object):
    mongo_key = lambda x: x

    def __init__(self, mongo_db, mongo_collection):
        mongo_db = ""

    def mongo_add(self):
        pass

    def mongo_get(self, id):
        pass

    def mongo_update(self):
        pass

    def mongo_remove(self):
        pass

    def mongo_update_fields(self):
        pass

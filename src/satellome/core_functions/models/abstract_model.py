#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.06.2011
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import json
import logging

logger = logging.getLogger(__name__)


class AbstractModel(object):
    """
    Abstract base class for data model objects with automatic attribute handling.

    Provides automatic type conversion and serialization for model classes.
    Subclasses define their schema via class attributes, and this base class
    handles initialization, type conversion, and serialization/deserialization.

    Class Attributes:
        dumpable_attributes (list): List of attribute names to include in serialization
        int_attributes (list): Attributes that should be converted to integers
        float_attributes (list): Attributes that should be converted to floats
        list_attributes (list): Attributes that should be converted to lists
        list_attributes_types (dict): Mapping of list attribute names to element types
        other_attributes (dict): Additional attributes with custom handling

    Methods:
        set_with_dict: Initialize attributes from a dictionary
        get_as_dict: Serialize to dictionary
        get_as_json: Serialize to JSON string
        preprocess_data: Hook for data preprocessing before serialization
        preprocess_pair: Hook for data preprocessing during initialization

    Example:
        >>> class MyModel(AbstractModel):
        ...     dumpable_attributes = ["name", "count"]
        ...     int_attributes = ["count"]
        >>> obj = MyModel()
        >>> obj.set_with_dict({"name": "test", "count": "42"})
        >>> obj.count  # Automatically converted to int
        42
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
        """Set object with dictionary, auto-detecting legacy format if applicable."""
        # Auto-detect and convert legacy format if model supports it
        if hasattr(self, 'legacy_to_new_field_mapping'):
            # Check if this looks like legacy format by checking for legacy-only fields
            legacy_only_fields = set(self.legacy_to_new_field_mapping.keys()) - set(self.dumpable_attributes)
            if any(field in dictionary for field in legacy_only_fields):
                logger.debug(f"Detected legacy format for {self.__class__.__name__}, converting...")
                dictionary = self._convert_legacy_dict(dictionary)

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
                logger.error(
                    f"ValueError while parsing data in {self.__class__.__name__}. "
                    f"Expected schema: {self.dumpable_attributes}. "
                    f"Input data: {dict(list(dictionary.items())[:5])}...",  # Show first 5 items
                    exc_info=True
                )
                raise ValueError(
                    f"Failed to parse data for {self.__class__.__name__}. "
                    f"Check that all fields match expected types."
                ) from e
            except TypeError as e:
                logger.error(
                    f"TypeError while parsing data in {self.__class__.__name__}. "
                    f"Expected schema: {self.dumpable_attributes}. "
                    f"Input data: {dict(list(dictionary.items())[:5])}...",  # Show first 5 items
                    exc_info=True
                )
                raise TypeError(
                    f"Failed to parse data for {self.__class__.__name__}. "
                    f"Check that all fields have correct types."
                ) from e

    def _convert_legacy_dict(self, legacy_dict):
        """Convert legacy format dictionary to new format.

        Args:
            legacy_dict: Dictionary with legacy field names

        Returns:
            Dictionary with new field names
        """
        new_dict = {}
        for legacy_key, legacy_value in legacy_dict.items():
            new_key = self.legacy_to_new_field_mapping.get(legacy_key)
            if new_key is not None:  # Keep only mapped fields
                new_dict[new_key] = legacy_value
        return new_dict

    def set_with_list(self, data):
        """Set object with list, auto-detecting legacy format if applicable."""
        n = len(data)
        dumpable_attributes = self.dumpable_attributes

        # Auto-detect legacy format
        if n != len(self.dumpable_attributes):
            if (
                hasattr(self, "legacy_dumpable_attributes")
                and len(self.legacy_dumpable_attributes) == n
            ):
                logger.debug(f"Detected legacy format ({n} fields), converting to new format...")
                dumpable_attributes = self.legacy_dumpable_attributes
                # Convert to dict first, then use legacy mapping
                legacy_dict = {key: value for key, value in zip(dumpable_attributes, data)}
                new_dict = self._convert_legacy_dict(legacy_dict)
                self.set_with_dict(new_dict)
                return
            elif (
                hasattr(self, "alt_dumpable_attributes")
                and len(self.alt_dumpable_attributes) == n
            ):
                dumpable_attributes = self.alt_dumpable_attributes
            else:
                logger.error(f"Wrong number of fields in data: expected {len(self.dumpable_attributes)}, got {n}. Data: {data}")
                raise ValueError(
                        f"Expected {len(self.dumpable_attributes)} fields in {self.__class__.__name__}, "
                        f"got {n}. Sample data: {data[:50]}..."
                        )

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
        """
        Get dictionary representation of the model.

        Alias for get_as_dict(). Calls preprocess_data() before serialization
        to allow subclasses to perform any necessary data transformations.

        Returns:
            dict: Dictionary with keys from dumpable_attributes

        See Also:
            get_as_dict: The underlying implementation
        """
        return self.get_as_dict()

    def get_as_dict(self):
        """
        Get dictionary representation with fields defined in dumpable_attributes.

        Calls preprocess_data() to allow subclass-specific transformations,
        then builds a dictionary containing only the attributes listed in
        the dumpable_attributes class variable.

        Returns:
            dict: Dictionary with keys from dumpable_attributes and their current values

        Example:
            >>> class MyModel(AbstractModel):
            ...     dumpable_attributes = ["name", "value"]
            >>> obj = MyModel()
            >>> obj.name, obj.value = "test", 42
            >>> obj.get_as_dict()
            {'name': 'test', 'value': 42}
        """
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
        """
        Hook for preprocessing data before serialization.

        Override this method in subclasses to perform any necessary data
        transformations before the object is serialized to dict/JSON.
        Called automatically by get_as_dict() and get_as_json().

        Common use cases:
        - Computing derived fields
        - Formatting values for output
        - Cleaning up temporary attributes

        Note:
            Base implementation does nothing. Subclasses should override as needed.

        Example:
            >>> class MyModel(AbstractModel):
            ...     def preprocess_data(self):
            ...         # Ensure uppercase before serialization
            ...         if hasattr(self, 'name'):
            ...             self.name = self.name.upper()
        """
        pass

    def preprocess_pair(self, key, value):
        """
        Hook for preprocessing key-value pairs during initialization.

        Override this method in subclasses to transform or validate data
        before it's assigned to object attributes. Called by set_with_dict()
        for each key-value pair in the input dictionary.

        Args:
            key (str): Attribute name
            value: Attribute value (any type)

        Returns:
            tuple: (processed_key, processed_value)

        Note:
            Base implementation returns the pair unchanged.
            Subclasses can override to:
            - Rename keys
            - Transform values
            - Validate input
            - Skip unwanted attributes

        Example:
            >>> class MyModel(AbstractModel):
            ...     def preprocess_pair(self, key, value):
            ...         # Convert all keys to lowercase
            ...         return key.lower(), value
        """
        return key, value

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        return setattr(self, key, value)


#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 17.08.2026
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com
"""Turn outside text into a safe single path component.

Organism names are free text, not identifiers. NCBI strain designations
routinely contain path separators::

    Leishmania braziliensis MHOM/BR/75/M2904
    Chlorella vulgaris CCAP 1055/1
    Leishmania panamensis MHOM/GT/2001/U1103

and quotes (``'Nostoc azollae' 0708`` and other Candidatus names). Interpolated
straight into an output path, such a name stops being a file name and becomes a
directory hierarchy that nobody created — the write then fails with
``FileNotFoundError`` at the very end of the pipeline, after every data file has
already been produced.

So any value that reaches a path from outside the pipeline goes through
:func:`safe_filename_component` first. Callers are expected to *report* when the
name had to be changed (the sanitizer itself cannot: a logged-only rename is
invisible), which is why it is a pure function returning the new name.
"""

import re

# What survives untouched: ASCII letters, digits, and the few punctuation marks
# that every filesystem, shell and downstream script handles without quoting.
_UNSAFE_RUN = re.compile(r"[^A-Za-z0-9._+-]+")
_LEADING_JUNK = re.compile(r"^[._-]+")
_TRAILING_JUNK = re.compile(r"[._-]+$")

DEFAULT_NAME = "unnamed"

# Component length cap. ext4/XFS allow 255 bytes per component, and the taxon
# name is only a prefix — ".karyo.gaps.1000bp.enhanced.svg" and friends are
# appended to it — so leave generous room for the suffixes.
MAX_COMPONENT_LENGTH = 120


def safe_filename_component(name, default=DEFAULT_NAME, max_length=MAX_COMPONENT_LENGTH):
    """Return *name* reduced to one safe path component.

    Runs of unsafe characters (path separators included) collapse into a single
    underscore; leading/trailing dots, dashes and underscores are dropped so the
    result is never a hidden file, never looks like a command-line option, and is
    never ``.`` or ``..``. A name that has nothing usable left becomes *default*.

    Args:
        name: any value (``None`` accepted) to be used as a file name part.
        default: returned when *name* holds no usable character.
        max_length: maximum length of the returned component.

    Returns:
        str: a non-empty string containing no path separator.

    Examples:
        >>> safe_filename_component("Leishmania_braziliensis_MHOM/BR/75/M2904")
        'Leishmania_braziliensis_MHOM_BR_75_M2904'
        >>> safe_filename_component("'Nostoc azollae' 0708")
        'Nostoc_azollae_0708'
        >>> safe_filename_component("   ")
        'unnamed'
    """
    if name is None:
        return default

    collapsed = _UNSAFE_RUN.sub("_", str(name).strip())
    trimmed = _TRAILING_JUNK.sub("", _LEADING_JUNK.sub("", collapsed))

    if len(trimmed) > max_length:
        trimmed = _TRAILING_JUNK.sub("", trimmed[:max_length])

    return trimmed or default


def is_safe_filename_component(name):
    """True when *name* can be used as a path component with no change."""
    return bool(name) and str(name) == safe_filename_component(name)

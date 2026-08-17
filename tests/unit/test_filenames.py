"""Tests for filesystem-safe names built from outside text.

Regression guard for the drawing crash on organisms whose strain designation
contains a slash::

    Leishmania braziliensis MHOM/BR/75/M2904
    Chlorella vulgaris CCAP 1055/1

The taxon name went into the karyotype file name unchanged, so the slash was read
as a path separator: the output "file name" became a directory hierarchy nobody
created and drawing died with FileNotFoundError — after the whole pipeline had
already written its data files.
"""

import os

import pytest

from satellome.core_functions.tools.filenames import (
    DEFAULT_NAME,
    MAX_COMPONENT_LENGTH,
    is_safe_filename_component,
    safe_filename_component,
)


# The five organism names that actually failed, as satellome passes them on
# (spaces already replaced with underscores in main.py).
FAILING_TAXA = [
    "Leishmania_braziliensis_MHOM/BR/75/M2904",
    "Chlorella_vulgaris_CCAP_1055/1",
    "Leishmania_panamensis_MHOM/GT/2001/U1103",
    "Fusarium_MF3/22",
    "Streptomyces_NIH/UT8656",
]


@pytest.mark.parametrize("taxon", FAILING_TAXA)
def test_slash_in_taxon_name_does_not_survive_into_a_file_name(taxon):
    safe = safe_filename_component(taxon)

    assert "/" not in safe
    assert os.sep not in safe
    # The name is a single path component, not a path.
    assert os.path.basename(safe) == safe
    assert os.path.dirname(safe) == ""


def test_the_real_name_is_kept_recognizable():
    assert (
        safe_filename_component("Leishmania_braziliensis_MHOM/BR/75/M2904")
        == "Leishmania_braziliensis_MHOM_BR_75_M2904"
    )
    assert safe_filename_component("Chlorella_vulgaris_CCAP_1055/1") == "Chlorella_vulgaris_CCAP_1055_1"


def test_names_that_are_already_safe_pass_through_untouched():
    for name in ["Homo_sapiens", "GCF_000005845.2_ASM584v2", "Escherichia-coli-K12", "Unknown"]:
        assert safe_filename_component(name) == name
        assert is_safe_filename_component(name)


def test_a_run_of_unsafe_characters_collapses_into_one_underscore():
    assert safe_filename_component("Genus sp. // strain ?? 1") == "Genus_sp._strain_1"


def test_quotes_and_shell_metacharacters_are_removed():
    # NCBI Candidatus names carry single quotes; they used to reach a shell
    # command line as -t '...' and break the quoting.
    assert safe_filename_component("'Nostoc azollae' 0708") == "Nostoc_azollae_0708"
    assert safe_filename_component("Genus sp. $(rm -rf /)") == "Genus_sp._rm_-rf"
    assert safe_filename_component("Genus\tsp.\nstrain") == "Genus_sp._strain"


def test_windows_separators_and_drive_letters_are_neutralized():
    safe = safe_filename_component(r"C:\taxa\Genus species")
    assert "\\" not in safe
    assert ":" not in safe
    assert safe == "C_taxa_Genus_species"


def test_result_is_never_hidden_never_an_option_and_never_a_directory_alias():
    assert safe_filename_component(".hidden") == "hidden"
    assert safe_filename_component("--force") == "force"
    assert safe_filename_component(".") == DEFAULT_NAME
    assert safe_filename_component("..") == DEFAULT_NAME
    assert safe_filename_component("/") == DEFAULT_NAME
    assert safe_filename_component("../../etc/passwd") == "etc_passwd"


def test_nothing_usable_falls_back_to_the_default():
    for name in [None, "", "   ", "???", "___"]:
        assert safe_filename_component(name) == DEFAULT_NAME
    assert safe_filename_component(None, default="Unknown") == "Unknown"


def test_non_string_values_are_accepted():
    assert safe_filename_component(9606) == "9606"


def test_long_names_are_capped_with_room_left_for_suffixes():
    long_name = "Candidatus_" + "Verylongstrainname" * 20
    safe = safe_filename_component(long_name)

    assert len(safe) == MAX_COMPONENT_LENGTH
    assert safe.startswith("Candidatus_")
    # A component plus the longest suffix the drawing step appends stays inside
    # the 255-byte limit every common filesystem enforces.
    assert len(safe + ".karyo.gaps.1000bp.enhanced.svg") < 255


def test_sanitizing_is_idempotent():
    for taxon in FAILING_TAXA + ["'Nostoc azollae' 0708", "...", "a" * 500]:
        once = safe_filename_component(taxon)
        assert safe_filename_component(once) == once
        assert is_safe_filename_component(once)


def test_is_safe_filename_component_rejects_the_failing_names():
    for taxon in FAILING_TAXA:
        assert not is_safe_filename_component(taxon)
    assert not is_safe_filename_component("")
    assert not is_safe_filename_component(None)

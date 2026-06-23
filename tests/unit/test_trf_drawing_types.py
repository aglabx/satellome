"""Tests for numeric-field typing in read_trf_file (karyotype-drawing crash).

Regression guard for the karyotype crash:

    TypeError: '>' not supported between instances of 'int' and 'str'

Records read from a .sat file come through csv.DictReader, so every field is a
*string*. The drawing code does ``max(length, enhance)`` and the clustering code
takes a *median* of periods; on strings the former crashes and the latter sorts
lexicographically ("100" < "20"). read_trf_file now coerces the coordinate /
length / period fields to int at the single point where records are built, so we
assert the concrete types and values, not merely "it ran".
"""

import logging

from satellome.core_functions.trf_drawing import read_trf_file, _to_int


FIELDS = [
    "project", "trf_id", "trf_head", "trf_l_ind", "trf_r_ind", "trf_period",
    "trf_n_copy", "trf_pmatch", "trf_pvar", "trf_entropy", "trf_consensus",
    "trf_array", "trf_array_gc", "trf_consensus_gc", "trf_array_length",
    "trf_joined", "trf_family", "trf_ref_annotation",
]
HEADER = "\t".join(FIELDS)


def _row(**overrides):
    base = {
        "project": "proj", "trf_id": "1", "trf_head": "chr1",
        "trf_l_ind": "1000", "trf_r_ind": "2000", "trf_period": "100",
        "trf_n_copy": "10", "trf_pmatch": "95", "trf_pvar": "5",
        "trf_entropy": "1.9", "trf_consensus": "AT", "trf_array": "ATATAT",
        "trf_array_gc": "40", "trf_consensus_gc": "50",
        "trf_array_length": "1000", "trf_joined": "0", "trf_family": "",
        "trf_ref_annotation": "",
    }
    base.update(overrides)
    return "\t".join(str(base[k]) for k in FIELDS)


def _write_sat(tmp_path, rows):
    path = tmp_path / "test.sat"
    path.write_text(HEADER + "\n" + "\n".join(rows) + "\n")
    return str(path)


def test_numeric_fields_are_ints_with_correct_values(tmp_path):
    sat = _write_sat(tmp_path, [_row(trf_l_ind="1000", trf_r_ind="2000",
                                     trf_period="100", trf_array_length="1500")])
    (record,) = read_trf_file(sat)

    assert record["start"] == 1000 and isinstance(record["start"], int)
    assert record["end"] == 2000 and isinstance(record["end"], int)
    assert record["period"] == 100 and isinstance(record["period"], int)
    assert record["length"] == 1500 and isinstance(record["length"], int)


def test_drawing_numeric_op_does_not_crash(tmp_path):
    """The exact crash path: max(record['length'], enhance) on a string raised."""
    sat = _write_sat(tmp_path, [_row(trf_array_length="1000")])
    (record,) = read_trf_file(sat)

    # Pre-fix this raised TypeError because length was the string "1000".
    enhanced = max(record["length"], 1_000_000)
    assert enhanced == 1_000_000


def test_period_sorts_numerically_not_lexicographically(tmp_path):
    """Median-period clustering must order periods numerically, not as strings."""
    sat = _write_sat(tmp_path, [
        _row(trf_id="1", trf_period="9"),
        _row(trf_id="2", trf_period="20"),
        _row(trf_id="3", trf_period="100"),
    ])
    periods = sorted(r["period"] for r in read_trf_file(sat))

    assert periods == [9, 20, 100]          # numeric order
    assert periods != ["100", "20", "9"]    # would be the lexicographic result


def test_missing_length_falls_back_to_zero(tmp_path):
    sat = _write_sat(tmp_path, [_row(trf_array_length="")])
    (record,) = read_trf_file(sat)
    assert record["length"] == 0


def test_to_int_coerces_and_surfaces_malformed(caplog):
    assert _to_int("123") == 123
    assert _to_int("12.0") == 12
    assert _to_int("") == 0
    assert _to_int(None) == 0
    assert _to_int(None, default=-1) == -1

    # A non-empty malformed value must be surfaced (visible warning), not silent.
    with caplog.at_level(logging.WARNING):
        assert _to_int("not-a-number") == 0
    assert any("not-a-number" in rec.message for rec in caplog.records)

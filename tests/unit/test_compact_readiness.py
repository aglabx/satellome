"""What counts as a directory it is safe to compact.

The obvious markers do not work, and the tests say why so that nobody
reintroduces them: 1,020 finished catalogues carry no ``.rc_complete``, and
1,015 carry a ``fasta/`` directory left over from an earlier campaign.
"""

import gzip
import os

from satellome.compact import readiness


def make_dir(tmp_path, prefix="ASM_v1", master=b"ok\n", gz=True, extras=True):
    run_dir = tmp_path / "tree" / "ACC.1"
    (run_dir / "fastan").mkdir(parents=True)
    if gz:
        with gzip.open(run_dir / f"{prefix}.sat.gz", "wb") as fh:
            fh.write(master)
    else:
        (run_dir / f"{prefix}.sat").write_bytes(master)
    if extras:
        (run_dir / "fastan" / f"{prefix}.bed").write_text("chr1\t0\t10\t5\t900\n")
        (run_dir / "results.yaml").write_text("pid: project\n")
    return str(run_dir)


def test_a_finished_directory_is_accepted(tmp_path):
    ready = readiness.check(make_dir(tmp_path))
    assert ready.ok, ready.reasons
    assert ready.prefix == "ASM_v1"


def test_the_master_is_found_by_suffix_not_by_the_directory_name(tmp_path):
    run_dir = make_dir(tmp_path, prefix="GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic")
    ready = readiness.check(run_dir)
    assert ready.ok, ready.reasons
    assert ready.prefix == "GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic"
    assert os.path.basename(run_dir) == "ACC.1"


def test_derived_views_are_not_mistaken_for_the_master(tmp_path):
    run_dir = make_dir(tmp_path)
    for suffix in ["1kb", "10kb", "micro", "pmicro", "tssr", "complex"]:
        with gzip.open(os.path.join(run_dir, f"ASM_v1.{suffix}.sat.gz"), "wb") as fh:
            fh.write(b"view\n")
    ready = readiness.check(run_dir)
    assert ready.ok, ready.reasons
    assert ready.prefix == "ASM_v1"


def test_a_staging_directory_is_refused_whatever_it_contains(tmp_path):
    staging = tmp_path / "_incoming"
    run_dir = make_dir(staging)
    ready = readiness.check(run_dir)
    assert not ready.ok
    assert "staging" in " ".join(ready.reasons)


def test_a_truncated_master_is_refused_with_the_reason(tmp_path):
    run_dir = make_dir(tmp_path)
    path = os.path.join(run_dir, "ASM_v1.sat.gz")
    data = open(path, "rb").read()
    open(path, "wb").write(data[:-6])
    ready = readiness.check(run_dir)
    assert not ready.ok
    assert any("does not decompress to the last byte" in r for r in ready.reasons)


def test_an_uncompressed_sat_beside_its_gz_means_a_run_is_still_writing(tmp_path):
    run_dir = make_dir(tmp_path)
    open(os.path.join(run_dir, "ASM_v1.sat"), "w").write("growing\n")
    ready = readiness.check(run_dir)
    assert not ready.ok
    assert any("sits uncompressed beside" in r for r in ready.reasons)


def test_a_leftover_fasta_directory_is_not_a_refusal(tmp_path):
    # 1,015 tree directories carry one from the earlier campaign.
    run_dir = make_dir(tmp_path)
    os.makedirs(os.path.join(run_dir, "fasta"))
    assert readiness.check(run_dir).ok


def test_a_missing_rc_complete_is_not_a_refusal(tmp_path):
    # 1,020 finished catalogues predate the marker.
    run_dir = make_dir(tmp_path)
    assert not os.path.exists(os.path.join(run_dir, ".rc_complete"))
    assert readiness.check(run_dir).ok


def test_a_leftover_partial_file_is_refused(tmp_path):
    run_dir = make_dir(tmp_path)
    open(os.path.join(run_dir, "ASM_v1.sat.gz.partial"), "w").write("x")
    ready = readiness.check(run_dir)
    assert not ready.ok
    assert any(".partial" in r for r in ready.reasons)


def test_a_missing_bed_is_refused_by_name(tmp_path):
    run_dir = make_dir(tmp_path, extras=False)
    ready = readiness.check(run_dir)
    assert not ready.ok
    assert any("fastan/*.bed" in r for r in ready.reasons)
    assert any("results.yaml" in r for r in ready.reasons)


def test_an_empty_master_is_refused(tmp_path):
    run_dir = make_dir(tmp_path, gz=False, master=b"")
    ready = readiness.check(run_dir)
    assert not ready.ok
    assert any("is empty" in r for r in ready.reasons)


def test_no_master_at_all_is_refused(tmp_path):
    run_dir = tmp_path / "tree" / "ACC.1"
    (run_dir / "fastan").mkdir(parents=True)
    ready = readiness.check(str(run_dir))
    assert not ready.ok
    assert any("no master" in r for r in ready.reasons)


def test_two_different_masters_are_refused_rather_than_guessed(tmp_path):
    run_dir = make_dir(tmp_path)
    with gzip.open(os.path.join(run_dir, "OTHER_v2.sat.gz"), "wb") as fh:
        fh.write(b"ok\n")
    ready = readiness.check(run_dir)
    assert not ready.ok
    assert any("cannot tell which .sat is the master" in r for r in ready.reasons)

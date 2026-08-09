"""Tests for the output-directory lock.

Two satellome runs pointed at the same -o truncate and interleave each other's
files: the second run's write reopens an artifact the first already finished. The
lock makes that an explicit refusal naming the holder, instead of a directory
whose contents belong to neither run.
"""

import json
import os

import pytest

from satellome.core_functions.tools.dir_lock import (
    LOCK_NAME,
    OutputDirLocked,
    acquire_output_lock,
    release_output_lock,
)


STARTED = "2026-08-09T12:00:00"


def test_lock_file_records_the_holder(tmp_path):
    path = acquire_output_lock(str(tmp_path), STARTED)

    assert path == str(tmp_path / LOCK_NAME)
    payload = json.loads((tmp_path / LOCK_NAME).read_text())
    assert payload["pid"] == os.getpid()
    assert payload["started"] == STARTED
    assert payload["host"]


def test_second_run_is_refused_and_told_who_holds_it(tmp_path):
    acquire_output_lock(str(tmp_path), STARTED)

    with pytest.raises(OutputDirLocked) as excinfo:
        acquire_output_lock(str(tmp_path), "2026-08-09T12:05:00")

    message = str(excinfo.value)
    assert str(os.getpid()) in message
    assert str(tmp_path) in message
    assert "--ignore-lock" in message


def test_release_allows_the_next_run(tmp_path):
    path = acquire_output_lock(str(tmp_path), STARTED)
    release_output_lock(path)

    assert not (tmp_path / LOCK_NAME).exists()
    assert acquire_output_lock(str(tmp_path), STARTED) == path


def test_release_is_idempotent(tmp_path):
    path = acquire_output_lock(str(tmp_path), STARTED)
    release_output_lock(path)
    release_output_lock(path)  # no raise
    release_output_lock(None)


def test_stale_lock_from_dead_pid_is_reclaimed(tmp_path, caplog):
    """A crashed run must not block the directory forever."""
    dead_pid = _find_dead_pid()
    (tmp_path / LOCK_NAME).write_text(json.dumps({
        "pid": dead_pid,
        "host": os.uname().nodename,
        "started": "2026-08-09T09:00:00",
        "argv": "-i genome.fasta -o out",
    }))

    with caplog.at_level("WARNING"):
        path = acquire_output_lock(str(tmp_path), STARTED)

    assert path == str(tmp_path / LOCK_NAME)
    assert json.loads((tmp_path / LOCK_NAME).read_text())["pid"] == os.getpid()
    assert any("stale output lock" in r.message for r in caplog.records)


def test_lock_from_another_host_is_not_reclaimed(tmp_path):
    """We cannot check a remote pid, so we refuse rather than guess it is dead."""
    (tmp_path / LOCK_NAME).write_text(json.dumps({
        "pid": 4242,
        "host": "some-other-cluster-node",
        "started": "2026-08-09T09:00:00",
        "argv": "-i genome.fasta -o out",
    }))

    with pytest.raises(OutputDirLocked, match="some-other-cluster-node"):
        acquire_output_lock(str(tmp_path), STARTED)


def test_unreadable_lock_is_refused_not_silently_overwritten(tmp_path):
    (tmp_path / LOCK_NAME).write_text("{ truncated json")

    with pytest.raises(OutputDirLocked, match="unreadable lock file"):
        acquire_output_lock(str(tmp_path), STARTED)

    # The foreign lock file is left exactly as it was.
    assert (tmp_path / LOCK_NAME).read_text() == "{ truncated json"


def test_ignore_lock_proceeds_and_warns(tmp_path, caplog):
    acquire_output_lock(str(tmp_path), STARTED)

    with caplog.at_level("WARNING"):
        path = acquire_output_lock(str(tmp_path), STARTED, ignore_lock=True)

    # Nothing for the ignoring run to release: the other run still owns the lock.
    assert path is None
    assert any("--ignore-lock" in r.message for r in caplog.records)
    assert json.loads((tmp_path / LOCK_NAME).read_text())["started"] == STARTED


def test_release_does_not_steal_a_reclaimed_lock(tmp_path, caplog):
    """If our lock was reclaimed as stale, releasing must not evict the new owner."""
    path = acquire_output_lock(str(tmp_path), STARTED)
    (tmp_path / LOCK_NAME).write_text(json.dumps({
        "pid": os.getpid() + 1,
        "host": os.uname().nodename,
        "started": "2026-08-09T13:00:00",
        "argv": "second run",
    }))

    with caplog.at_level("WARNING"):
        release_output_lock(path)

    assert (tmp_path / LOCK_NAME).exists()
    assert any("Not releasing" in r.message for r in caplog.records)


def _find_dead_pid():
    """A pid that is not currently running (searching downward from a high value)."""
    for pid in range(99999, 1, -1):
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return pid
        except PermissionError:
            continue
    pytest.skip("could not find an unused pid on this system")

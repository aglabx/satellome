""" ``.compact.json`` — the record without which compaction is deletion.

Every compacted directory carries one. It states, per file the compaction
touched: what class the policy gave it, what was done, the md5 and byte length
of the content **before**, the md5 and byte length **after**, and — for anything
dropped — the recipe that puts it back. ``expand`` reads nothing else.

Two properties are deliberate:

* The recorded digests are of the *content*, not of the bytes on disk. A file
  that was ``.tsv.gz`` before and is restored as ``.tsv.gz`` with a different
  gzip level still verifies, because gzip framing is not the thing we promised
  to preserve — the table is.
* A file whose recipe could not be proven at compaction time is never recorded
  as dropped, because it is never dropped. The engine regenerates and compares
  before it unlinks; the record therefore describes only reversible acts.
"""

import datetime
import json
import os
import platform
import socket

from satellome.compact.policy import POLICY_VERSION
from satellome.core_functions.tools.atomic_io import atomic_output

RECORD_NAME = ".compact.json"
SCHEMA_VERSION = 1


class RecordError(Exception):
    """The record is missing, unreadable, or not of a shape we can act on."""


def record_path(run_dir):
    return os.path.join(run_dir, RECORD_NAME)


def is_compacted(run_dir):
    return os.path.exists(record_path(run_dir))


def new_record(run_dir, config, satellome_version, prefix=None):
    """Start an empty record for *run_dir*."""
    return {
        "schema_version": SCHEMA_VERSION,
        "policy_version": POLICY_VERSION,
        "satellome_version": satellome_version,
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "run_dir": os.path.abspath(run_dir),
        "prefix": prefix,
        "config": {
            "min_array_length": config.min_array_length,
            "level": config.level,
            "codec": "zstd",
            "verify_drops": config.verify_drops,
        },
        "files": [],
        "kept_unclassified": [],
        "refusals": [],
        "totals": {"before": 0, "after": 0, "freed": 0},
    }


def add_entry(record, entry):
    record["files"].append(entry)
    record["totals"]["before"] += entry.get("before", {}).get("stored_bytes", 0)
    record["totals"]["after"] += entry.get("after", {}).get("stored_bytes", 0)
    record["totals"]["freed"] = record["totals"]["before"] - record["totals"]["after"]
    return entry


def write_record(run_dir, record):
    """Write ``.compact.json`` atomically."""
    path = record_path(run_dir)
    with atomic_output(path) as fh:
        json.dump(record, fh, indent=2, sort_keys=True)
        fh.write("\n")
    return path


def load_record(run_dir):
    """Load and shape-check ``.compact.json``.

    Raises:
        RecordError: absent, unreadable, invalid JSON, wrong shape, or written
            by a schema this build does not understand. A record we cannot read
            is never treated as "not compacted" — that would let ``expand``
            report success on a directory it did not restore.
    """
    path = record_path(run_dir)
    if not os.path.exists(path):
        raise RecordError(
            f"No {RECORD_NAME} in {run_dir}: this directory was not compacted "
            f"by satellome (or the record was removed)"
        )
    try:
        with open(path, "r") as fh:
            record = json.load(fh)
    except json.JSONDecodeError as e:
        raise RecordError(f"Corrupt {RECORD_NAME} in {run_dir}: {e}") from e
    except OSError as e:
        raise RecordError(f"Unreadable {RECORD_NAME} in {run_dir}: {e}") from e

    if not isinstance(record, dict) or not isinstance(record.get("files"), list):
        raise RecordError(
            f"Corrupt {RECORD_NAME} in {run_dir}: expected an object with a "
            f"'files' list, got {type(record).__name__}"
        )
    schema = record.get("schema_version")
    if schema != SCHEMA_VERSION:
        raise RecordError(
            f"{RECORD_NAME} in {run_dir} has schema_version {schema}, this "
            f"satellome understands {SCHEMA_VERSION}"
        )
    return record


def entries_by_action(record, action):
    return [e for e in record["files"] if e.get("action") == action]

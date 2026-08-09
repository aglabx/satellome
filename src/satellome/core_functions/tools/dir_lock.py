"""Exclusive lock on an output directory.

Why this exists
---------------
Two satellome runs pointed at the same ``-o`` directory (a pool worker and a
retry of the same genome, say) overwrite each other's files mid-write. The
second run's ``open(path, 'w')`` truncates an artifact the first run already
finished, and both processes then interleave writes into the same paths. On
disk the result is a file that is neither run's output, with nothing to indicate
it. Atomic publication fixes the *within-run* window; only a lock can stop two
runs from fighting over the same names.

The lock is a file created with ``O_EXCL``, holding the owner's pid, host and
command line so a collision reports *who* holds it rather than just "busy".

A lock left by a process that no longer exists (crash, ``kill -9``) is stale.
We only ever reclaim such a lock when we can prove it: same host, and the
recorded pid is gone. A lock from another host cannot be checked from here, so
it is reported and the run stops — the escape hatch is explicit
(``--ignore-lock``), never implicit.
"""

import json
import logging
import os
import socket
import sys

logger = logging.getLogger("satellome")

LOCK_NAME = ".satellome.lock"


class OutputDirLocked(Exception):
    """Another satellome run holds the lock on this output directory."""


def _lock_path(output_dir):
    return os.path.join(output_dir, LOCK_NAME)


def _read_lock(path):
    """Return the lock's payload dict, or None if it is unreadable/garbage."""
    try:
        with open(path, "r") as fh:
            payload = json.load(fh)
    except (json.JSONDecodeError, OSError) as e:
        logger.warning(f"Lock file {path} is unreadable or malformed: {e}")
        return None
    return payload if isinstance(payload, dict) else None


def _pid_alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        # Exists but owned by another user.
        return True
    except (OverflowError, TypeError, ValueError):
        return False
    return True


def _describe(payload):
    if not payload:
        return "an unreadable lock file"
    return (
        f"pid {payload.get('pid', '?')} on {payload.get('host', '?')} "
        f"started {payload.get('started', 'at an unknown time')} "
        f"({payload.get('argv', '?')})"
    )


def acquire_output_lock(output_dir, started, ignore_lock=False):
    """Take the lock on ``output_dir``.

    Args:
        output_dir: directory to lock (must exist).
        started: ISO timestamp string recorded in the lock, supplied by the
            caller so this module stays free of clock access.
        ignore_lock: proceed even if another live run holds the lock. Logged as
            a warning, because concurrent runs into one directory can still
            corrupt each other's outputs.

    Returns:
        Path of the lock file, or None when an existing lock was ignored (there
        is then nothing for us to release — the other run owns it).

    Raises:
        OutputDirLocked: a live (or unverifiable) run holds the lock.
    """
    path = _lock_path(output_dir)
    payload = {
        "pid": os.getpid(),
        "host": socket.gethostname(),
        "started": started,
        "argv": " ".join(sys.argv[1:]),
    }

    while True:
        try:
            fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            pass
        else:
            with os.fdopen(fd, "w") as fh:
                json.dump(payload, fh, indent=2, sort_keys=True)
                fh.write("\n")
            logger.debug(f"Acquired output lock: {path}")
            return path

        existing = _read_lock(path)
        same_host = bool(existing) and existing.get("host") == payload["host"]
        pid = existing.get("pid") if existing else None
        stale = same_host and isinstance(pid, int) and not _pid_alive(pid)

        if stale:
            logger.warning(
                f"Reclaiming stale output lock in {output_dir}: "
                f"{_describe(existing)} is no longer running"
            )
            try:
                os.unlink(path)
            except FileNotFoundError:
                pass
            continue  # retry the O_EXCL create

        message = (
            f"Output directory {output_dir} is locked by {_describe(existing)}. "
            f"Two runs writing the same output directory overwrite each other's "
            f"files. Use a different -o, wait for that run, or pass --ignore-lock "
            f"if you are certain it is safe."
        )
        if ignore_lock:
            logger.warning(f"--ignore-lock: proceeding anyway. {message}")
            return None
        raise OutputDirLocked(message)


def release_output_lock(lock_path):
    """Release a lock previously taken by this process. Best effort, idempotent."""
    if not lock_path:
        return
    payload = _read_lock(lock_path)
    if payload and payload.get("pid") != os.getpid():
        # Someone reclaimed it as stale while we were alive: deleting it would
        # yank the lock out from under its new owner.
        logger.warning(
            f"Not releasing {lock_path}: it is now held by {_describe(payload)}"
        )
        return
    try:
        os.unlink(lock_path)
    except FileNotFoundError:
        pass
    except OSError as e:
        logger.warning(f"Could not remove lock file {lock_path}: {e}")

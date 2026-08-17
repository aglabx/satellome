"""Atomic publication of pipeline output files.

Why this exists
---------------
Every satellome output used to be written *in place*: ``open(path, 'w')`` on the
final path, then streamed for minutes (a ``.sat`` for a T2T genome is
gigabytes). That has two consequences, both observed in production:

1. A concurrent reader — a driver script that gzips the output directory as soon
   as satellome's *previous* stage looks done — sees a file that is still
   growing. ``gzip`` warns ``file size changed while zipping`` and produces a
   *structurally valid* ``.gz`` of a *truncated* file: ``gzip -t`` passes, the
   file exists, and the truncation is invisible to any existence check.
2. A rerun into the same output directory truncates an already-valid file at
   ``open(..., 'w')`` and rebuilds it. Between those two moments the
   authoritative artifact on disk is partial.

Writing to ``<path>.partial`` and ``os.replace``-ing it into place at the end
makes the final name appear only when the file is complete: ``os.replace`` is
atomic within a filesystem, so a reader either sees the previous complete
version or the new complete version, never a half-written one. Existence of the
final path becomes a meaningful signal again.

A failed write leaves no final file at all — the ``.partial`` is removed and the
exception propagates. We never publish what we could not finish writing.
"""

import contextlib
import logging
import os

logger = logging.getLogger("satellome")

PARTIAL_SUFFIX = ".partial"


def partial_path(path):
    """Return the in-progress path used while writing ``path``."""
    return f"{path}{PARTIAL_SUFFIX}"


def publish(path):
    """Move ``<path>.partial`` onto ``path`` atomically.

    Used by the paths where an external binary (e.g. the Rust ``bed-extract``)
    writes the file for us: we hand it the ``.partial`` name and publish after
    it exits successfully.
    """
    tmp = partial_path(path)
    if not os.path.exists(tmp):
        raise FileNotFoundError(
            f"Cannot publish {path}: expected in-progress file {tmp} is missing"
        )
    _fsync_file(tmp)
    os.replace(tmp, path)
    return path


def discard(path):
    """Remove a leftover ``<path>.partial``, if any. Best effort."""
    tmp = partial_path(path)
    try:
        os.unlink(tmp)
    except FileNotFoundError:
        pass
    except OSError as e:
        # Not fatal for the caller, but must not be invisible: a stale
        # .partial is a real thing someone will find on disk later.
        logger.warning(f"Could not remove in-progress file {tmp}: {e}")


@contextlib.contextmanager
def atomic_output(path, mode="w", **open_kwargs):
    """Context manager yielding a handle to ``<path>.partial``.

    On clean exit the data is flushed, fsynced and ``os.replace``-d onto
    ``path``. On exception the partial file is removed and ``path`` is left
    untouched (a previous complete version survives), then the exception
    propagates — a write failure is never swallowed.

    ``mode`` must be a writing mode; anything else is a programming error.
    """
    if "r" in mode and "+" not in mode:
        raise ValueError(f"atomic_output is for writing, got mode={mode!r}")

    tmp = partial_path(path)
    _ensure_parent(tmp)

    fh = open(tmp, mode, **open_kwargs)
    try:
        yield fh
    except BaseException:
        fh.close()
        discard(path)
        raise
    else:
        fh.flush()
        try:
            os.fsync(fh.fileno())
        except OSError as e:
            # fsync can legitimately fail on some filesystems (e.g. certain
            # network mounts). Durability is degraded, atomicity of the rename
            # is not — say so instead of pretending it happened.
            logger.warning(f"fsync failed for {tmp}: {e} (data may not be durable)")
        fh.close()
        os.replace(tmp, path)


@contextlib.contextmanager
def atomic_outputs(*paths, mode="w", **open_kwargs):
    """Like :func:`atomic_output` but for several files published together.

    Yields a tuple of handles, one per path (``None`` for a ``None`` path, so
    callers with an optional output do not need a special case). All files are
    published only after the body completes; if the body raises, none of them
    are, and every partial is removed. This keeps sibling outputs — a ``.sat``
    and the ``.fasta`` extracted in the same pass — from disagreeing with each
    other on disk.
    """
    handles = []
    opened = []
    try:
        for path in paths:
            if path is None:
                handles.append(None)
                continue
            tmp = partial_path(path)
            _ensure_parent(tmp)
            fh = open(tmp, mode, **open_kwargs)
            handles.append(fh)
            opened.append((path, fh))
        yield tuple(handles)
    except BaseException:
        for path, fh in opened:
            fh.close()
            discard(path)
        raise
    else:
        for path, fh in opened:
            fh.flush()
            try:
                os.fsync(fh.fileno())
            except OSError as e:
                logger.warning(f"fsync failed for {partial_path(path)}: {e}")
            fh.close()
        for path, _ in opened:
            os.replace(partial_path(path), path)


def _ensure_parent(tmp):
    """Make sure the directory of ``tmp`` exists, reporting if we had to build it.

    The pipeline creates all of its output subdirectories up front, so a missing
    parent here means the *file name* contained a path separator — an organism
    name like ``MHOM/BR/75/M2904`` interpolated into a path, for instance. Making
    the directory silently would scatter outputs into a hierarchy nobody looks
    in, which is worse than the crash it replaces: say so.
    """
    parent = os.path.dirname(os.path.abspath(tmp))
    if parent and not os.path.isdir(parent):
        logger.warning(
            f"Creating missing output directory {parent} for {os.path.basename(tmp)} "
            f"— check the file name for stray path separators"
        )
        os.makedirs(parent, exist_ok=True)


def _fsync_file(path):
    """fsync a closed file by reopening it; warns rather than raising."""
    try:
        fd = os.open(path, os.O_RDONLY)
    except OSError as e:
        logger.warning(f"Could not open {path} to fsync: {e}")
        return
    try:
        os.fsync(fd)
    except OSError as e:
        logger.warning(f"fsync failed for {path}: {e} (data may not be durable)")
    finally:
        os.close(fd)

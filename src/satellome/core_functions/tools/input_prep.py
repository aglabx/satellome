"""Accept compressed and .2bit genomes, and hand downstream tools plain FASTA.

Satellome's own readers cope with a few compressed forms, but the pipeline
shells out to FasTAN, TRF and the Rust helpers with a *path*, and those read
plain FASTA. So the rule is: detect what the user gave us, and if it is not
something every downstream tool can open, materialise a plain FASTA once and
use that everywhere.

Deliberately not silent: the conversion is announced, its result is cached
between runs, and a corrupt or unreadable archive is an error rather than a
fallback to "treat the bytes as FASTA" — which would otherwise surface much
later as a genome with zero tandem repeats.
"""

import bz2
import gzip
import logging
import lzma
import os
import shutil
import zipfile
from pathlib import Path
from typing import Optional, Tuple

from satellome.core_functions.io.twobit_file import TwoBitError, is_twobit, twobit_to_fasta
from satellome.core_functions.tools.atomic_io import atomic_output

logger = logging.getLogger(__name__)

# Magic numbers, because an extension is a hint and the bytes are the truth:
# a .fa that is really gzip, or a .2bit named .fasta, both happen in practice.
MAGIC = {
    b"\x1f\x8b": "gzip",
    b"BZh": "bzip2",
    b"\xfd7zXZ": "xz",
    b"\x28\xb5\x2f\xfd": "zstd",
    b"PK\x03\x04": "zip",
}

PLAIN = "fasta"
CONVERTED_SUFFIX = ".satellome.fasta"


class InputFormatError(Exception):
    """The input is not a genome we can read."""


def sniff_format(path) -> str:
    """Identify the container by content: gzip/bzip2/xz/zstd/zip/2bit/fasta."""
    path = str(path)
    try:
        with open(path, "rb") as handle:
            head = handle.read(8)
    except OSError as error:
        raise InputFormatError(f"cannot read {path}: {error}") from error

    if not head:
        raise InputFormatError(f"{path} is empty")

    if is_twobit(path):
        return "2bit"
    for magic, name in MAGIC.items():
        if head.startswith(magic):
            return name
    return PLAIN


def needs_conversion(path) -> bool:
    """True when downstream tools cannot be handed this path directly."""
    return sniff_format(path) != PLAIN


def _open_compressed(path: str, kind: str):
    if kind == "gzip":
        return gzip.open(path, "rb")
    if kind == "bzip2":
        return bz2.open(path, "rb")
    if kind == "xz":
        return lzma.open(path, "rb")
    if kind == "zstd":
        try:
            import zstandard
        except ImportError as error:
            raise InputFormatError(
                f"{path} is zstd-compressed but the 'zstandard' package is not "
                "installed. Install it (pip install zstandard) or decompress "
                "the file first."
            ) from error
        return zstandard.open(path, "rb")
    raise InputFormatError(f"{path}: unsupported container '{kind}'")


def _single_fasta_member(archive: zipfile.ZipFile, path: str) -> str:
    """The one FASTA inside a zip, or a clear error naming what was found."""
    members = [n for n in archive.namelist() if not n.endswith("/")]
    candidates = [
        n for n in members
        if n.lower().endswith((".fa", ".fasta", ".fna", ".fas", ".seq"))
    ]
    if len(candidates) == 1:
        return candidates[0]
    if not candidates and len(members) == 1:
        return members[0]
    if not candidates:
        raise InputFormatError(
            f"{path}: no FASTA inside the archive (members: {', '.join(members[:5])})"
        )
    raise InputFormatError(
        f"{path}: the archive holds {len(candidates)} FASTA files "
        f"({', '.join(candidates[:5])}); unpack the one you want and pass it directly"
    )


def converted_path(path, work_dir) -> Path:
    """Where the plain-FASTA form of ``path`` is cached."""
    stem = Path(str(path)).name
    for suffix in (".gz", ".bz2", ".xz", ".zst", ".zip", ".2bit"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return Path(work_dir) / f"{stem}{CONVERTED_SUFFIX}"


def ensure_plain_fasta(path, work_dir, force: bool = False) -> Tuple[str, Optional[str]]:
    """Return a path every downstream tool can read, converting if needed.

    Returns ``(usable_path, source_format)`` where ``source_format`` is None
    when the input was already plain FASTA and nothing was written.

    The conversion is cached in ``work_dir`` and reused on the next run unless
    the source is newer (or ``force``), so re-running a pipeline on a
    compressed genome does not pay for decompression twice.
    """
    path = str(path)
    kind = sniff_format(path)
    if kind == PLAIN:
        return path, None

    target = converted_path(path, work_dir)
    if target.exists() and not force:
        try:
            fresh = target.stat().st_mtime >= os.path.getmtime(path)
        except OSError:
            fresh = False
        if fresh and target.stat().st_size > 0:
            logger.info(f"Using cached decompressed input: {target}")
            return str(target), kind
        logger.info(f"Cached input is stale, re-converting: {target}")

    Path(work_dir).mkdir(parents=True, exist_ok=True)
    logger.info(f"Input is {kind}; converting to plain FASTA at {target}")

    if kind == "2bit":
        try:
            written = twobit_to_fasta(path, target)
        except TwoBitError as error:
            raise InputFormatError(str(error)) from error
        logger.info(f"Converted {written} sequences from .2bit")
        return str(target), kind

    try:
        if kind == "zip":
            with zipfile.ZipFile(path) as archive:
                member = _single_fasta_member(archive, path)
                logger.info(f"Extracting '{member}' from the archive")
                with archive.open(member) as source, atomic_output(str(target), "wb") as out:
                    shutil.copyfileobj(source, out, length=1 << 20)
        else:
            with _open_compressed(path, kind) as source, atomic_output(str(target), "wb") as out:
                shutil.copyfileobj(source, out, length=1 << 20)
    except InputFormatError:
        raise
    except (OSError, EOFError, zipfile.BadZipFile, gzip.BadGzipFile, lzma.LZMAError) as error:
        # A truncated or corrupt archive must stop the run here. Falling back to
        # reading the raw bytes would produce a genome with no tandem repeats
        # and no visible reason why.
        raise InputFormatError(f"{path}: cannot decompress ({kind}): {error}") from error

    if target.stat().st_size == 0:
        raise InputFormatError(f"{path}: decompressed to an empty file")

    with open(target, "rb") as handle:
        if not handle.read(1).startswith(b">"):
            raise InputFormatError(
                f"{path}: decompressed content does not start with '>' - "
                "it does not look like FASTA"
            )

    logger.info(f"Decompressed to {target} ({target.stat().st_size / 1e9:.2f} GB)")
    return str(target), kind

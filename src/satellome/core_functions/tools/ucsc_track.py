"""Emit a UCSC Genome Browser track from a finished .sat file.

Produced on request (``--ucsc-track``), because most runs never need it and it
costs a full pass over the repeats.

Output is BED9 with per-item colour, which is what the browser needs to show
tandem repeats as something other than an undifferentiated black carpet: the
colour carries the period class, so microsatellites and satellite-sized arrays
are distinguishable at a glance at chromosome zoom.

Two files are written:

* ``<project>.ucsc.bed`` - a custom track, ready to paste into "add custom
  track" or to serve over HTTP.
* ``<project>.chrom.sizes`` - alongside, so the BED can be turned into a
  bigBed. Whole-genome BED custom tracks hit the browser's size limit; bigBed
  does not, and it is what a track hub needs. If ``bedToBigBed`` is on PATH we
  run it and write ``<project>.ucsc.bb`` as well.

Coordinates: the .sat file stores 0-based half-open [l_ind, r_ind), which is
already BED's convention, so positions pass through unchanged. (The GFF3
writer is the one that shifts, because GFF is 1-based inclusive.)
"""

import logging
import os
import shutil
import subprocess
from typing import Dict, Iterable, List, Optional, Tuple

from satellome.core_functions.io.tab_file import sc_iter_tab_file
from satellome.core_functions.models.trf_model import TRModel
from satellome.core_functions.tools.atomic_io import atomic_output

logger = logging.getLogger(__name__)

# Period classes and their colours. Chosen to stay distinguishable both from
# each other and against the browser's own black features, and to survive the
# grey-blue of a dense view.
PERIOD_CLASSES: List[Tuple[int, int, str, str]] = [
    (1, 6, "230,25,75", "microsatellite (1-6 bp)"),
    (7, 9, "245,130,48", "short SSR (7-9 bp)"),
    (10, 100, "60,180,75", "minisatellite (10-100 bp)"),
    (101, 1000, "0,130,200", "satellite (101-1000 bp)"),
    (1001, 10 ** 9, "145,30,180", "macrosatellite (>1000 bp)"),
]

DEFAULT_TRACK_NAME = "satellome_TRs"


def colour_for_period(period: int) -> str:
    for low, high, colour, _ in PERIOD_CLASSES:
        if low <= period <= high:
            return colour
    return "128,128,128"


def class_label(period: int) -> str:
    for low, high, _, label in PERIOD_CLASSES:
        if low <= period <= high:
            return label
    return "unclassified"


def _score(trf) -> int:
    """BED score 0-1000. Percent identity is the honest signal we have."""
    try:
        pmatch = float(getattr(trf, "trf_pmatch", 0) or 0)
    except (TypeError, ValueError):
        return 0
    return max(0, min(1000, int(round(pmatch * 10))))


def _chrom(trf) -> Optional[str]:
    """The sequence name exactly as the assembly spells it.

    This must be trf_head, not the trf_chr property: trf_chr normalises names
    for analysis ("chr1" -> "1"), which is right for grouping and wrong here.
    A track whose names do not match chrom.sizes shows nothing in the browser,
    and bedToBigBed fails outright with "chrom not found" - so the same string
    that goes into chrom.sizes (the first token of the FASTA header) is the
    only correct choice.
    """
    value = getattr(trf, "trf_head", None)
    if value and str(value) not in ("?", "Unknown"):
        # a FASTA header may carry a description after the first space
        head = str(value).lstrip(">").split()
        if head:
            return head[0]
    return None


def write_chrom_sizes(fasta_file: str, output_file: str) -> Dict[str, int]:
    """Write a UCSC chrom.sizes for the genome. Returns name -> length."""
    from satellome.core_functions.io.fasta_file import sc_iter_fasta_brute

    sizes: Dict[str, int] = {}
    for header, sequence in sc_iter_fasta_brute(fasta_file):
        name = header.lstrip(">").split()[0]
        sizes[name] = len(sequence)

    with atomic_output(output_file, "w") as handle:
        for name, length in sizes.items():
            handle.write(f"{name}\t{length}\n")
    return sizes


def iter_track_rows(trf_file: str, min_length: int = 0) -> Iterable[Tuple]:
    """Yield BED9 tuples for every repeat in the .sat file."""
    for trf in sc_iter_tab_file(trf_file, TRModel):
        chrom = _chrom(trf)
        if chrom is None:
            continue

        start, end = int(trf.trf_l_ind), int(trf.trf_r_ind)
        if start > end:
            start, end = end, start
        if end - start < min_length:
            continue

        period = int(getattr(trf, "trf_period", 0) or 0)
        family = getattr(trf, "trf_family", None) or f"period_{period}"
        name = f"{family}_{period}bp".replace(" ", "_")

        yield (
            chrom, start, end, name, _score(trf), "+",
            start, end, colour_for_period(period),
        )


def write_ucsc_bed(
    trf_file: str,
    output_file: str,
    track_name: str = DEFAULT_TRACK_NAME,
    description: Optional[str] = None,
    min_length: int = 0,
    with_track_line: bool = True,
) -> int:
    """Write a sorted, coloured BED9 custom track. Returns the row count.

    Rows are sorted by (chrom, start): the browser tolerates unsorted input but
    ``bedToBigBed`` refuses it, and writing it sorted here means the same file
    serves both uses.
    """
    rows = sorted(iter_track_rows(trf_file, min_length=min_length),
                  key=lambda r: (r[0], r[1], r[2]))

    with atomic_output(output_file, "w") as handle:
        if with_track_line:
            text = description or "Tandem repeats found by Satellome, coloured by period class"
            handle.write(
                f'track name="{track_name}" description="{text}" '
                'itemRgb="On" visibility=dense useScore=0\n'
            )
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")

    logger.info(f"UCSC track: {len(rows)} features -> {output_file}")
    return len(rows)


def make_bigbed(bed_file: str, chrom_sizes_file: str, output_file: str) -> bool:
    """Convert to bigBed with UCSC's bedToBigBed, if it is available.

    Returns False when the tool is absent (an expected, optional dependency) or
    when it fails - in which case its own stderr is reported, because "bigBed
    generation failed" alone tells the user nothing they can act on.
    """
    tool = shutil.which("bedToBigBed")
    if not tool:
        logger.info(
            "bedToBigBed not found; wrote the BED track only. Install it from "
            "http://hgdownload.soe.ucsc.edu/admin/exe/ to also get a bigBed "
            "(needed for whole-genome tracks and track hubs)."
        )
        return False

    # bedToBigBed rejects a track line, so feed it the rows only.
    stripped = bed_file + ".notrack"
    try:
        with open(bed_file) as source, open(stripped, "w") as target:
            for line in source:
                if not line.startswith(("track", "browser", "#")):
                    target.write(line)

        result = subprocess.run(
            [tool, "-type=bed9", stripped, chrom_sizes_file, output_file],
            capture_output=True, text=True, timeout=3600,
        )
        if result.returncode != 0:
            logger.error(f"bedToBigBed failed: {result.stderr.strip()[-1000:]}")
            return False
    except (OSError, subprocess.SubprocessError) as error:
        logger.error(f"bedToBigBed could not be run: {error}")
        return False
    finally:
        if os.path.exists(stripped):
            os.remove(stripped)

    logger.info(f"UCSC bigBed: {output_file}")
    return True


def build_ucsc_track(
    trf_file: str,
    fasta_file: str,
    output_dir: str,
    project: str,
    min_length: int = 0,
) -> Dict[str, str]:
    """Write the track files. Returns {kind: path} for what was produced."""
    produced: Dict[str, str] = {}

    bed_path = os.path.join(output_dir, f"{project}.ucsc.bed")
    count = write_ucsc_bed(
        trf_file, bed_path,
        track_name=f"{project}_TRs",
        description=f"Satellome tandem repeats for {project}, coloured by period class",
        min_length=min_length,
    )
    produced["bed"] = bed_path
    if count == 0:
        # An empty track is a real answer for a genome with no repeats above the
        # cutoff, but it must not look like a successful full track.
        logger.warning(
            f"UCSC track is empty: no tandem repeats in {trf_file} "
            f"passed min_length={min_length}"
        )

    sizes_path = os.path.join(output_dir, f"{project}.chrom.sizes")
    try:
        sizes = write_chrom_sizes(fasta_file, sizes_path)
        produced["chrom_sizes"] = sizes_path
        logger.info(f"chrom.sizes: {len(sizes)} sequences -> {sizes_path}")
    except OSError as error:
        logger.error(f"Could not write chrom.sizes: {error}")
        return produced

    bb_path = os.path.join(output_dir, f"{project}.ucsc.bb")
    if make_bigbed(bed_path, sizes_path, bb_path):
        produced["bigbed"] = bb_path

    return produced

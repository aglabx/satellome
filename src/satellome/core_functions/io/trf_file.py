#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 05.06.2011
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

"""
TRF (Tandem Repeat Finder) file I/O and processing module.

Provides parsers, filters, and converters for TRF output files. Handles both
raw TRF .dat format and tab-delimited parsed format. Includes consensus
canonicalization, overlap resolution, and format conversion utilities.

Classes:
    TRFFileIO: Main parser for raw TRF .dat files with filtering capabilities

Functions:
    join_overlapped: Merge two overlapping tandem repeats if similar
    get_int_gc: Calculate GC content as integer percentage
    get_shifts_variants: Generate all circular rotations of a sequence
    sort_dictionary_by_value: Sort dictionary by values
    remove_consensus_redundancy: Canonicalize consensus sequences to minimal form
    sc_parse_raw_trf_folder: Batch parse all TRF files in a folder
    sc_trf_to_fasta: Convert parsed TRF file to FASTA format

Key Features:
    - Memory-efficient streaming parser for large genomes
    - Automatic overlap detection and merging of similar repeats
    - Consensus sequence canonicalization (minimal lexicographic form)
    - Filtering of nested, duplicate, and dissimilar overlaps
    - Multimer detection (e.g., ACTACTACT → ACT)
    - Tab-delimited and FASTA output formats

Example:
    >>> # Parse TRF output and write to tab-delimited file
    >>> reader = TRFFileIO()
    >>> reader.parse_to_file("genome.dat", "genome_repeats.tab", project="hg38")
    >>>
    >>> # Convert to FASTA
    >>> sc_trf_to_fasta("genome_repeats.tab", "repeats.fasta")
    >>>
    >>> # Batch process folder
    >>> sc_parse_raw_trf_folder("trf_output/", "all_repeats.tab", project="hg38")

See Also:
    satellome.core_functions.models.trf_model: TRModel data model
    satellome.core_functions.io.tab_file: Tab-delimited file I/O
"""
import logging
import os
from collections import defaultdict

from satellome.core_functions.exceptions import SequenceError
from satellome.core_functions.io.abstract_reader import WiseOpener

logger = logging.getLogger(__name__)

from satellome.core_functions.io.block_file import AbstractBlockFileIO
from satellome.core_functions.io.file_system import iter_filepath_folder
from satellome.core_functions.io.tab_file import sc_iter_tab_file
from satellome.core_functions.models.trf_model import TRModel
from satellome.core_functions.tools.parsers import refine_name
from satellome.core_functions.tools.processing import get_revcomp, get_gc_content
from satellome.core_functions.trf_embedings import get_cosine_distance


def join_overlapped(obj1, obj2, cutoff_distance=0.1):
    """
    Join two overlapping tandem repeat objects if they are similar.

    Determines if two TR objects should be merged based on:
    1. Genomic overlap (obj1.trf_r_ind > obj2.trf_l_ind)
    2. Sequence similarity via cosine distance of 5-mer vectors
    3. Intersection fraction threshold

    Objects are joined if either:
    - Cosine distance < cutoff_distance (default 0.1), OR
    - Intersection fraction > 0.2

    Args:
        obj1 (TRModel): First tandem repeat object (left position)
        obj2 (TRModel): Second tandem repeat object (right position)
        cutoff_distance (float, optional): Maximum cosine distance for merging.
                                          Defaults to 0.1.

    Returns:
        bool: True if objects were joined (obj1 extended with obj2 data),
              False if not overlapping or too dissimilar

    Note:
        - Modifies obj1 in-place via set_form_overlap(obj2)
        - Uses 5-mer embedding vectors for similarity comparison
        - Requires obj1.trf_r_ind > obj2.trf_l_ind (overlap condition)
    """
    # a ------
    # b    -----
    if obj1.trf_r_ind > obj2.trf_l_ind and obj1.trf_r_ind < obj2.trf_r_ind:
        vector1 = obj1.get_vector()
        vector2 = obj2.get_vector()
        dist = get_cosine_distance(vector1, vector2)
        left_part = obj2.trf_l_ind - obj1.trf_l_ind
        right_part = obj2.trf_r_ind - obj1.trf_r_ind
        middle_part = obj1.trf_r_ind - obj2.trf_l_ind

        intersect_fraction = middle_part / (left_part + right_part + middle_part)

        if dist < cutoff_distance or intersect_fraction > 0.2:
            obj1.set_form_overlap(obj2)
            return True
    return False


def get_int_gc(sequence):
    """
    Calculate GC content as integer percentage (0-100).

    Args:
        sequence (str): DNA sequence

    Returns:
        int: GC content as integer from 0 to 100

    Example:
        >>> get_int_gc("ATGC")
        50
        >>> get_int_gc("AAAA")
        0
    """
    gc = get_gc_content(sequence)
    return int(100 * round(gc, 2))


def get_shifts_variants(sequence):
    """
    Generate all circular rotations of a sequence.

    Args:
        sequence (str): Input sequence (e.g., DNA monomer)

    Returns:
        list: All unique circular rotations of the sequence

    Raises:
        SequenceError: If sequence cannot generate shift variants

    Example:
        >>> get_shifts_variants("ACGT")
        ['ACGT', 'CGTA', 'GTAC', 'TACG']
    """
    shifts = set()
    for i in range(len(sequence)):
        shifts.add(sequence[i:] + sequence[:i])
    return list(shifts)


def sort_dictionary_by_value(d, reverse=False):
    """
    Sort dictionary by values in ascending or descending order.

    Args:
        d (dict): Dictionary to sort
        reverse (bool, optional): If True, sort descending. Defaults to False.

    Returns:
        list: List of (value, key) tuples sorted by value

    Example:
        >>> sort_dictionary_by_value({"a": 3, "b": 1, "c": 2})
        [(1, 'b'), (2, 'c'), (3, 'a')]
    """
    result = [(v, k) for k, v in list(d.items())]
    result.sort(reverse=reverse, key=lambda x: x[0])
    return result


def remove_consensus_redundancy(trf_objs):
    """
    Canonicalize tandem repeat consensus sequences to minimal lexicographic form.

    For each consensus sequence, generates all circular rotations and reverse
    complements, then selects the lexicographically minimal variant as the
    canonical form. Also detects and collapses multimers (e.g., "GTAGTAGTA"
    becomes "ACT" if ACT is the minimal form).

    Canonicalization ensures that equivalent repeat units (rotations and reverse
    complements) share the same consensus representation, enabling proper
    clustering and frequency counting.

    Args:
        trf_objs (list of TRModel): Tandem repeat objects to canonicalize

    Returns:
        tuple: (canonical_trf_objs, consensus_frequencies) where:
            - canonical_trf_objs (list): TRF objects with updated trf_consensus
                                        (None values removed)
            - consensus_frequencies (list): List of (count, consensus) tuples
                                           sorted by count (descending)

    Raises:
        SequenceError: If a monomer cannot generate shift variants (indicates
                      invalid non-ACGT characters or data corruption)

    Example:
        >>> # Input: Multiple representations of same repeat
        >>> trf1 = TRModel()
        >>> trf1.trf_consensus = "ACT"  # Original
        >>> trf2 = TRModel()
        >>> trf2.trf_consensus = "CTA"  # Rotation
        >>> trf3 = TRModel()
        >>> trf3.trf_consensus = "AGT"  # Reverse complement
        >>> trf4 = TRModel()
        >>> trf4.trf_consensus = "ACTACTACT"  # Multimer
        >>> objs, freqs = remove_consensus_redundancy([trf1, trf2, trf3, trf4])
        >>> all(obj.trf_consensus == "ACT" for obj in objs)
        True
        >>> freqs[0]  # (count, consensus)
        (4, 'ACT')

    Note:
        - Modifies trf_consensus attribute of input objects in-place
        - Objects with empty consensus are filtered out (set to None)
        - Multimer detection works by GC content matching and period division
        - Uses lexicographic ordering: e.g., ACT < AGT < CTA < GTA < TAC < TAG
    """
    # sort by length
    consensuses = [x.trf_consensus for x in trf_objs]
    consensuses = list(set(consensuses))
    consensuses.sort(key=lambda x: len(x))
    max_consensus_length = max(len(c) for c in consensuses) if consensuses else 2020
    length2consensuses = {}
    for i, monomer in enumerate(consensuses):
        n = len(monomer)
        length2consensuses.setdefault(n, {})
        gc = get_int_gc(monomer)
        length2consensuses[n].setdefault(gc, [])
        length2consensuses[n][gc].append(i)
    result_rules = {}
    for i, monomer in enumerate(consensuses):
        if not monomer:
            continue
        if monomer in result_rules:
            continue
        gc = get_int_gc(monomer)
        base = len(monomer)
        n = base
        variants = set(
            get_shifts_variants(monomer) + get_shifts_variants(get_revcomp(monomer))
        )
        if not variants:
            raise SequenceError(
                f"Invalid monomer sequence '{monomer}': cannot generate shift variants. "
                f"This may indicate non-ACGT characters or issues with sequence rotation/reverse-complement. "
                f"Valid monomers must contain only A, C, G, T nucleotides. "
                f"Check TRF output for data corruption or invalid consensus sequences."
            )
        lex_consensus = min(variants)
        result_rules[monomer] = lex_consensus
        while n <= max_consensus_length:
            if n in length2consensuses and gc in length2consensuses[n]:
                for k in length2consensuses[n][gc]:
                    monomer_b = consensuses[k]
                    if monomer_b in result_rules:
                        continue
                    s = n // base
                    v = set()
                    for p in range(s):
                        v.add(monomer_b[p * base : (p + 1) * base])
                    if len(v) > 1:
                        continue
                    item = v.pop()
                    if item in variants:
                        result_rules[consensuses[k]] = lex_consensus

            n += base
    variants2df = defaultdict(int)
    for i, trf_obj in enumerate(trf_objs):
        if not trf_obj.trf_consensus:
            trf_objs[i] = None
            continue
        if trf_obj.trf_consensus in result_rules:
            variants2df[result_rules[trf_obj.trf_consensus]] += 1
        else:
            variants = set(
                get_shifts_variants(trf_obj.trf_consensus) + get_shifts_variants(get_revcomp(trf_obj.trf_consensus))
            )
            lex_consensus = min(variants) if variants else trf_obj.trf_consensus
            result_rules[trf_obj.trf_consensus] = lex_consensus
            variants2df[lex_consensus] += 1
        trf_obj.trf_consensus = result_rules[trf_obj.trf_consensus]
    variants2df = sort_dictionary_by_value(variants2df, reverse=True)
    trf_objs = [x for x in trf_objs if x is not None]
    return trf_objs, variants2df


class TRFFileIO(AbstractBlockFileIO):
    """
    Parser and processor for raw TRF (Tandem Repeat Finder) output files.

    Handles TRF .dat format output, which consists of blocks starting with
    "Sequence:" lines. Each block contains tandem repeat records for one
    sequence (chromosome/scaffold). Provides filtering, overlap resolution,
    and consensus canonicalization.

    Key Features:
        - Streaming parser for memory-efficient chromosome-by-chromosome processing
        - Automatic overlap detection and merging of similar repeats
        - Redundancy removal via consensus canonicalization
        - Filtering of nested and duplicate repeats
        - Tab-delimited output format compatible with downstream tools

    Attributes:
        use_mongodb (bool): Enable MongoDB integration (inherited)
        token (str): Block start token ("Sequence:"), set in __init__

    Public Methods:
        iter_parse: Iterate over TRF file yielding filtered TRModel objects
        parse_to_file: Parse TRF file and write to tab-delimited output
        refine_old_to_file: Legacy parsing method (same as parse_to_file)

    Private Methods:
        _gen_data_line: Extract data lines from TRF block body
        _filter_obj_set: Filter overlapping/nested repeats from a sequence
        _join_overlapped: Wrapper for join_overlapped() function

    Inherited Attributes:
        data: Iterable of (head, body) tuples representing TRF blocks
        N: Number of blocks in data

    Inherited Methods:
        read_from_file, read_online, get_block_sequence, get_blocks,
        gen_block_sequences, write_to_file, iterate, clear, and others
        (see AbstractBlockFileIO for full API)

    Example:
        >>> reader = TRFFileIO()
        >>> # Parse and iterate
        >>> for trf_obj_set in reader.iter_parse("genome.fasta.2.7.7.80.10.50.500.dat"):
        ...     for trf_obj in trf_obj_set:
        ...         print(f"{trf_obj.trf_head}: {trf_obj.trf_consensus}")
        >>> # Parse to file
        >>> reader.parse_to_file("input.dat", "output.tab", project="my_genome")

    Note:
        - TRF output format: 15 space-separated fields per repeat
        - Filtering removes exact duplicates, nested repeats, and dissimilar overlaps
        - Similar overlapping repeats (cosine distance < 0.1) are merged
        - Consensus sequences are canonicalized to minimal lexicographic form
    """

    def __init__(self):
        """
        Initialize TRF file reader with "Sequence:" block delimiter.

        Sets the block start token to "Sequence:" which marks the beginning
        of each chromosome/scaffold block in TRF .dat output files.
        """
        token = "Sequence:"
        super(TRFFileIO, self).__init__(token)

    def iter_parse(self, trf_file, filter=True):
        """
        Parse TRF file and yield filtered tandem repeat objects per chromosome.

        Processes TRF .dat format file chromosome-by-chromosome, parsing each
        repeat, optionally filtering overlaps/duplicates, canonicalizing consensus
        sequences, and assigning unique IDs.

        Args:
            trf_file (str): Path to TRF output file (.dat format)
            filter (bool, optional): Apply overlap/duplicate filtering.
                                    Defaults to True.

        Yields:
            list of TRModel: Filtered and canonicalized tandem repeat objects
                            for each chromosome/scaffold block

        Example:
            >>> reader = TRFFileIO()
            >>> for chromosome_repeats in reader.iter_parse("genome.dat"):
            ...     print(f"Chromosome has {len(chromosome_repeats)} repeats")
            ...     for tr in chromosome_repeats:
            ...         print(f"  {tr.trf_consensus} at {tr.trf_l_ind}-{tr.trf_r_ind}")

        Note:
            - Assigns sequential trf_id starting from 1
            - Filtering includes duplicate removal, nested repeat removal,
              and overlap merging for similar sequences
            - All consensus sequences are canonicalized to minimal lexicographic form
        """
        trf_id = 1
        for ii, (head, body, start, next) in enumerate(self.read_online(trf_file)):
            head = head.replace("\t", " ")
            obj_set = []
            n = body.count("\n")
            for i, line in enumerate(self._gen_data_line(body)):
                trf_obj = TRModel()
                trf_obj.set_raw_trf(head, None, line)
                obj_set.append(trf_obj)
            if filter:
                # Filter object set
                trf_obj_set = self._filter_obj_set(obj_set)
                obj_set = [x for x in trf_obj_set if x]
            ### set trf_id
            for trf_obj in obj_set:
                trf_obj.trf_id = trf_id
                trf_id += 1
            obj_set, variants2df = remove_consensus_redundancy(obj_set)
            yield obj_set

    def parse_to_file(
        self, file_path, output_path, trf_id=0, project=None, verbose=True
    ):
        """
        Parse TRF file and write filtered repeats to tab-delimited output.

        Processes entire TRF .dat file, applies filtering and canonicalization,
        assigns sequential IDs, optionally sets project metadata, and writes
        results in tab-delimited format compatible with downstream analysis tools.

        Args:
            file_path (str): Path to input TRF .dat file
            output_path (str): Path to output tab-delimited file
            trf_id (int, optional): Starting ID number for repeat numbering.
                                   If 0, creates new file (mode 'w').
                                   If >0, appends to existing file (mode 'a').
                                   Defaults to 0.
            project (str, optional): Project name to add to metadata fields.
                                    Defaults to None.
            verbose (bool, optional): Enable verbose logging (currently unused).
                                     Defaults to True.

        Returns:
            int: Next available trf_id (total repeats processed + 1)

        Example:
            >>> reader = TRFFileIO()
            >>> # Process first file
            >>> next_id = reader.parse_to_file("chr1.dat", "all_repeats.tab", project="hg38")
            >>> # Append second file with continuing IDs
            >>> next_id = reader.parse_to_file("chr2.dat", "all_repeats.tab", trf_id=next_id, project="hg38")

        Note:
            - Output format: tab-delimited with header defined by TRModel.dumpable_attributes
            - Calls refine_name() to generate standardized IDs and normalize sequences
            - Creates output file if trf_id=0, appends otherwise
        """
        if trf_id == 0:
            mode = "w"
        else:
            mode = "a"

        with WiseOpener(output_path, mode) as fw:
            for trf_obj_set in self.iter_parse(file_path):
                for trf_obj in trf_obj_set:
                    trf_obj.trf_id = trf_id

                    if project:
                        trf_obj.set_project_data(project)
                    refine_name(trf_id, trf_obj)

                    fw.write(str(trf_obj))

                    trf_id += 1
        return trf_id

    def refine_old_to_file(
        self, file_path, output_path, trf_id=0, project=None, verbose=True
    ):
        """
        Legacy method for parsing TRF file to tab-delimited output.

        Identical to parse_to_file(). Retained for backward compatibility.

        Args:
            file_path (str): Path to input TRF .dat file
            output_path (str): Path to output tab-delimited file
            trf_id (int, optional): Starting ID number. Defaults to 0.
            project (str, optional): Project name. Defaults to None.
            verbose (bool, optional): Enable verbose logging. Defaults to True.

        Returns:
            int: Next available trf_id

        See Also:
            parse_to_file: Recommended method with identical functionality
        """
        if trf_id == 0:
            mode = "w"
        else:
            mode = "a"

        with WiseOpener(output_path, mode) as fw:
            for trf_obj_set in self.iter_parse(file_path):
                for trf_obj in trf_obj_set:
                    trf_obj.trf_id = trf_id

                    if project:
                        trf_obj.set_project_data(project)
                    refine_name(trf_id, trf_obj)

                    fw.write(str(trf_obj))

                    trf_id += 1
        return trf_id

    def _gen_data_line(self, data):
        """
        Extract TRF data lines from block body, skipping headers and blank lines.

        Filters out "Sequence:" headers, "Parameters:" lines, and empty lines,
        yielding only the 15-field tandem repeat data lines.

        Args:
            data (str): Raw block body text from TRF output

        Yields:
            str: Individual TRF data lines (15 space-separated fields each)

        Note:
            - Called internally during iter_parse() to extract repeat records
            - Each yielded line represents one tandem repeat annotation
        """
        for line in data.split("\n"):
            line = line.strip()
            if line.startswith("Sequence"):
                continue
            if line.startswith("Parameters"):
                continue
            if not line:
                continue
            yield line

    def _filter_obj_set(self, obj_set):
        """
        Filter and merge overlapping/nested tandem repeats from a chromosome.

        Implements complex filtering logic to handle various overlap patterns:
        - Exact duplicates: keeps higher percent match
        - Nested repeats: keeps outer (longer) repeat
        - Partial overlaps: merges if similar (via _join_overlapped)

        Filtering happens in two passes:
        1. First pass: remove exact duplicates, nested repeats, merge similar overlaps
        2. Second pass: iteratively merge remaining similar overlaps (disabled)

        Args:
            obj_set (list of TRModel): Unfiltered tandem repeat objects for one
                                      chromosome, sorted by (trf_l_ind, trf_r_ind)

        Returns:
            list of TRModel: Filtered repeat objects with duplicates/nested removed
                            and similar overlaps merged (None values removed)

        Note:
            - Modifies objects in-place via _join_overlapped()
            - Overlap patterns handled:
                a ------ (obj1)
                b ------ (exact duplicate → keep higher pmatch)

                a ------ ------ (obj1 contains obj2)
                b ---    ---    (nested → remove obj2)

                a ------
                b    ----- (partial overlap → merge if similar)

            - Second pass is currently disabled (is_overlapping always False)
              due to "suspicious results" in original implementation
        """
        # NB: I removed the overlaping part due to suspicious results.
        # Complex filter
        is_overlapping = False
        n = len(obj_set)

        obj_set.sort(key=lambda x: (x.trf_l_ind, x.trf_r_ind))
        for a in range(0, n):
            obj1 = obj_set[a]
            if not obj1:
                continue
            for b in range(a + 1, n):
                obj2 = obj_set[b]
                if not obj2:
                    continue
                # a ------
                # b ------
                if (
                    obj1.trf_l_ind == obj2.trf_l_ind
                    and obj1.trf_r_ind == obj2.trf_r_ind
                ):
                    # Check period
                    if obj1.trf_pmatch >= obj2.trf_pmatch:
                        obj_set[b] = None
                    else:
                        obj_set[a] = None
                    continue
                # a ------ ------  -------
                # b ---       ---    ---
                if (
                    obj1.trf_l_ind <= obj2.trf_l_ind
                    and obj1.trf_r_ind >= obj2.trf_r_ind
                ):
                    obj_set[b] = None
                    continue
                # a ---       ---    ---
                # b ------ ------  -------
                if (
                    obj2.trf_l_ind <= obj1.trf_l_ind
                    and obj2.trf_r_ind >= obj1.trf_r_ind
                ):
                    obj_set[a] = None
                    continue
                # a ------
                # b    -----
                if obj1.trf_r_ind > obj2.trf_l_ind and obj1.trf_r_ind < obj2.trf_r_ind:
                    if self._join_overlapped(obj1, obj2, cutoff_distance=0.1):
                        obj_set[b] = None
                    continue
                # a ------
                # b                -----
                if obj1.trf_r_ind < obj2.trf_l_ind:
                    break
                # a               ------
                # b -----
                if obj2.trf_r_ind < obj1.trf_l_ind:
                    break
        obj_set = [a for a in obj_set if not a is None]
        n = len(obj_set)

        while is_overlapping:
            is_overlapping = False

            for a in range(0, n):
                obj1 = obj_set[a]
                if not obj1:
                    continue
                for b in range(a + 1, n):
                    obj2 = obj_set[b]
                    if not obj2:
                        continue
                    # a ------
                    # b               -----
                    if obj1.trf_r_ind < obj2.trf_l_ind:
                        break
                    # a              ------
                    # b -----
                    if obj2.trf_r_ind < obj1.trf_l_ind:
                        break
                    # a ------
                    # b    -----
                    if (
                        obj1.trf_r_ind > obj2.trf_l_ind
                        and obj1.trf_r_ind < obj2.trf_r_ind
                    ):

                        overlap = float(abs(obj1.trf_r_ind - obj2.trf_l_ind))
                        min_length = min(obj1.trf_array_length, obj2.trf_array_length)
                        overlap_proc_diff = overlap * 1.0 / min_length
                        gc_dif = abs(obj1.trf_array_gc - obj2.trf_array_gc)

                        if (
                            overlap_proc_diff
                            >= 30
                            and gc_dif
                            <= 0.05
                        ):
                            is_overlapping = True
                            if self._join_overlapped(obj1, obj2):
                                obj2 = None
                        continue
                    # a ------
                    # b ------
                    if (
                        obj1.trf_l_ind == obj2.trf_l_ind
                        and obj1.trf_r_ind == obj2.trf_r_ind
                    ):
                        # Check period
                        if obj1.trf_pmatch >= obj2.trf_pmatch:
                            obj_set[b] = None
                            continue
                        else:
                            obj_set[a] = None
                            continue
                    # a ------ ------  -------
                    # b ---       ---     ---
                    if (
                        obj1.trf_l_ind <= obj2.trf_l_ind
                        and obj1.trf_r_ind >= obj2.trf_r_ind
                    ):
                        obj_set[b] = None
                        continue
                    # a ---       ---            ---
                    # b ------ ------  -------
                    if (
                        obj2.trf_l_ind <= obj1.trf_l_ind
                        and obj2.trf_r_ind >= obj1.trf_r_ind
                    ):
                        obj_set[a] = None
                        continue

            obj_set = [a for a in obj_set if not a is None]

        return obj_set

    def _join_overlapped(self, obj1, obj2, cutoff_distance=0.1):
        """
        Wrapper method for join_overlapped() function.

        Args:
            obj1 (TRModel): First tandem repeat object
            obj2 (TRModel): Second tandem repeat object
            cutoff_distance (float, optional): Maximum cosine distance for merging.
                                              Defaults to 0.1.

        Returns:
            bool: True if objects were joined, False otherwise

        See Also:
            join_overlapped: The underlying implementation function
        """
        return join_overlapped(obj1, obj2, cutoff_distance=cutoff_distance)


def sc_parse_raw_trf_folder(trf_raw_folder, output_trf_file, project=None):
    """
    Parse all TRF .dat files in a folder and write to single tab-delimited output.

    Processes all .dat files found in the specified folder, applies filtering
    and canonicalization, and writes results to a single output file with
    continuous ID numbering across all files.

    Args:
        trf_raw_folder (str): Path to folder containing TRF .dat output files
        output_trf_file (str): Path to output tab-delimited file
        project (str, optional): Project name for metadata. Defaults to None.

    Example:
        >>> sc_parse_raw_trf_folder("/path/to/trf_output/", "all_repeats.tab", project="hg38")

    Note:
        - Removes output file if it already exists (fresh write)
        - Only processes files ending with .dat extension
        - IDs are assigned sequentially starting from 1 across all files
        - Logs progress for each file processed
    """
    reader = TRFFileIO()
    trf_id = 1
    if os.path.isfile(output_trf_file):
        os.remove(output_trf_file)
    for file_path in iter_filepath_folder(trf_raw_folder, mask="dat"):
        if not file_path.endswith(".dat"):
            continue
        logger.info("Start parse file %s..." % file_path)
        trf_id = reader.parse_to_file(
            file_path, output_trf_file, trf_id=trf_id, project=project
        )


def sc_trf_to_fasta(trf_file, fasta_file):
    """
    Convert tab-delimited TRF file to FASTA format.

    Reads parsed TRF data from tab-delimited file and writes tandem repeat
    sequences in FASTA format. Uses the fasta property of TRModel objects
    which includes proper headers and sequence data.

    Args:
        trf_file (str): Path to input tab-delimited TRF file
        fasta_file (str): Path to output FASTA file

    Example:
        >>> sc_trf_to_fasta("genome_repeats.tab", "repeats.fasta")

    Note:
        - Input must be tab-delimited format from TRFFileIO.parse_to_file()
        - Output FASTA headers include repeat ID and genomic coordinates
        - Sequence data is the tandem repeat array (trf_array field)
    """
    with open(fasta_file, "w") as fw:
        for trf_obj in sc_iter_tab_file(trf_file, TRModel):
            fw.write(trf_obj.fasta)

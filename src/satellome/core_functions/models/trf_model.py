#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.06.2011
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import logging
import re

from satellome.core_functions.exceptions import SequenceError
from satellome.core_functions.models.abstract_model import AbstractModel

logger = logging.getLogger(__name__)

from satellome.core_functions.tools.parsers import (parse_chromosome_name,
                                                   parse_fasta_head,
                                                   trf_parse_head,
                                                   trf_parse_line)
from satellome.core_functions.tools.processing import get_gc_content
from satellome.core_functions.trf_embedings import create_vector, token2id, token2revtoken


def clear_sequence(sequence):
    """
    Normalize and clean DNA/RNA sequence to uppercase IUPAC nucleotides.

    Removes whitespace and non-nucleotide characters, keeping only valid
    IUPAC nucleotide codes (including ambiguity codes and gaps).

    Args:
        sequence (str): Raw DNA/RNA sequence string

    Returns:
        str: Cleaned sequence in uppercase with:
            - All whitespace removed
            - Only valid IUPAC codes retained: ACTGNUWSMKRYBDHV-
            - Uppercase format

    Example:
        >>> clear_sequence("  acgt nnnn  ")
        'ACGTNNNN'
        >>> clear_sequence("acgt123xyz")
        'ACGT'
    """
    sequence = sequence.strip().upper()
    sequence = re.sub(r"\s+", "", sequence)
    return re.sub(r"[^actgnuwsmkrybdhvACTGNUWSMKRYBDHV\-]", "", sequence)


class TRModel(AbstractModel):
    """Class for tandem repeat wrapping.

    Core Attributes:

    Identifiers:
    - project: Project name
    - trf_id: Unique tandem repeat ID
    - trf_head: Sequence header from FASTA

    Coordinates:
    - trf_l_ind (int): Left index (start position)
    - trf_r_ind (int): Right index (end position)

    Tandem Repeat Properties:
    - trf_period (int): Period size (monomer length)
    - trf_n_copy (float): Number of copies
    - trf_pmatch (float): Percent match
    - trf_pvar (float): Percent variation
    - trf_entropy (float): Sequence entropy
    - trf_array_length (int): Total array length in bp

    Sequences:
    - trf_consensus: Consensus monomer sequence
    - trf_array: Full tandem repeat array sequence

    GC Content:
    - trf_array_gc (float): GC content of array
    - trf_consensus_gc (float): GC content of consensus

    Optional Annotation Attributes:
    - trf_joined (int): 1 if created by merging overlapping TRs, 0 otherwise
    - trf_family: Family classification (e.g., "(AT)n", "tSSR_AT", "fSSR_ACG")
    - trf_ref_annotation: Annotations from GFF/RepeatMasker (pipe-separated)

    Dynamic Properties (parsed from trf_head on access):
    - trf_chr: Chromosome name (via @property)
    - trf_gi: GI identifier (via @property)

    """

    dumpable_attributes = [
        "project",
        "trf_id",
        "trf_head",
        "trf_l_ind",
        "trf_r_ind",
        "trf_period",
        "trf_n_copy",
        "trf_pmatch",
        "trf_pvar",
        "trf_entropy",
        "trf_consensus",
        "trf_array",
        "trf_array_gc",
        "trf_consensus_gc",
        "trf_array_length",
        "trf_joined",
        "trf_family",
        "trf_ref_annotation",
    ]

    # Legacy format (v1.4.2 and earlier) - 37 fields
    # Kept for backward compatibility when reading old TRF files
    legacy_dumpable_attributes = [
        "project",
        "id",  # REMOVED - was duplicate of trf_id
        "trf_id",
        "trf_type",  # REMOVED - never used
        "trf_family",
        "trf_family_prob",  # REMOVED - never used
        "trf_l_ind",
        "trf_r_ind",
        "trf_period",
        "trf_n_copy",
        "trf_pmatch",
        "trf_pvar",
        "trf_entropy",
        "trf_consensus",
        "trf_array",
        "trf_array_gc",
        "trf_consensus_gc",
        "trf_gi",  # REMOVED - now computed from trf_head
        "trf_head",
        "trf_param",  # REMOVED - always 0
        "trf_array_length",
        "trf_chr",  # REMOVED - now computed from trf_head
        "trf_joined",
        "trf_superfamily",  # REMOVED - never used
        "trf_superfamily_ref",  # REMOVED - never used
        "trf_superfamily_self",  # REMOVED - never used
        "trf_subfamily",  # REMOVED - never used
        "trf_subsubfamily",  # REMOVED - never used
        "trf_family_network",  # REMOVED - never used
        "trf_family_self",  # REMOVED - never used
        "trf_family_ref",  # REMOVED - never used
        "trf_hor",  # REMOVED - never used
        "trf_n_chrun",  # REMOVED - never used
        "trf_ref_annotation",
        "trf_bands_refgenome",  # REMOVED - never used
        "trf_repbase",  # REMOVED - never used
        "trf_strand",  # REMOVED - never used
    ]

    # Mapping from legacy fields to new fields (for fields that were kept)
    legacy_to_new_field_mapping = {
        "project": "project",
        "id": None,  # Ignored
        "trf_id": "trf_id",
        "trf_type": None,  # Ignored
        "trf_family": "trf_family",
        "trf_family_prob": None,  # Ignored
        "trf_l_ind": "trf_l_ind",
        "trf_r_ind": "trf_r_ind",
        "trf_period": "trf_period",
        "trf_n_copy": "trf_n_copy",
        "trf_pmatch": "trf_pmatch",
        "trf_pvar": "trf_pvar",
        "trf_entropy": "trf_entropy",
        "trf_consensus": "trf_consensus",
        "trf_array": "trf_array",
        "trf_array_gc": "trf_array_gc",
        "trf_consensus_gc": "trf_consensus_gc",
        "trf_gi": None,  # Ignored - computed from trf_head
        "trf_head": "trf_head",
        "trf_param": None,  # Ignored
        "trf_array_length": "trf_array_length",
        "trf_chr": None,  # Ignored - computed from trf_head
        "trf_joined": "trf_joined",
        "trf_superfamily": None,  # Ignored
        "trf_superfamily_ref": None,  # Ignored
        "trf_superfamily_self": None,  # Ignored
        "trf_subfamily": None,  # Ignored
        "trf_subsubfamily": None,  # Ignored
        "trf_family_network": None,  # Ignored
        "trf_family_self": None,  # Ignored
        "trf_family_ref": None,  # Ignored
        "trf_hor": None,  # Ignored
        "trf_n_chrun": None,  # Ignored
        "trf_ref_annotation": "trf_ref_annotation",
        "trf_bands_refgenome": None,  # Ignored
        "trf_repbase": None,  # Ignored
        "trf_strand": None,  # Ignored
    }

    int_attributes = [
        "trf_l_ind",
        "trf_r_ind",
        "trf_period",
        "trf_array_length",
        "trf_joined",
    ]

    float_attributes = [
        "trf_n_copy",
        "trf_entropy",
        "trf_pmatch",
        "trf_pvar",
        "trf_array_gc",
        "trf_consensus_gc",
    ]

    @property
    def trf_chr(self):
        """Get chromosome name parsed from trf_head."""
        return parse_chromosome_name(self.trf_head)

    @property
    def trf_gi(self):
        """Get GI identifier parsed from trf_head."""
        return parse_fasta_head(self.trf_head)[0]

    def set_project_data(self, project):
        """
        Set the project identifier for this tandem repeat.

        Args:
            project (str): Project name or identifier to associate with this TR
        """
        self.project = project

    def set_raw_trf(self, head, body, line):
        """
        Initialize TR object from raw TRF (Tandem Repeat Finder) output.

        Parses TRF output format into structured TRModel attributes, including:
        - Sequence header parsing for chromosome/GI extraction
        - TRF data line parsing for coordinates, period, copy number, etc.
        - Sequence cleaning and normalization
        - GC content computation

        Args:
            head (str): TRF sequence header line (e.g., "Sequence: chr1")
            body (str): TRF output body (currently unused, kept for compatibility)
            line (str): TRF data line with 15 space-separated fields:
                       [start, end, period, copies, consensus_size, %match, %indels,
                        score, %A, %C, %G, %T, entropy, consensus, sequence]

        Note:
            - Automatically cleans sequences (uppercase, removes invalid chars)
            - Computes GC% for both consensus and array sequences
            - Handles parsing errors gracefully (logs and sets defaults)
        """
        parsed_head = trf_parse_head(head)
        self.trf_head = parsed_head.strip() if parsed_head else "Unknown"

        (
            self.trf_l_ind,
            self.trf_r_ind,
            self.trf_period,
            self.trf_n_copy,
            self.trf_l_cons,
            self.trf_pmatch,
            self.trf_indels,
            self.trf_score,
            self.trf_n_a,
            self.trf_n_c,
            self.trf_n_g,
            self.trf_n_t,
            self.trf_entropy,
            self.trf_consensus,
            self.trf_array,
        ) = trf_parse_line(line)

        self.trf_pmatch = float(self.trf_pmatch)
        self.trf_pvar = int(100 - float(self.trf_pmatch))

        try:
            self.trf_l_ind = int(self.trf_l_ind)
        except (ValueError, TypeError) as e:
            logger.error(f"Error converting trf_l_ind to int: {e}")
            logger.error(f"Object data: {self}")
            # Set a default value or re-raise the exception
            self.trf_l_ind = 0  # or raise

        self.trf_r_ind = int(self.trf_r_ind)
        self.trf_period = int(self.trf_period)
        self.trf_n_copy = float(self.trf_n_copy)

        self.trf_consensus = clear_sequence(self.trf_consensus)
        self.trf_array = clear_sequence(self.trf_array)

        self.trf_array_gc = get_gc_content(self.trf_array)
        self.trf_consensus_gc = get_gc_content(self.trf_consensus)
        self.trf_array_length = len(self.trf_array)

    def set_form_overlap(self, obj2):
        """Init object with data from overlap with another TRFObj located right to self."""

        use_first = False
        # Convert scores to float for comparison, handle None values
        score1 = float(self.trf_score) if self.trf_score is not None else 0
        score2 = float(obj2.trf_score) if obj2.trf_score is not None else 0
        if score1 > score2:
            use_first = True

        self.trf_pmatch = float(
            (
                self.trf_pmatch * self.trf_array_length
                + obj2.trf_pmatch * obj2.trf_array_length
            )
            / (self.trf_array_length + obj2.trf_array_length)
        )

        # Calculate overlap based on coordinates, not array length
        overlap_length = self.trf_r_ind - obj2.trf_l_ind + 1
        if overlap_length > 0:
            # Repeats overlap - skip the overlapping portion from obj2
            self.trf_array = self.trf_array + obj2.trf_array[overlap_length:]
        else:
            # No overlap - just concatenate
            self.trf_array = self.trf_array + obj2.trf_array

        self.trf_r_ind = obj2.trf_r_ind
        self.trf_array_length = len(self.trf_array)

        # Verify array length matches coordinate range
        expected_length = self.trf_r_ind - self.trf_l_ind + 1
        if self.trf_array_length != expected_length:
            raise SequenceError(
                f"Array length mismatch after merging tandem repeats: "
                f"expected {expected_length} bp (from coordinates {self.trf_l_ind}-{self.trf_r_ind}), "
                f"but got {self.trf_array_length} bp. "
                f"This indicates a coordinate or sequence extraction error during TR merging. "
                f"Check overlap calculation and array concatenation logic."
            )

        if not use_first:
            self.trf_period = obj2.trf_period
            self.trf_consensus = obj2.trf_consensus
            self.trf_entropy = obj2.trf_entropy
            
        self.trf_n_copy = self.trf_array_length / self.trf_period
        self.trf_indels = None
        self.trf_score = None
        self.trf_n_a = self.trf_array.count("A")
        self.trf_n_c = self.trf_array.count("C")
        self.trf_n_g = self.trf_array.count("G")
        self.trf_n_t = self.trf_array.count("T")
        self.trf_pvar = float(100 - float(self.trf_pmatch))
        self.trf_array_gc = get_gc_content(self.trf_array)
        self.trf_consensus_gc = get_gc_content(self.trf_consensus)
        self.trf_joined = 1

    def get_vector(self):
        """Get vector representation of TRF object."""
        return create_vector(token2id, token2revtoken, self.trf_array, k=5)

    def get_string_repr(self):
        """Get string representation for TRF file."""
        return str(self)

    def get_index_repr(self):
        """Get string for index file."""
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            self.trf_id,
            self.trf_period,
            self.trf_array_length,
            self.trf_array_gc,
            self.trf_pvar,
            self.trf_gi,
            self.trf_l_ind,
            self.trf_r_ind,
            self.trf_chr,
        )

    def get_header_string(self):
        """Get header string for tsv file."""
        data = '\t'.join(self.dumpable_attributes)
        return f"#{data}\n"

    def get_numerical_repr(self):
        """Get str for Mathematica."""
        return "%s\t%s\t%.2f\n" % (
            self.trf_period,
            self.trf_array_length,
            self.trf_array_gc,
        )

    def get_fasta_repr(self, add_project=False):
        """Get array fasta representation, head - trf_id."""
        if add_project:
            return ">%s_%s\n%s\n" % (self.trf_id, self.project, self.trf_array)
        else:
            return ">%s\n%s\n" % (self.trf_id, self.trf_array)

    def get_monomer_fasta_repr(self):
        """Get monomer fasta representation, head - trf_id."""
        return ">%s\n%s\n" % (self.trf_id, self.trf_consensus)

    def get_family_repr(self):
        """Get str for family index (legacy method - simplified)."""
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            self.trf_id,
            self.trf_period,
            self.trf_array_length,
            self.trf_array_gc,
            self.trf_pvar,
            self.trf_gi,  # computed property
            self.trf_l_ind,
            self.trf_r_ind,
            self.trf_chr,  # computed property
            self.trf_family,
        )

    @property
    def fasta(self):
        return self.get_fasta_repr()

    def get_gff3_string(
        self,
        chromosome=True,
        trs_type="complex_tandem_repeat",
        probability=1000,
        tool="PySatDNA",
        prefix=None,
        properties=None,
        force_header=False,
    ):
        """
        Format tandem repeat as GFF3 (Generic Feature Format version 3) line.

        Generates a standard GFF3 line with 9 tab-separated columns for genome
        annotation visualization and analysis.

        Args:
            chromosome (bool, optional): Use chromosome name instead of GI for seqid.
                                        Defaults to True.
            trs_type (str, optional): GFF3 feature type (column 3).
                                     Defaults to "complex_tandem_repeat".
            probability (int, optional): Score value (column 6).
                                        Defaults to 1000.
            tool (str, optional): Source/tool name (column 2).
                                 Defaults to "PySatDNA".
            prefix (str, optional): Prefix to add to seqid.
                                   Defaults to None.
            properties (dict, optional): Dictionary mapping property names to
                                       attribute names for custom GFF3 attributes.
                                       Defaults to None.
            force_header (bool, optional): Force use of trf_head as seqid.
                                          Defaults to False.

        Returns:
            str: GFF3-formatted line with 9 columns:
                seqid, source, type, start, end, score, strand, phase, attributes

        Note:
            - Automatically determines strand from coordinate order
            - Swaps coordinates if trf_l_ind > trf_r_ind (reverse strand)
            - Falls back to trf_gi or trf_head if chromosome unavailable
        """
        if chromosome and self.trf_chr and self.trf_chr != "?":
            seqid = self.trf_chr
        elif self.trf_gi and self.trf_gi != "Unknown":
            seqid = self.trf_gi
        else:
            seqid = self.trf_head
        features = []
        if not properties:
            properties = {}
        for name, attr in properties.items():
            features.append("%s=%s" % (name, getattr(self, attr)))
        features = ";".join(features) if features else "."
        if prefix:
            seqid = prefix + seqid
        if self.trf_l_ind < self.trf_r_ind:
            strand = "+"
        else:
            strand = "-"
            self.trf_l_ind, self.trf_r_ind = self.trf_r_ind, self.trf_l_ind

        if force_header:
            seqid = self.trf_head
        d = (
            seqid,
            tool,
            trs_type,
            self.trf_l_ind,
            self.trf_r_ind,
            probability,
            strand,
            ".",
            features,
        )
        return "%s\n" % "\t".join(map(str, d))

    def get_bed_string(self):
        """
        Format tandem repeat as BED (Browser Extensible Data) format line.

        Generates a minimal 3-column BED format line for genome browser visualization.
        BED format uses 0-based, half-open coordinates [start, end).

        Returns:
            str: Tab-delimited BED line with 3 columns: chr, start, end

        Note:
            - Automatically swaps coordinates if trf_l_ind > trf_r_ind (reverse strand)
            - Uses trf_head as chromosome/sequence identifier
            - Strand information is computed but not included in output (minimal BED3)

        Example:
            >>> tr = TRModel()
            >>> tr.trf_head, tr.trf_l_ind, tr.trf_r_ind = "chr1", 1000, 1100
            >>> tr.get_bed_string()
            'chr1\\t1000\\t1100\\n'
        """
        if self.trf_l_ind < self.trf_r_ind:
            strand = "+"
        else:
            strand = "-"
            self.trf_l_ind, self.trf_r_ind = self.trf_r_ind, self.trf_l_ind

        seqid = self.trf_head
        d = (
            seqid,
            self.trf_l_ind,
            self.trf_r_ind,
        )
        return "%s\n" % "\t".join(map(str, d))


class NetworkSliceModel(TRModel):
    """Class for network slice data."""

    def __init__(self):
        self.dumpable_attributes = ["gid"] + self.dumpable_attributes
        self.int_attributes = ["gid"] + self.int_attributes
        super(TRModel, self).__init__()


class TRsClassificationModel(AbstractModel):
    """Model for keeping classification data."""

    dumpable_attributes = [
        "project",
        "id",
        "trf_id",
        "trf_period",
        "trf_array_length",
        "trf_array_gc",
        "trf_type",
        "trf_family",
        "trf_subfamily",
        "trf_family_prob",
        "trf_family_kmer",
        "trf_subfamily_kmer",
        "trf_family_self",
        "class_ssr",  # from cf_find_ssr
        "class_tssr",  # cf_find_ssr
        "class_sl",  # cf_find_ssr, a number of loci
        "class_good",  # cf_trs_classififcaiton_by_type
        "class_micro",  # cf_trs_classififcaiton_by_type
        "class_100bp",  # cf_trs_classififcaiton_by_type
        "class_perfect",  # cf_trs_classififcaiton_by_type
        "class_x4",  # cf_trs_classififcaiton_by_type
        "class_entropy",  # cf_trs_classififcaiton_by_type
        "class_gc",  # cf_trs_classififcaiton_by_type
        "trf_consensus",
    ]

    alt_dumpable_attributes = [
        "project",
        "id",
        "trf_id",
        "trf_period",
        "trf_array_length",
        "trf_array_gc",
        "trf_type",
        "trf_family",
        "trf_subfamily",
        "trf_family_prob",
        "class_ssr",
        "class_sl",
        "class_good",
        "class_micro",
        "class_100bp",
        "class_perfect",
        "class_x4",
        "class_entropy",
        "class_gc",
        "trf_consensus",
    ]

    int_attributes = [
        "trf_id",
        "id",
        "trf_period",
        "trf_array_length",
    ]

    float_attributes = [
        "trf_family_prob",
    ]

    def set_with_trs(self, trf_obj):
        """Set data with trf_obj."""
        for attr in self.dumpable_attributes:
            if hasattr(trf_obj, attr):
                setattr(self, attr, getattr(trf_obj, attr))

    @property
    def network_head(self):
        """Return TRs representation for network slices."""
        dumpable_attributes = [
            "trf_id",
            "trf_period",
            "trf_array_length",
            "trf_array_gc",
            "trf_type",
            "trf_family",
            "trf_subfamily",
            "trf_family_prob",
            "trf_family_kmer",
            "trf_subfamily_kmer",
            "trf_family_self",
            "class_ssr",
            "class_tssr",
            "class_sl",
            "class_good",
            "class_micro",
            "class_100bp",
            "class_perfect",
            "class_x4",
            "class_entropy",
            "class_gc",
            "trf_consensus",
        ]
        return self.get_as_string(dumpable_attributes)

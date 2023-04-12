#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 14.06.2011
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import re
from PyExp import AbstractModel
from parsers import parse_fasta_head, parse_chromosome_name
from sequence_tools import clear_sequence, get_revcomp, get_gc, check_gapped


class SequenceModel(AbstractModel):
    """Class for sequence wrapping

    Dumpable attributes:

    - seq_gi (int)
    - seq_ref
    - seq_description
    - seq_sequence
    - seq_length (int)
    - seq_gc (float)
    - seq_revcom
    - seq_gaped (int)
    - seq_chr
    - seq_head
    - seq_start_position (int)
    - seq_end_position (int)

    """

    dumpable_attributes = [
        "seq_gi",
        "seq_ref",
        "seq_description",
        "seq_sequence",
        "seq_length",
        "seq_gc",
        "seq_revcom",
        "seq_gaped",
        "seq_chr",
        "seq_head",
        "seq_start_position",
        "seq_end_position",
    ]

    int_attributes = [
        "seq_gi",
        "seq_length",
        "seq_gaped" "seq_start_position",
        "seq_end_position",
    ]

    float_attributes = [
        "seq_gc",
    ]

    def __init__(self, lower=True, protein=False):
        super(SequenceModel, self).__init__()
        self.lower = lower
        self.protein = protein

    @property
    def length(self):
        """Return a sequence length."""
        return self.seq_length

    @property
    def sequence(self):
        """Return a sequence."""
        return self.seq_sequence

    @property
    def fasta(self):
        """Return a fasta representation."""
        if not self.seq_head:
            self.seq_head = ">%s" % self.seq_ref
        return "%s\n%s\n" % (self.seq_head.strip(), self.seq_sequence.strip())

    def fastq(self, i):
        """Return a fastq representation."""
        self.seq_head = "%%m%s" % i
        return "%s\n%s\n%s\n%s\n" % (
            self.seq_head,
            self.seq_sequence.strip(),
            "+",
            "A" * self.length,
        )

    @property
    def fasta60(self):
        """Return a fasta representation."""
        if not self.seq_head:
            self.seq_head = ">%s" % self.seq_ref
        seq = self.seq_sequence.strip()
        s = []
        for i in range(0, len(seq), 60):
            s.append(seq[i : i + 60])
        seq = "\n".join(s)
        return "%s\n%s\n" % (self.seq_head.strip(), seq.strip())

    @property
    def header(self):
        return self.seq_head[1:].strip()

    def set_header(self, header):
        self.seq_head = ">%s" % header

    @property
    def contige_coverage(self):
        if "_cov_" in self.header:
            return float(self.header.split("_")[0])
        return None

    @property
    def sa_input(self):
        """Return a sequence left flanked by $."""
        return "%s$" % self.sequence

    @property
    def ncbi_fasta(self):
        """Return fasta in NCBI format with gi and ref."""
        return ">gi|%s|ref|%s|%s\n%s\n" % (
            self.seq_gi,
            self.seq_ref,
            self.seq_description,
            self.seq_sequence,
        )

    def add_sequence_revcom(self):
        """Add reverse complement."""
        self.seq_revcom = get_revcomp(self.seq_sequence)

    def set_dna_sequence(self, title, sequence, description=None, skip_clean=False):
        """Set sequence object."""
        self.seq_ref = title
        self.seq_description = description
        self.seq_sequence = re.sub("\s+", "", sequence)
        if not self.protein:
            if not skip_clean:
                self.seq_sequence = clear_sequence(sequence, lower=self.lower)
            self.seq_gc = get_gc(self.seq_sequence)
            self.seq_gaped = check_gapped(self.seq_sequence)
        self.seq_length = len(self.seq_sequence)

    def set_ncbi_sequence(self, head, sequence, skip_clean=False):
        """Set NCBI fasta sequence object."""
        self.seq_head = head
        self.set_dna_sequence(None, sequence, skip_clean=skip_clean)
        (self.seq_gi, self.seq_ref, self.seq_description) = parse_fasta_head(head)

        chr = parse_chromosome_name(head)
        if chr != "?":
            self.seq_chr = chr

    def set_gbff_sequence(self, head, sequence):
        """Set NCBI fasta sequence object."""
        self.seq_head = "\t".join(
            [
                head["gi"],
                head["ref"],
                head["description"],
            ]
        )
        self.set_dna_sequence(None, sequence)
        (self.seq_gi, self.seq_ref, self.seq_description) = (
            head["gi"],
            head["ref"],
            head["description"],
        )

        chr = parse_chromosome_name(head["description"])
        if chr != "?":
            self.seq_chr = chr

    def set_by_genbank(self, genbank_id):
        """ """
        from Bio import Entrez

        Entrez.email = "youremail@gmail.com"

        handle = Entrez.efetch(
            db="nucleotide", id=genbank_id, rettype="fasta", retmode="text"
        )
        fasta = handle.read().split("\n")
        head = fasta[0]
        sequence = "".join(fasta[1:])
        self.set_ncbi_sequence(head, sequence)



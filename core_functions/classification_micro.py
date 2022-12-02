#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 07.06.2011
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com
"""
Core functions related to classification of microsatellites tandem repeats.

- cf_separate_perfect_microsatellites(settings, project)
- cf_compute_perfect_micro_kmers(settings, project)
- cf_separate_microsatellites(settings, project)
- cf_compute_micro_kmers(settings, project)
- cf_separate_true_ssr(settings, project)
- cf_separate_fuzzy_ssr(settings, project)
- cf_separate_complex_trs(settings, project)
- cf_compute_complex_kmers(settings, project)
- cf_separate_1kb(settings, project)
- cf_separate_3kb(settings, project)
- cf_separate_10kb(settings, project)
- scf_basic_trs_classification(settings, project)

"""
import os
from collections import defaultdict

from PyExp import core_logger
from trseeker.seqio.gff_file import sc_gff3_reader
from trseeker.models.trf_model import TRModel
from trseeker.seqio.tab_file import sc_iter_tab_file
from trseeker.seqio.tr_file import save_trs_as_fasta
from trseeker.tools.ngrams_tools import compute_kmer_index_for_trf_file
from trseeker.tools.statistics import get_simple_statistics


class RepeatCountStatsModel(object):

    def __init__(self):
        self.max_length = 0
        self.min_length = 0
        self.n = 0
        self.lengths = []
        self.pmatch = []
        self.gc = []
        self.name = None

    def __str__(self):
        '''"family\tn\ttotal_length\tmin_length\tmax_length\tmean_length\tstd_length\tmin_pmatch\tmax_pmatch\tmean_pmatch\tstd_pmatch\tmin_gc\tmax_gc\tmean_gc\n"'''

        round_2 = lambda x: round(x, 2)

        self.pmatch_stats = get_simple_statistics(self.pmatch)
        self.lengths_stats = get_simple_statistics(self.lengths)
        self.gc_stats = get_simple_statistics(self.gc)
        s = "\t".join(
            map(
                str,
                [
                    self.name,
                    self.n,
                    sum(self.lengths),
                    self.min_length,
                    self.max_length,
                    round_2(self.lengths_stats["mean"]),
                    round_2(self.lengths_stats["standard_deviation"]),
                    min(self.pmatch),
                    max(self.pmatch),
                    round_2(self.pmatch_stats["mean"]),
                    round_2(self.pmatch_stats["standard_deviation"]),
                    round_2(min(self.gc)),
                    round_2(max(self.gc)),
                    round_2(self.gc_stats["mean"]),
                    round_2(self.gc_stats["standard_deviation"]),
                ],
            )
        )

        return "%s\n" % s


def _save_families_to_file(stats, family_table_file):
    """Save statistics about family frequencies into file."""
    values = list(stats.values())
    values.sort(reverse=True, key=lambda x: x.n)
    with open(family_table_file, "w") as fh:
        s = "family\tn\ttotal_length\tmin_length\tmax_length\tmean_length\tstd_length\tmin_pmatch\tmax_pmatch\tmean_pmatch\tstd_pmatch\tmin_gc\tmax_gc\tmean_gc\n"
        fh.write(s)
        for count_obj in values:
            fh.write(str(count_obj))


def _trs_separate_something(
    input_trf_file,
    output_trf_file,
    output_gff_file,
    filter_func,
    name_func,
    output_math_file=None,
    family_table_file=None,
):
    """Helper function for extracting subset of TRs from full dataset.
    Additionaly it writes gff and math files (optional).
    """
    if output_math_file:
        fh_math = open(output_math_file, "w")
    trf_objs = []
    for i, trf_obj in enumerate(sc_iter_tab_file(input_trf_file, TRModel)):
        trf_objs.append(trf_obj)
    trf_objs.sort(key=lambda x: x.trf_head)
    N = len(trf_objs)
    core_logger.info("Iterate TRs")

    stats = defaultdict(RepeatCountStatsModel)

    with open(output_trf_file, "w") as fh:
        with open(output_gff_file, "w") as fh_gff:
            selected = 0
            total_length = 0
            for i, trf_obj in enumerate(trf_objs):
                if filter_func(trf_obj):
                    selected += 1
                    total_length += trf_obj.trf_array_length
                    trf_obj.trf_family, gff_string, mathstr = name_func(trf_obj)
                    fh.write(str(trf_obj))
                    fh_gff.write(gff_string)
                    if output_math_file:
                        fh_math.write(mathstr)

                    stats[trf_obj.trf_family].n += 1
                    stats[trf_obj.trf_family].name = trf_obj.trf_family
                    stats[trf_obj.trf_family].max_length = max(
                        stats[trf_obj.trf_family].max_length, trf_obj.trf_array_length
                    )
                    stats[trf_obj.trf_family].min_length = min(
                        stats[trf_obj.trf_family].max_length, trf_obj.trf_array_length
                    )
                    stats[trf_obj.trf_family].lengths.append(trf_obj.trf_array_length)
                    stats[trf_obj.trf_family].pmatch.append(trf_obj.trf_pmatch)
                    stats[trf_obj.trf_family].gc.append(100 * trf_obj.trf_array_gc)
        if output_math_file:
            fh_gff.close()
        for key in stats:
            core_logger.info("%s\t%s\t%s" % (key, stats[key].n, stats[key].max_length))
        core_logger.info("selected %s from %s " % (selected, N))

        if family_table_file:
            _save_families_to_file(stats, family_table_file)

        return {"filtered": selected, "dataset": N, "total_length": total_length}


def scf_basic_trs_classification(settings, project):
    """
    Classify TRs into perfect microsatellites, microsatellites, true SSR,
    fuzzy SSR and complex TRs.
    """
    # cf_set_ref_trf_file(settings, project)
    cf_separate_perfect_microsatellites(settings, project)
    cf_separate_microsatellites(settings, project)
    cf_separate_true_ssr(settings, project)
    cf_separate_fuzzy_ssr(settings, project)
    cf_separate_complex_trs(settings, project)
    cf_separate_1kb(settings, project)
    cf_separate_3kb(settings, project)
    cf_separate_10kb(settings, project)
    cf_get_micro_summary_table(settings, project)

    return settings, project


def cf_separate_perfect_microsatellites(settings, project):
    """Split all TRs into perfect microsatellites and other.

    @settings:files trf_parsed_folder: folder with parsed trf_all.trf files
    @project ref_dataset: name of reference dataset
    @settings:files trf_all_file: trf_all_file (if not ref_dataset available)
    @settings:files trf_perfect_micro_file: file with perfect microsatellites TRs (monomer less than 5bp)
    @settings:files trf_work_file: file with remaining TRs
    """
    if not os.path.isdir(settings["folders"]["data_gff3"]):
        os.makedirs(settings["folders"]["data_gff3"])
    if not os.path.isdir(settings["folders"]["reports"]):
        os.makedirs(settings["folders"]["reports"])

    trf_parsed_folder = settings["folders"]["trf_parsed_folder"]
    if "ref_trf_file" in project["work_files"]:
        trf_all_file = project["work_files"]["ref_trf_file"]
    else:
        trf_all_file = settings["files"]["trf_all_file"]

    dataset = project["work_files"]["ref_assembly_name_for_trf"]

    trf_pmicro_file = settings["files"]["trf_perfect_micro_file"]
    gff_pmicro_file = settings["files"]["gff_pmicro_file"]
    family_table_file = settings["files"]["report_pmicro_file"]

    filter_func = lambda x: x.trf_period < 6 and x.trf_pmatch == 100

    def name_func(trf_obj):
        name = "(%s)n" % trf_obj.trf_consensus.upper()
        trf_obj.trf_family = name
        gff = trf_obj.get_gff3_string(
            chromosome=False,
            trs_type="perfect microsatellite",
            probability=1000,
            tool="Trevis",
            prefix=None,
            properties={
                "id": "trf_id",
                "name": "trf_family",
            },
        )
        return name, gff, None

    r = _trs_separate_something(
        trf_all_file,
        trf_pmicro_file,
        gff_pmicro_file,
        filter_func,
        name_func,
        family_table_file=family_table_file,
    )
    dataset = project["work_files"]["ref_assembly_name_for_trf"]
    n = r["filtered"]
    pgenome = round(
        100.0
        * float(r["total_length"])
        / project["work_files"]["assembly_stats"][dataset]["total_length"],
        3,
    )

    project["work_files"].setdefault("repeats", {})
    project["work_files"]["repeats"].setdefault(dataset, {})
    project["work_files"]["repeats"][dataset].setdefault("trevis", {})
    project["work_files"]["repeats"][dataset]["trevis"].setdefault("pmicro", {})
    project["work_files"]["repeats"][dataset]["trevis"]["pmicro"] = {
        "trf_file": trf_pmicro_file,
        "gff_file": gff_pmicro_file,
        "report_file": family_table_file,
        "n": n,
        "pgenome": pgenome,
    }

    return r


def cf_separate_microsatellites(settings, project):
    """Split all TRs into microsatellites and other.

    @settings:files trf_work_file: file with remaining TRs
    @settings:files trf_micro_file: file with notperfect microsatellites TRs (monomer less than 5bp)
    """
    trf_parsed_folder = settings["folders"]["trf_parsed_folder"]
    if "ref_dataset" in project:
        trf_all_file = os.path.join(
            trf_parsed_folder, project["ref_dataset"] + ".trf_all.trf"
        )
    else:
        trf_all_file = settings["files"]["trf_all_file"]

    trf_micro_file = settings["files"]["trf_micro_file"]
    gff_micro_file = settings["files"]["gff_micro_file"]
    family_table_file = settings["files"]["report_micro_file"]

    filter_func = lambda x: x.trf_period < 6

    def name_func(trf_obj):
        name = "(%s)n" % trf_obj.trf_consensus.upper()
        trf_obj.trf_family = name
        gff = trf_obj.get_gff3_string(
            chromosome=False,
            trs_type="microsatellite",
            probability=1000,
            tool="Trevis",
            prefix=None,
            properties={
                "id": "trf_id",
                "name": "trf_family",
                "pmatch": "trf_pmatch",
            },
        )
        return name, gff, None

    r = _trs_separate_something(
        trf_all_file,
        trf_micro_file,
        gff_micro_file,
        filter_func,
        name_func,
        family_table_file=family_table_file,
    )
    dataset = project["work_files"]["ref_assembly_name_for_trf"]
    n = r["filtered"]
    pgenome = round(
        100.0
        * float(r["total_length"])
        / project["work_files"]["assembly_stats"][dataset]["total_length"],
        3,
    )

    project["work_files"]["repeats"][dataset].setdefault("trevis", {})
    project["work_files"]["repeats"][dataset]["trevis"].setdefault("micro", {})
    project["work_files"]["repeats"][dataset]["trevis"]["micro"] = {
        "trf_file": trf_micro_file,
        "gff_file": gff_micro_file,
        "n": n,
        "pgenome": pgenome,
        "report_file": family_table_file,
    }

    return r


def cf_separate_true_ssr(settings, project):
    """SSR - simple sequence repeat that contains not all nucleotide.
    e.g. aaaaatataa -> SSR-AT

    @settings:files trf_work_file: file with remaining TRs
    @settings:files trf_tssr_file: file with true SSR TRs
    """
    trf_parsed_folder = settings["folders"]["trf_parsed_folder"]
    if "ref_dataset" in project:
        trf_all_file = os.path.join(
            trf_parsed_folder, project["ref_dataset"] + ".trf_all.trf"
        )
    else:
        trf_all_file = settings["files"]["trf_all_file"]

    trf_file = settings["files"]["trf_tssr_file"]
    gff_file = settings["files"]["gff_tssr_file"]
    family_table_file = settings["files"]["report_tssr_file"]

    def filter_func(x):
        array = x.trf_array.lower()
        n = float(len(array))
        a = array.count("a")
        t = array.count("t")
        c = array.count("c")
        g = array.count("g")
        if a == 0 or t == 0 or c == 0 or g == 0:
            return True
        return False

    def name_func(trf_obj):
        array = trf_obj.trf_array.lower()
        n = float(len(array))
        a = array.count("a")
        t = array.count("t")
        c = array.count("c")
        g = array.count("g")
        if (a and c and g) or (t and g and c):
            name = "tSSR_ACG"
        elif (a and c and t) or (t and g and a):
            name = "tSSR_ACT"
        elif a and g and t:
            name = "tSSR_AGT"
        elif (a and c) or (t and g):
            name = "tSSR_AC"
        elif (a and g) or (t and c):
            name = "tSSR_AG"
        elif (a and t) or (t and a):
            name = "tSSR_AT"
        elif (c and g) or (g and c):
            name = "tSSR_CG"
        elif c or g:
            name = "tSSR_C"
        elif a or t:
            name = "tSSR_A"
        else:
            raise Exception([a, c, t, g])

        trf_obj.trf_family = name
        gff = trf_obj.get_gff3_string(
            chromosome=False,
            trs_type="tSSR",
            probability=1000,
            tool="Trevis",
            prefix=None,
            properties={
                "id": "trf_id",
                "name": "trf_family",
                "pmatch": "trf_pmatch",
            },
        )
        return name, gff, None

    r = _trs_separate_something(
        trf_all_file,
        trf_file,
        gff_file,
        filter_func,
        name_func,
        family_table_file=family_table_file,
    )
    dataset = project["work_files"]["ref_assembly_name_for_trf"]
    n = r["filtered"]
    pgenome = round(
        100.0
        * float(r["total_length"])
        / project["work_files"]["assembly_stats"][dataset]["total_length"],
        3,
    )

    project["work_files"]["repeats"][dataset].setdefault("trevis", {})
    project["work_files"]["repeats"][dataset]["trevis"].setdefault("tSSR", {})
    project["work_files"]["repeats"][dataset]["trevis"]["tSSR"] = {
        "trf_file": trf_file,
        "gff_file": gff_file,
        "n": n,
        "pgenome": pgenome,
        "report_file": family_table_file,
    }
    return r


def cf_separate_fuzzy_ssr(settings, project):
    """SSR - simple sequence repeat that contains not all nucleotide.
    e.g. aaaaatataa -> SSR-AT

    @settings:files trf_work_file: file with remaining TRs
    @settings:files trf_fssr_file: file with fuzzy SSR TRs
    """
    trf_parsed_folder = settings["folders"]["trf_parsed_folder"]
    if "ref_dataset" in project:
        trf_all_file = os.path.join(
            trf_parsed_folder, project["ref_dataset"] + ".trf_all.trf"
        )
    else:
        trf_all_file = settings["files"]["trf_all_file"]

    trf_file = settings["files"]["trf_micro_file"]
    gff_file = settings["files"]["trf_fssr_file"]
    family_table_file = settings["files"]["report_fssr_file"]

    def filter_func(x):
        array = x.trf_array.lower()
        n = float(len(array))
        a = array.count("a")
        t = array.count("t")
        c = array.count("c")
        g = array.count("g")
        if a / n < 0.01 or a < 4:
            return True
        if c / n < 0.01 or c < 4:
            return True
        if g / n < 0.01 or g < 4:
            return True
        if t / n < 0.01 or t < 4:
            return True
        return False

    def name_func(trf_obj):
        array = trf_obj.trf_array.lower()
        n = float(len(array))
        a = array.count("a")
        t = array.count("t")
        c = array.count("c")
        g = array.count("g")
        if a / n < 0.01 or a < 4:
            a = False
        if c / n < 0.01 or c < 4:
            c = False
        if g / n < 0.01 or g < 4:
            g = False
        if t / n < 0.01 or t < 4:
            t = False
        if (a and c and g) or (t and g and c):
            name = "fSSR_ACG"
        elif (a and c and t) or (t and g and a):
            name = "fSSR_ACT"
        elif a and g and t:
            name = "fSSR_AGT"
        elif (a and c) or (t and g):
            name = "fSSR_AC"
        elif (a and g) or (t and c):
            name = "fSSR_AG"
        elif (a and t) or (t and a):
            name = "fSSR_AT"
        elif (c and g) or (g and c):
            name = "fSSR_CG"
        elif c or g:
            name = "fSSR_C"
        elif a or t:
            name = "fSSR_A"
        else:
            raise Exception([a, c, t, g])

        trf_obj.trf_family = name
        gff = trf_obj.get_gff3_string(
            chromosome=False,
            trs_type="fSSR",
            probability=1000,
            tool="Trevis",
            prefix=None,
            properties={
                "id": "trf_id",
                "name": "trf_family",
                "pmatch": "trf_pmatch",
            },
        )
        return name, gff, None

    r = _trs_separate_something(
        trf_all_file,
        trf_file,
        gff_file,
        filter_func,
        name_func,
        family_table_file=family_table_file,
    )
    dataset = project["work_files"]["ref_assembly_name_for_trf"]
    n = r["filtered"]
    pgenome = round(
        100.0
        * float(r["total_length"])
        / project["work_files"]["assembly_stats"][dataset]["total_length"],
        3,
    )

    project["work_files"]["repeats"][dataset].setdefault("trevis", {})
    project["work_files"]["repeats"][dataset]["trevis"].setdefault("fSSR", {})
    project["work_files"]["repeats"][dataset]["trevis"]["fSSR"] = {
        "trf_file": trf_file,
        "gff_file": gff_file,
        "n": n,
        "pgenome": pgenome,
        "report_file": family_table_file,
    }
    return r


def cf_separate_complex_trs(settings, project):
    """Separate complex tandem repeats types.

    @settings:files trf_work_file: file with remaining TRs
    @settings:files trf_complex_file: file with complex TRs
    """
    trf_parsed_folder = settings["folders"]["trf_parsed_folder"]
    if "ref_dataset" in project:
        trf_all_file = os.path.join(
            trf_parsed_folder, project["ref_dataset"] + ".trf_all.trf"
        )
    else:
        trf_all_file = settings["files"]["trf_all_file"]

    if not os.path.isdir(settings["folders"]["mathematica"]):
        os.makedirs(settings["folders"]["mathematica"])

    trf_complex_file = (
        settings["files"]["trf_complex_file"]
        % project["work_files"]["ref_assembly_name_for_trf"]
    )
    gff_complex_file = (
        settings["files"]["gff_complex_file"]
        % project["work_files"]["ref_assembly_name_for_trf"]
    )
    clouds_complex_file = (
        settings["files"]["clouds_complex_file"]
        % project["work_files"]["ref_assembly_name_for_trf"]
    )
    clouds_png_complex_file = (
        settings["files"]["clouds_png_complex_file"]
        % project["work_files"]["ref_assembly_name_for_trf"]
    )

    filter_func = (
        lambda x: len(x.trf_consensus) > 4
        and x.trf_array_length > 100
        and x.trf_pmatch < 100.0
        and x.trf_array_gc > 0.2
        and x.trf_array_gc < 0.8
        and x.trf_n_copy > 4
        and x.trf_entropy > 1.82
    )

    def name_func(trf_obj):
        gff = trf_obj.get_gff3_string(
            chromosome=False,
            trs_type="complex TRs",
            probability=1000,
            tool="Trevis",
            prefix=None,
            properties={
                "id": "trf_id",
                "name": "trf_family",
                "pmatch": "trf_pmatch",
                "gc": "trf_array_gc",
                "length": "trf_array_length",
                "period": "trf_period",
            },
        )
        d = (
            trf_obj.trf_period,
            trf_obj.trf_array_gc,
            trf_obj.trf_pmatch,
            trf_obj.trf_family,
        )
        mathstr = "%s\n" % "\t".join(map(str, d))
        return trf_obj.trf_family, gff, mathstr

    r = _trs_separate_something(
        trf_all_file,
        trf_complex_file,
        gff_complex_file,
        filter_func,
        name_func,
        output_math_file=clouds_complex_file,
    )
    dataset = project["work_files"]["ref_assembly_name_for_trf"]
    n = r["filtered"]
    pgenome = round(
        100.0
        * float(r["total_length"])
        / project["work_files"]["assembly_stats"][dataset]["total_length"],
        3,
    )

    project["work_files"]["repeats"][dataset].setdefault("trevis", {})
    project["work_files"]["repeats"][dataset]["trevis"].setdefault("compex", {})
    project["work_files"]["repeats"][dataset]["trevis"]["compex"] = {
        "trf_file": trf_complex_file,
        "gff_file": gff_complex_file,
        "n": n,
        "pgenome": pgenome,
        "clouds_complex_file": clouds_complex_file,
    }
    return r


def cf_separate_1kb(settings, project):
    """Split all TRs by length greater 1kb."""

    trf_parsed_folder = settings["folders"]["trf_parsed_folder"]
    if "ref_trf_file" in project["work_files"]:
        trf_all_file = project["work_files"]["ref_trf_file"]
    else:
        trf_all_file = settings["files"]["trf_all_file"]

    trf_file = settings["files"]["trf_1k_file"]
    gff_file = settings["files"]["gff_1k_file"]
    fasta_file = settings["files"]["trf_1k_fasta_file"]

    filter_func = lambda x: x.trf_array_length > 1000

    def name_func(trf_obj):
        gff = trf_obj.get_gff3_string(
            chromosome=False,
            trs_type="complex TRs",
            probability=1000,
            tool="Trevis",
            prefix=None,
            properties={
                "id": "trf_id",
                "name": "trf_family",
                "pmatch": "trf_pmatch",
                "gc": "trf_array_gc",
                "length": "trf_array_length",
                "period": "trf_period",
            },
        )
        return trf_obj.trf_family, gff, None

    r = _trs_separate_something(
        trf_all_file, trf_file, gff_file, filter_func, name_func
    )

    save_trs_as_fasta(trf_file, fasta_file, project)

    dataset = project["work_files"]["ref_assembly_name_for_trf"]
    n = r["filtered"]
    pgenome = round(
        100.0
        * float(r["total_length"])
        / project["work_files"]["assembly_stats"][dataset]["total_length"],
        3,
    )

    project["work_files"]["repeats"][dataset].setdefault("trevis", {})
    project["work_files"]["repeats"][dataset]["trevis"].setdefault("1kb", {})
    project["work_files"]["repeats"][dataset]["trevis"]["1kb"] = {
        "trf_file": trf_file,
        "gff_file": gff_file,
        "fasta_file": fasta_file,
        "n": n,
        "pgenome": pgenome,
    }

    return r


def cf_separate_3kb(settings, project):
    """Split all TRs by length greater 3kb."""

    trf_parsed_folder = settings["folders"]["trf_parsed_folder"]
    if "ref_trf_file" in project["work_files"]:
        trf_all_file = project["work_files"]["ref_trf_file"]
    else:
        trf_all_file = settings["files"]["trf_all_file"]

    trf_file = settings["files"]["trf_3k_file"]
    gff_file = settings["files"]["gff_3k_file"]
    fasta_file = settings["files"]["trf_3k_fasta_file"]

    filter_func = lambda x: x.trf_array_length > 3000

    def name_func(trf_obj):
        gff = trf_obj.get_gff3_string(
            chromosome=False,
            trs_type="complex TRs",
            probability=1000,
            tool="Trevis",
            prefix=None,
            properties={
                "id": "trf_id",
                "name": "trf_family",
                "pmatch": "trf_pmatch",
                "gc": "trf_array_gc",
                "length": "trf_array_length",
                "period": "trf_period",
            },
        )
        return trf_obj.trf_family, gff, None

    r = _trs_separate_something(
        trf_all_file, trf_file, gff_file, filter_func, name_func
    )

    save_trs_as_fasta(trf_file, fasta_file, project)

    dataset = project["work_files"]["ref_assembly_name_for_trf"]
    n = r["filtered"]
    pgenome = round(
        100.0
        * float(r["total_length"])
        / project["work_files"]["assembly_stats"][dataset]["total_length"],
        3,
    )

    project["work_files"]["repeats"][dataset].setdefault("trevis", {})
    project["work_files"]["repeats"][dataset]["trevis"].setdefault("3kb", {})
    project["work_files"]["repeats"][dataset]["trevis"]["3kb"] = {
        "trf_file": trf_file,
        "gff_file": gff_file,
        "fasta_file": fasta_file,
        "n": n,
        "pgenome": pgenome,
    }

    return r


def cf_separate_10kb(settings, project):
    """Split all TRs by length greater 10kb."""

    trf_parsed_folder = settings["folders"]["trf_parsed_folder"]
    if "ref_trf_file" in project["work_files"]:
        trf_all_file = project["work_files"]["ref_trf_file"]
    else:
        trf_all_file = settings["files"]["trf_all_file"]

    trf_file = settings["files"]["trf_10k_file"]
    gff_file = settings["files"]["gff_10k_file"]
    fasta_file = settings["files"]["trf_10k_fasta_file"]

    filter_func = lambda x: x.trf_array_length > 10000

    def name_func(trf_obj):
        gff = trf_obj.get_gff3_string(
            chromosome=False,
            trs_type="complex TRs",
            probability=1000,
            tool="Trevis",
            prefix=None,
            properties={
                "id": "trf_id",
                "name": "trf_family",
                "pmatch": "trf_pmatch",
                "gc": "trf_array_gc",
                "length": "trf_array_length",
                "period": "trf_period",
            },
        )
        return trf_obj.trf_family, gff, None

    r = _trs_separate_something(
        trf_all_file, trf_file, gff_file, filter_func, name_func
    )

    save_trs_as_fasta(trf_file, fasta_file, project)

    dataset = project["work_files"]["ref_assembly_name_for_trf"]
    n = r["filtered"]
    pgenome = round(
        100.0
        * float(r["total_length"])
        / project["work_files"]["assembly_stats"][dataset]["total_length"],
        3,
    )

    project["work_files"]["repeats"][dataset].setdefault("trevis", {})
    project["work_files"]["repeats"][dataset]["trevis"].setdefault("10kb", {})
    project["work_files"]["repeats"][dataset]["trevis"]["10kb"] = {
        "trf_file": trf_file,
        "gff_file": gff_file,
        "fasta_file": fasta_file,
        "n": n,
        "pgenome": pgenome,
    }

    return r


def cf_compute_micro_kmers(settings, project):
    """
    Compute kmer for noperfect microsatellites

    @settings:files trf_micro_file: file with noperfect microsatellites TRs (monomer less than 5bp)
    @settings:files trf_micro_kmers_file: file with noperfect microsatellites kmers
    """
    trf_micro_file = settings["files"]["trf_micro_file"]
    trf_micro_kmers_file = settings["files"]["trf_micro_kmers_file"]
    index = compute_kmer_index_for_trf_file(
        trf_micro_file,
        trf_micro_kmers_file,
        k=23,
        max_complexity=None,
        min_complexity=None,
    )


def cf_compute_perfect_micro_kmers(settings, project):
    """
    Compute kmer for perfect microsatellites

    @settings:files trf_perfect_micro_file: file with perfect microsatellites TRs (monomer less than 5bp)
    @settings:files trf_pmicro_kmers_file: file with perfect microsatellites kmers
    """
    trf_micro_file = settings["files"]["trf_perfect_micro_file"]
    trf_micro_kmers_file = settings["files"]["trf_pmicro_kmers_file"]
    index = compute_kmer_index_for_trf_file(
        trf_pmicro_file,
        trf_micro_kmers_file,
        k=23,
        max_complexity=None,
        min_complexity=None,
    )


def cf_get_micro_summary_table(settings, project):
    """ """
    dataset = project["work_files"]["ref_assembly_name_for_trf"]

    input_gff = project["work_files"]["repeats"][dataset]["trevis"]["micro"]["gff_file"]
    output_tsv = os.path.join(
        settings["reports"],
        project["pid"],
        "microsatellites.summary.tsv",
    )

    report_folder = os.path.join(
        settings["reports"], project["pid"]
    )
    if not os.path.isdir(report_folder):
        os.makedirs(report_folder)

    micro_s = defaultdict(int)
    pmicro_s = defaultdict(int)
    nmicro_s = defaultdict(int)

    lengths = defaultdict(int)

    keys = set()

    for gff_obj in sc_gff3_reader(input_gff):

        name = gff_obj.attributes["name"]

        keys.add(name)

        pmatch = float(gff_obj.attributes["pmatch"])
        if pmatch == 100.0:
            pmicro_s[name] += 1
        else:
            nmicro_s[name] += 1
        micro_s[name] += 1
        lengths[name] += abs(int(gff_obj.end) - int(gff_obj.start))

    data = []
    for name in keys:
        s = (
            name,
            micro_s[name],
            lengths[name],
            round(100.0 * nmicro_s[name] / micro_s[name], 2),
            round(100.0 * pmicro_s[name] / micro_s[name], 2),
        )
        data.append(s)
    data.sort(reverse=True, key=lambda x: x[1])

    with open(output_tsv, "w") as fh:
        s = "#Name\t#\tLength (bp)\t%unperfect\t%perfect\n"
        fh.write(s)
        for d in data:
            s = "%s\t%s\t%s\t%s\t%s\n" % d
            print(s.strip())
            fh.write(s)

#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 03.03.2023
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse
import re

from satellome.core_functions.io.fasta_file import sc_iter_fasta_brute

# https://convertingcolors.com/hex-color-969696.html?search=Hex(969696)
# spectral
# http://www.di.fc.ul.pt/~jpn/r/GraphicalTools/colorPalette.html
COLOR_PALETTE_HR = {
    "dpink": "#9E0142",
    "dred": "#D53E4F",
    "orange": "#F46D43",
    "gold": "#FDAE61",
    "yellow": "#FEE08B",
    "lyellow": "#FFFFBF",
    "lgreen": "#E6F598",
    "green": "#ABDDA4",
    "dgreen": "#66C2A5",
    "blue": "#3288BD",
    "dblue": "#5E4FA2",
    "grey": "#969696",
}

# spectral
COLOR_PALETTE = {
    "w600": "#9E0142",
    "w500": "#D53E4F",
    "w400": "#F46D43",
    "w300": "#FDAE61",
    "w200": "#FEE08B",
    "w100": "#FFFFBF",
    "c100": "#E6F598",
    "c200": "#ABDDA4",
    "c300": "#66C2A5",
    "c400": "#3288BD",
    "c500": "#5E4FA2",
}

# https://rdrr.io/cran/ggsci/man/pal_jco.html
# https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html
COLOR_PALETTE_JCO = {
    "blue": "#0073C2FF",
    "yellow": "#EFC000FF",
    "grey": "#868686FF",
    "red": "#CD534CFF",
    "lblue": "#7AA6DCFF",
    "dblue": "#003C67FF",
    "poop": "#8F7700FF",
    "dgrey": "#3B3B3BFF",
    "cherry": "#A73030FF",
    "rain": "#4A6990FF",
}


def main(args):
    """Output file format:
    tr_class name color
    """

    header2name = {}
    header2size = {}
    for i, (header, seq) in enumerate(sc_iter_fasta_brute(args.fasta)):
        if len(seq) < args.cutoff:
            continue
        print(header)
        name = "chr%s" % i
        suggest = re.findall(r"chromosome (\d+)", header)
        if suggest:
            name = "chr%s" % suggest[0]
            if args.manual:
                name = input(f"Set name ({name})?") or name
        else:
            name = input(f"Set name ({name})?") or name
        header2name[header] = name
        header2size[header] = len(seq)
        print(f"Set name {name} with length ({len(seq)} bp)")

    with open(args.output, "w") as fh:
        for tr_class, (name, color) in header2name.items():
            fh.write("%s\t%s\t%s\n" % (tr_class, name, color))


def get_args():
    parser = argparse.ArgumentParser(
        description="Set TRs names and colors and create config file"
    )
    parser.add_argument("-i", "--trf", help="TRF file", required=True)
    parser.add_argument("-o", "--output", help="Output file", required=True)
    parser.add_argument("--manual", help="Ask each TRs", default=False)
    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = get_args()
    main(args)

#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 10.03.2019
# @author: Aleksey Komissarov
# @contact: ad3002@gmail.com

import argparse
from collections import defaultdict

from satellome.core_functions.models.gff3_model import sc_gff3_reader

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Print micro satellite stats.")
    parser.add_argument("-i", "--input", help="Input gff3 file", required=True)
    parser.add_argument("-o", "--output", help="Output tsv file", required=True)
    args = vars(parser.parse_args())

    input_gff = args["input"]
    output_tsv = args["output"]

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

    with open(output_tsv, "w") as fh:
        for name in keys:
            s = "%s\t%s\t%s\t%s\t%s\n" % (
                name,
                micro_s[name],
                lengths[name],
                round(100.0 * nmicro_s[name] / micro_s[name], 2),
                round(100.0 * pmicro_s[name] / micro_s[name], 2),
            )
            print(s.strip())
            fh.write(s)

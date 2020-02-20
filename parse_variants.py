#!/usr/bin/env python

# Import packages
from itertools import compress
import pybedtools as py
import numpy as np

flpath = "/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Variant_Calling/Run_Eliprotocol_12_02_18/All_SNP.vcf"

vcf = py.BedTool(flpath)

## Print output headers ##


## Some functions ##

pm = vcf[0].fields
print(pm)
len(pm)


def getRatioDepth(GF):
    if GF[0] == "./.":
        ratio = np.nan
        dp = 0
    else:
        rf = int(GF[1].split(",")[0])
        alt = int(GF[1].split(",")[1])
        dp = rf+alt

        if dp == 0:
            ratio = np.nan
        else:
            ratio = round(rf / dp, 1)

    return(ratio, dp)

## Main Loop ##


for pm in vcf:

    ref = pm.fields[3]
    alt = pm.fields[4]
    pos = pm.start
    chrom = pm.chrom

    in12b = pm.fields[9].split(":")
    in10g = pm.fields[10].split(":")
    inA7 = pm.fields[11].split(":")
    inC2 = pm.fields[12].split(":")
    inE5 = pm.fields[13].split(":")

    r1, d1 = getRatioDepth(in12b)
    r2, d2 = getRatioDepth(in10g)
    r3, d3 = getRatioDepth(inA7)
    r4, d4 = getRatioDepth(inC2)
    r5, d5 = getRatioDepth(inE5)

    contrasts = [abs(r1-r2), abs(r1-r3), abs(r1-r4), abs(r1-r5),
                 abs(r2-r3), abs(r2-r4), abs(r2-r5),
                 abs(r3-r4), abs(r3-r5),
                 abs(r4-r5)]

    names = ["1.2B_10G", "1.2B_A7", "1.2B_C2", "1.2B_E5",
             "10G_A7", "10G_C2", "10G_E5",
             "A7_C2", "A7_E5",
             "C2_E5"]

    depths = [d1, d2, d3, d4, d5]
    freqs = [r1, r2, r3, r4, r5]

    if any(i > 0.8 for i in contrasts) and all(d > 20 for d in depths):
        mask = [x > 0.8 for x in contrasts]
        output = [ref, alt, chrom, pos,
                  r1, r2, r3, r4, r5,
                  d1, d2, d3, d4, d5]

        strout = [str(x) for x in output]

        print(list(compress(names, mask)), "\t", "\t".join(strout))

    # if (any(i > 0.6 for i in [contrasts[y] for y in [1, 2, 4, 9, 11, 13]])
    #         and all(d > 20 for d in [depths[y] for y in[0, 2, 3, 5]])):

    #     mask = [x > 0.6 for x in [contrasts[y] for y in [1, 2, 4, 9, 11, 13]]]
    #     output = [ref, alt, chrom, pos,
    #               r1, r2, r3, r4, r5, r6,
    #               d1, d2, d3, d4, d5, d6]

    #     strout = [str(x) for x in output]

    #     print(list(compress(n, mask)), "\t", "\t".join(strout))

        # if any(f > 0.8 for f in freqs) and all(d > 100 for d in depths):
        #     output = [ref, alt, chrom, pos,
        #               r1, r2, r3, r4, r5, r6,
        #               d1, d2, d3, d4, d5, d6]

        #     strout = [str(x) for x in output]

        #     print("\t".join(strout))

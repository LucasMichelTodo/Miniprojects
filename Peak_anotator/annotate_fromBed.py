#!/usr/bin/env python

# Import packages
import sys
import pybedtools as py
import subprocess as sp

# Load directly edited gff with 5' and 3' regions of each gene.

ref = py.BedTool(
    "/home/lucas/ISGlobal/Gen_Referencies/anotation_PlasmoDB-41_sorted.gff")


# Anotating main function.

def annotate_from_bed(bed_file):

    bed = py.BedTool(bed_file)
    bedlen = len(bed[0].fields)

    with open(bed_file.replace(".bed", "_annotated.bed"), "w+") as outfile:

        # Intersect each gene in the gff with the peaks in the bed file.
        intersect = bed.intersect(ref, wao=True)

        for line in intersect:

            peak = line[bedlen-1]

            if line[bedlen] != ".":
                info = line[bedlen+8].split(";")
                gid = info[0].replace("ID=", "")
                if line[bedlen+2] == "gene":
                    anot = [x for x in info if "description=" in x][0].replace(
                        "description=", "")
                else:
                    anot = [x for x in info if "Parent" in x][0]
                out = line.fields[0:bedlen] + [peak, gid, anot]
            else:
                out = line.fields[0:bedlen] + ["no-hit", "no-hit"]

            outfile.write("\t".join(out)+"\n")
            #print("\t".join(out))


if __name__ == "__main__":
    filenames = sys.argv[1:]
    for i in filenames:
        print("Anotating: ", i)
        annotate_from_bed(i)

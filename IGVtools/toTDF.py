import subprocess as sp
import os
import sys


def calltoTDF(infile):

    outfile = infile+".tdf"
    cmd = ("~/Programs/IGV_2.4.10/IGVTools/igvtools toTDF "
           f"{infile} {outfile} "
           "~/Programs/IGV_2.4.10/Custom_Genomes/PlasmoDB-41_Pfalciparum3D7.genome")

    sp.call(cmd, shell=True)


if __name__ == "__main__":
    file_list = sys.argv[1:]
    for file in file_list:
        calltoTDF(file)

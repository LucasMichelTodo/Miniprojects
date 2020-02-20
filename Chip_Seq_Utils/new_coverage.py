import os
import subprocess as sp
import pandas as pd
import numpy as np


def remove_duplicates(bam):

    i = bam
    o = bam.replace(".bam", "_noDup.bam")
    m = bam.replace(".bam", "_metrics.txt")

    cmd = ("java -jar /home/lucas/Programs/picard.jar MarkDuplicates "
           "REMOVE_DUPLICATES=true I={} O={} M={}") .format(i, o, m)

    sp.call(cmd, shell=True)


def index_bam(bam):
    cmd = "samtools index {}" .format(bam)
    sp.call(cmd, shell=True)


def get_scaling(bam):
    nreads = {}
    cmd = "samtools view -c -f 1 -F 12 {}" .format(bam)
    nreads = int(sp.check_output(cmd, shell=True))
    scaling = nreads/2000000
    return(scaling)


def coverage(bam, binsize, process):

    out = bam.replace(".bam", "_cov.bdg")
    params = [bam, out, binsize, process]

    cmd = ("bamCoverage --bam {} --outFileName {} "
           "--binSize {} --numberOfProcessors {} "
           "--extendReads --outFileFormat bedgraph") .format(*params)

    print(cmd)
    sp.call(cmd, shell=True)


def normCoverage(bam, binsize, process):

    out = bam.replace(".bam", "_cov.bdg")
    params = [bam, out, binsize, process]

    cmd = ("bamCoverage --bam {} --outFileName {} "
           "--binSize {} --numberOfProcessors {} "
           "--extendReads --outFileFormat bedgraph "
           "--normalizeUsingRPKM") .format(*params)

    print(cmd)
    sp.call(cmd, shell=True)


def normalizeRPKM(bed, scaling):
    """ Divide coverage in each position by number of reads. 
    Requires scaling (from get_scaling), 
    the number of reads in the alignment. """

    with open(bed, "r+") as infile:
        with open(bed.replace(".bdg", "_RPKM.bdg"), "w+") as outfile:
            for line in infile:

                linelist = line.strip().split("\t")

                raw_val = int(linelist[3])
                kb = (int(linelist[2]) - int(linelist[1]))/1000

                rpkm = raw_val/(scaling*kb)
                outfile.write("\t".join(linelist[0:3])+"\t"+str(rpkm)+"\n")


def divideByInput(chip_bed, input_bed, pseudo=0.1, log=True):
    """ Divide chip_bed signal by chip_input signal.
        Add a pseudocount to both chip and input.
        Optionally log the result (log2, default: True). """

    chip = pd.read_csv(chip_bed, sep="\t", header=None)
    ctl = pd.read_csv(input_bed, sep="\t", header=None)

    chip.iloc[:, 3] = (chip.iloc[:, 3]+pseudo)/(ctl.iloc[:, 3]+pseudo)

    if log:
        chip.iloc[:, 3] = np.log2(chip.iloc[:, 3])

    chip.to_csv(path_or_buf=chip_bed.replace(
        ".bdg", "_normIN.bdg"), sep="\t", index=False, header=None)


def substractInput(chip_bed, input_bed,):

    treat = pd.read_csv(chip_bed, sep="\t", header=None)
    ctl = pd.read_csv(input_bed, sep="\t", header=None)

    treat.iloc[:, 3] = treat.iloc[:, 3]-ctl.iloc[:, 3]
    treat.to_csv(path_or_buf=chip_bed.replace(
        ".bdg", "_subsIN.bdg"), sep="\t", index=False, header=None)


def standarizeBed(bed):

    df = pd.read_csv(bed, sep="\t", header=None)
    df.iloc[:, 3:] = df.iloc[:, 3:].apply(
        lambda x: (x-x.mean())/x.std(), axis=0)
    df.to_csv(path_or_buf=bed.replace(
        ".bdg", "_STZ.bdg"), sep="\t", index=False, header=None)


def toTDF(infile):
    """ Input a .bdg file (bed) and create an IGV readable .tdf using igvtools. """

    outfile = infile+".tdf"
    cmd = ("~/Programs/IGV_2.4.10/IGVTools/igvtools toTDF "
           f"{infile} {outfile} "
           "~/Programs/IGV_2.4.10/Custom_Genomes/PlasmoDB-41_Pfalciparum3D7.genome")

    sp.call(cmd, shell=True)

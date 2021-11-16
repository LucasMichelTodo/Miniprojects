#### Functions ####
import subprocess
import os
import itertools
import pybedtools as pb
import subprocess as sp
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import pandas as pd
import itertools

## Functions
def get_RPKMs(bam, bs, smooth, out):

    cmd = (
        'bamCoverage -b {} '
        '--outFileFormat bedgraph '
        '--normalizeUsing RPKM '
        '-p 8 '
        '-bs {} '
        '--smoothLength {} '
        '-o {}'
    )

    subprocess.call(cmd .format(bam, bs, smooth, out), shell=True)


def get_RPKMs_normInput(bam_IP, bam_in, bs, smooth, out):

    cmd = (
        'bamCompare -b1 {} -b2 {} '
        '--outFileFormat bedgraph '
        '--scaleFactorsMethod None '
        '--normalizeUsing RPKM '
        '-p 8 '
        '-bs {} '
        '--smoothLength {} '
        '-o {}'
    )

    subprocess.call(cmd .format(bam_IP, bam_in, bs, smooth, out), shell=True)

def subtract_bams(bam1, bam2, outname):

    cmd = ('bamCompare -b1 {} -b2 {} '
           '--operation subtract '
           '-bs 100 '
           '-p 8 '
           '--outFileFormat bedgraph '
           '--scaleFactorsMethod None '
           '--normalizeUsing RPKM '
           '--smoothLength 500 '
           '-o {}')
    sp.call(cmd .format(bam1, bam2, outname), shell=True)


def getName(filename):
    if '/' in filename:
        filename = filename.rsplit("/", 1)[1] #Remove path
    out = filename.split("_", 1)[0] #Remove added names
    return(out)


def getDepth(peaksfile):
    dlist = []
    with open(peaksfile) as f:
        for line in f:
            if line.startswith("# total fragments "):
                d = line.split(":")[1].strip()
                dlist.append(int(d))
    depth = int(round(min(dlist) / 1000000))
    return(str(depth))


def macs2callpeak(t, c, params):

    # Make sure we are using apropiate MACS2 version (2.1.2)
    print("You are using MACS version:")
    cmd = "macs2 --version"
    subprocess.call(cmd, shell=True)
    print("\n")

    cmd = ("macs2 callpeak -t {} -c {} ") .format(t, c) + params

    subprocess.call(cmd, shell=True)


def macs2DifPeaks(t1, c1, t2, c2, g, l, c, outdir):

    # Make sure we are using apropiate MACS2 version (2.1.2)
    print("You are using MACS version:")
    cmd = "macs2 --version"
    subprocess.call(cmd, shell=True)
    print("\n")

    t1pile = getName(t1) + "_me_Macspeaks_treat_pileup.bdg"
    c1pile = getName(t1) + "_me_Macspeaks_control_lambda.bdg"

    t2pile = getName(t2) + "_me_Macspeaks_treat_pileup.bdg"
    c2pile = getName(t2) + "_me_Macspeaks_control_lambda.bdg"

    peaks1 = getName(t1) + "_me_Macspeaks_peaks.xls"
    peaks2 = getName(t2) + "_me_Macspeaks_peaks.xls"

    d1 = getDepth(peaks1)
    d2 = getDepth(peaks2)

    prefix = getName(t1)+"_vs_"+getName(t2)+"_g{}_l{}_c{}" .format(g,l,c)

    cmd = ("macs2 bdgdiff "
           "--t1 {} --c1 {} "
           "--t2 {} --c2 {} "
           "--d1 {} --d2 {} "
           "--outdir {} "
           "--o-prefix {} "
           "-g {} -l {} --cutoff {}") .format(t1pile, c1pile,
                                              t2pile, c2pile,
                                              d1, d2,
                                              outdir,
                                              prefix,
                                              g, l, c)

    print(cmd)
    subprocess.call(cmd, shell=True)

def manormDifPeaks(peaks1, peaks2, reads1, reads2, params):

    #Check which MaNorm version we are using:
    print("\n")
    print("You are using MANorm version:")
    cmd = "manorm --version"
    subprocess.call(cmd, shell=True)
    print("\n")

    cmd = ("manorm --p1 {} --p2 {} --r1 {} --r2 {} ") .format(peaks1, peaks2,
                                                              reads1, reads2)
    cmd = cmd + params
    subprocess.call(cmd, shell=True)

def create_log2_tracks(treat, ctl):
    outfile = treat.replace('.bam', '_log2.bdg')
    cmd = 'bamCompare -b1 {} -b2 {} -o {} -of bedgraph -p8' .format(treat, ctl, outfile)
    subprocess.call(cmd, shell=True)

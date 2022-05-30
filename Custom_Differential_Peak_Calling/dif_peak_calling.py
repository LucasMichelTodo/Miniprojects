import time
import os
import pybedtools as pb
import numpy as np
import subprocess as sp
from tqdm import tqdm
from scipy.stats import skewnorm
import matplotlib.pyplot as plt
import pandas as pd
from multiprocessing import Pool
from itertools import repeat

## Functions

def plot_skewnorm(cov_vector, params, plotname):

    # Plot the histogram.
    plt.hist(cov_vector, bins=25, density=True, alpha=0.6, color='g')

    # Plot the PDF.
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = skewnorm.pdf(x, *params)
    plt.plot(x, p, 'k', linewidth=2)
    plt.savefig(plotname)
    plt.clf()

def get_probs(coverage, params):
    return(skewnorm.cdf(coverage, *params))

def compare_cov(feat, skewnorm_params1, skewnorm_params2):
    c1 = float(feat.fields[3])
    c2 = float(feat.fields[4])
    prob1 = skewnorm.cdf(c1, *skewnorm_params1)
    prob2 = skewnorm.cdf(c2, *skewnorm_params2)
    probdif = prob1 - prob2
    return(probdif)

def get_differential_peaks(
        peakfile1, peakfile2,
        covfile1, covfile2,
        minprobdif, mergedist, minlen,
        outfld):

    start = time.time()
    os.makedirs(outfld, exist_ok = True)
    if not outfld.endswith('/'): outfld = outfld+'/'

    ## Generate skewnorm distribution model for coverage in peaks
    ## Peaks coverage
    print('Fitting distribution...')

    pf1 = pb.BedTool(peakfile1)
    pf2 = pb.BedTool(peakfile2)
    cf1 = pb.BedTool(covfile1)
    cf2 = pb.BedTool(covfile2)

    cov_peaks1 = np.array([float(f.fields[10]) for f in pf1.map(cf1, c = 4, o='mean')])
    cov_peaks2 = np.array([float(f.fields[10]) for f in pf2.map(cf2, c = 4, o='mean')])

    skewnorm_params1 = skewnorm.fit(cov_peaks1)
    skewnorm_params2 = skewnorm.fit(cov_peaks2)

    plot_skewnorm(cov_peaks1, skewnorm_params1, f'{outfld}peak_coverage_fit_1.png')
    plot_skewnorm(cov_peaks2, skewnorm_params2, f'{outfld}peak_coverage_fit_2.png')

    ## Create join coverage bdg
    print('Creating join coverage bed...')

    outfile = outfld+'common_coverage.bdg'
    cmd = ['bedtools', 'unionbedg', '-i', covfile1, covfile2, '>', outfile]
    sp.call(' '.join(cmd), shell = True)
    union_bed = pb.BedTool(outfile)

    ## Create common peaks bed
    print('Creating common peaks bed...')

    files_to_cross = []
    str_beds = peakfile1 + ' ' + peakfile2
    cmd = 'awk \'{print}\' '+f'{str_beds} > {outfld}common_peaks.bed'
    sp.call(cmd, shell = True)
    common_peaks = pb.BedTool(f'{outfld}common_peaks.bed').sort()

    ## Subset union_bed to common_peaks
    common_cov_common_peaks = union_bed.intersect(common_peaks)

    ## Call differential peaks
    print('Calling differential peaks (this might take a while)...')

    c1s = [float(x.fields[3]) for x in common_cov_common_peaks]
    c2s = [float(x.fields[4]) for x in common_cov_common_peaks]

    with Pool(8) as p:
        prob1s = p.starmap(get_probs, zip(c1s, repeat(skewnorm_params1)))
    with Pool(8) as p:
        prob2s = p.starmap(get_probs, zip(c2s, repeat(skewnorm_params2)))

    prob1s = np.array(prob1s)
    prob2s = np.array(prob2s)

    probdifs = prob1s - prob2s
    peaks1over2 = probdifs > minprobdif
    peaks2over1 = -probdifs > minprobdif
    chroms = [x.chrom for x in common_cov_common_peaks]
    starts = [x.start for x in common_cov_common_peaks]
    stops = [x.stop for x in common_cov_common_peaks]

    pd.DataFrame(common_cov_common_peaks)
    df = pd.read_table(
        common_cov_common_peaks.fn,
        names=['#chrom', 'start', 'stop', 'cov1', 'cov2']
               )
    rawout1 = f'{outfld}peaks_1over2_mpd{minprobdif}_merge{mergedist}_minlen{minlen}.bed'
    rawout2 = f'{outfld}peaks_2over1_mpd{minprobdif}_merge{mergedist}_minlen{minlen}.bed'

    df[peaks1over2].to_csv(rawout1, sep='\t', index=False)
    df[peaks2over1].to_csv(rawout2, sep='\t', index=False)

    rawbed1 = pb.BedTool(rawout1)
    rawbed2 = pb.BedTool(rawout2)

    out1 = rawbed1.sort().merge(d=mergedist).filter(lambda f: f.stop - f.start > minlen).saveas(f'{outfld}filtered_peaks_1over2.bed')
    out2 = rawbed2.sort().merge(d=mergedist).filter(lambda f: f.stop - f.start > minlen).saveas(f'{outfld}filtered_peaks_2over1.bed')

    end = time.time()
    print('Finished! Elapsed time:')
    print(end - start)
   

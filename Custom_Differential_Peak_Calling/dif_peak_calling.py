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
from itertools import combinations

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


wd = '/mnt/Disc4T/Projects/Miniprojects/Custom_Differential_Peak_Calling/'
os.chdir(wd)

## Call on 1 file
peakfile1 = 'A7K9_me_Macspeaks_peaks.narrowPeak'
peakfile2 = 'B11_me_Macspeaks_peaks.narrowPeak'
covfile1 = 'A7K9_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg'
covfile2 = 'B11_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg'
minprobdif = 0.4
mergedist = 500
minlen = 500
outfld = peakfile1.split('_')[0] + '_vs_' + peakfile2.split('_')[0]

get_differential_peaks(
    peakfile1, peakfile2,
    covfile1, covfile2,
    minprobdif, mergedist, minlen,
    outfld
)


## Call on all files
## Load MACS2 peaks file and normalized by input coverage

datadir = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/'

macs2_dir = datadir+'/Peak_Calling_MACS2/'
cov_dir = datadir+'/RPKMs_NormInput_noDup_bs10_smth_200_pseudo_10/'

macs2_fls = sorted([f for f in os.listdir(macs2_dir) if f.endswith('narrowPeak') and '_me_' in f])
cov_fls = sorted([f for f in os.listdir(cov_dir) if f.endswith('_pseudo10.bdg') and '_me_' in f])

macs2_fls = [
    '1.2B_me_Macspeaks_peaks.narrowPeak',
    '10G_me_Macspeaks_peaks.narrowPeak',
    'A7K9_me_Macspeaks_peaks.narrowPeak',
    'E5K9_me_Macspeaks_peaks.narrowPeak',
    'B11_me_Macspeaks_peaks.narrowPeak',
]

cov_fls = [
    '1.2B_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
    '10G_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
    'A7K9_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
    'E5K9_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
    'B11_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg',
]

out_generalfld = '/mnt/Disc4T/Projects/Miniprojects/Custom_Differential_Peak_Calling/DifPeaks/'
os.makedirs(out_generalfld)

for c in zip(combinations(macs2_fls, 2), combinations(cov_fls, 2)):
    print(c[0])
    print(c[1])

    peakfile1 = macs2_dir+c[0][0]
    peakfile2 = macs2_dir+c[0][1]
    covfile1 = cov_dir+c[1][0]
    covfile2 = cov_dir+c[1][1]
    minprobdif = 0.4
    mergedist = 500
    minlen = 500
    outfld = out_generalfld+c[0][0].split('_')[0] + '_vs_' + c[0][1].split('_')[0]
 
    get_differential_peaks(
        peakfile1, peakfile2,
        covfile1, covfile2,
        minprobdif, mergedist, minlen,
        outfld
    )


############## 'Old' stuff #####################

## Generate skewnorm distribution model for coverage in peaks
## Peaks coverage
cov_peaks1 = np.array([float(f.fields[10]) for f in macs2_1.map(cov_1, c = 4, o='mean')])
cov_peaks2 = np.array([float(f.fields[10]) for f in macs2_2.map(cov_2, c = 4, o='mean')])

skewnorm_params1 = skewnorm.fit(cov_peaks1)
skewnorm_params2 = skewnorm.fit(cov_peaks2)

## Create join coverage bdg
b1 = cov_dir+'10G_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg'
b2 = cov_dir+'E5K9_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg'
outfile = 'common_coverage.bdg'
cmd = ['bedtools', 'unionbedg', '-i', b1, b2, '>', outfile]
sp.call(' '.join(cmd), shell = True)
union_bed = pb.BedTool(outfile)

## Create common peaks bed
files_to_cross = [
    '10G_me_Macspeaks_peaks.narrowPeak',
    'E5K9_me_Macspeaks_peaks.narrowPeak'
]
str_beds = ' '.join([macs2_dir+f for f in files_to_cross])
cmd = 'awk \'{print}\' '+f'{str_beds} > common_peaks.bed'
sp.call(cmd, shell = True)
common_peaks = pb.BedTool('common_peaks.bed').sort()

## Subset union_bed to common_peaks
common_cov_common_peaks = union_bed.intersect(common_peaks)

## Call differential peaks
minprobdif = 0.2

## BedTools.filter approach
def compare_cov(feat):
    c1 = float(feat.fields[3])
    c2 = float(feat.fields[4])
    prob1 = skewnorm.cdf(c1, *skewnorm_params1)
    prob2 = skewnorm.cdf(c2, *skewnorm_params2)
    probdif = prob1 - prob2
    return(probdif)

peaks1over2 = common_cov_common_peaks.filter(lambda f: compare_cov(f) > minprobdif).saveas('peaks_1over2.bed')
peaks2over1 = common_cov_common_peaks.filter(lambda f: -compare_cov(f) > minprobdif).saveas('peaks_2over1.bed')

## For Loop approach
peaks1_2 = ''
peaks2_1 = ''
for cov_bin in tqdm(common_cov_common_peaks):
    c1, c2 = [float(x) for x in cov_bin.fields[3:5]]
    prob1 = skewnorm.cdf(c1, *skewnorm_params1)
    prob2 = skewnorm.cdf(c2, *skewnorm_params2)
    probdif = prob1 - prob2
    #print(probdif)
    if (abs(probdif) > minprobdif) and (prob1 > prob2):
        peaks1_2 += '\t'.join(cov_bin.fields[0:3]+[str(probdif)])+'\n'
    elif (abs(probdif) > minprobdif) and (prob2 > prob1):
        peaks2_1 += '\t'.join(cov_bin.fields[0:3]+[str(probdif)])+'\n'

with open('peaks_1_2.bed', 'w+') as outfile: outfile.write(peaks1_2)
with open('peaks_2_1.bed', 'w+') as outfile: outfile.write(peaks2_1)


print(common_cov_common_peaks[1])
for f in common_cov_common_peaks[0:5]:
    print(f)

## Vectorial approach


c1s = [float(x.fields[3]) for x in common_cov_common_peaks]
c2s = [float(x.fields[4]) for x in common_cov_common_peaks]

with Pool(8) as p:
    prob1s = p.starmap(get_probs, zip(c1s, repeat(skewnorm_params1)))
with Pool(8) as p:
    prob2s = p.starmap(get_probs, zip(c2s, repeat(skewnorm_params2)))

prob1s = np.array(prob1s)
prob2s = np.array(prob2s)
#prob1s = skewnorm.cdf(c1s, *skewnorm_params1)
#prob2s = skewnorm.cdf(c2s, *skewnorm_params2)
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
df[peaks1over2].to_csv('pandas_peak1over2.bed', sep='\t', index=False)
df[peaks2over1].to_csv('pandas_peak2over1.bed', sep='\t', index=False)


############### Plot Distribution ###################

a, loc, scale = skewnorm.fit(cov_peaks)

# Plot the histogram.
plt.hist(cov_peaks, bins=25, density=True, alpha=0.6, color='g')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = skewnorm.pdf(x, a, loc, scale)
plt.plot(x, p, 'k', linewidth=2)
#title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
plt.title(title)

plt.show()


skewnorm.cdf(0, a, loc, scale)


# ## Rest of the genome coverage
# whole_genome_bed = pb.BedTool('./whole_genome.bed')
# non_peak_bed = whole_genome_bed.subtract(macs2_10G)
# cov_nopeaks = np.array([float(f.fields[3]) for f in non_peak_bed.map(cov_10G, c = 4, o='mean')])

# a, loc, scale = skewnorm.fit(cov_nopeaks)

# # Plot the histogram.
# plt.hist(cov_nopeaks, bins=25, density=True, alpha=0.6, color='g')

# # Plot the PDF.
# xmin, xmax = plt.xlim()
# x = np.linspace(xmin, xmax, 100)
# p = skewnorm.pdf(x, a, loc, scale)
# plt.plot(x, p, 'k', linewidth=2)
# #title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
# plt.title(title)

# plt.show()


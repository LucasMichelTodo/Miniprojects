import numpy as np
import pybedtools as pb
from sklearn.mixture import GaussianMixture
from scipy.stats import skewnorm
import os
import matplotlib.pyplot as plt

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

#### Calls ####

wd = '/mnt/Disc4T/Projects/Miniprojects/Custom_Differential_Peak_Calling/'
os.chdir(wd)
outfld = './Distributions_Tests/'
os.makedirs(outfld, exist_ok=True)

peakfile1 = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/Peak_Calling_MACS2/1.2B_me_Macspeaks_peaks.narrowPeak'
peakfile2 = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/Peak_Calling_MACS2/10G_me_Macspeaks_peaks.narrowPeak'
covfile1 = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/RPKMs_NormInput_noDup_bs10_smth_200_pseudo_10/1.2B_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg'
covfile2 = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/RPKMs_NormInput_noDup_bs10_smth_200_pseudo_10/10G_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg'


pf1 = pb.BedTool(peakfile1)
pf2 = pb.BedTool(peakfile2)
cf1 = pb.BedTool(covfile1)
cf2 = pb.BedTool(covfile2)

cov_peaks1 = np.array([float(f.fields[10]) for f in pf1.map(cf1, c = 4, o='mean')])
cov_peaks2 = np.array([float(f.fields[10]) for f in pf2.map(cf2, c = 4, o='mean')])

skewnorm_params1 = skewnorm.fit(cov_peaks1)
skewnorm_params2 = skewnorm.fit(cov_peaks2)

suffix = '_peak_coverage_fit.png'
prefix1 = 'plot_12B'
prefix2 = 'plot_10G'

plot_skewnorm(cov_peaks1, skewnorm_params1, f'{outfld}{prefix1}{suffix}')
plot_skewnorm(cov_peaks2, skewnorm_params2, f'{outfld}{prefix2}{suffix}')

gm = GaussianMixture(n_components=2).fit(cov_peaks1.reshape(-1, 1))
cov_vector = cov_peaks1
plotname = f'{outfld}gaussian_mix.png'

# Plot the histogram.
plt.hist(cov_vector, bins=25, density=True, alpha=0.6, color='g')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
logprob = gm.predict_proba(x.reshape(-1, 1))
pdf = logprob[:,1]
plt.plot(x, pdf, 'k', linewidth=2)
plt.savefig(plotname)
plt.clf()

prob = gm.predict_proba(cov_peaks1.reshape(-1, 1))

prob
plt.hist(prob, bins=25, density=True, alpha=0.6, color='g')

cov_peaks1
x = np.array([[9.98565995e-01, 1.43400470e-03],
              [4.26389699e-03, 9.95736103e-01],
              [9.99963980e-01, 3.60198451e-05],
              [9.99695195e-01, 3.04805350e-04],
              [9.96041044e-01, 3.95895619e-03],
              [9.99753861e-01, 2.46139290e-04]])


x[:,0]

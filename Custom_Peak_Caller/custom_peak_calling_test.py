import pybedtools as pb
import numpy as np
from sklearn.mixture import GaussianMixture
import warnings
warnings.filterwarnings( "ignore", module = "matplotlib\..*" )
import matplotlib.pyplot as plt
import os

#### Functions ####

def plot_gaussian_mixture(cov_vector, gm, plotname):

    # Plot the histogram.
    plt.hist(cov_vector, bins=25, density=True, alpha=0.6, color='g')

    # Plot the PDF.
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    logprob = gm.score_samples(x.reshape(-1, 1))
    pdf = np.exp(logprob)
    plt.plot(x, pdf, 'k', linewidth=2)
    plt.savefig(plotname)
    plt.clf()

def plot_components(cov_vector, gm, plotname):
    # Plot the histogram.
    plt.hist(cov_vector, bins=25, density=True, alpha=0.6, color='g')

    # Plot the PDF.
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    logprob = gm.predict_proba(x.reshape(-1, 1))
    cdf1 = logprob[:,0]
    cdf2 = logprob[:,1]
    plt.plot(x, cdf1, 'k', linewidth=2, color='blue', label="Component 1")
    plt.plot(x, cdf2, 'k', linewidth=2, color='green', label="Component 2")
    plt.legend(loc='upper left')
    plt.savefig(plotname)
    plt.clf()

def get_probs(coverage, gm):
    return(gm.predict_proba(coverage.reshape(-1, 1)))

#### Calls ####

wd = '/mnt/Disc4T/Projects/Miniprojects/Custom_Peak_Caller/'
os.chdir(wd)

datadir = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/RPKMs_NormInput_noDup_bs10_smth_200_pseudo_10/'
covfile = '1.2B_me_sort_q5_noDup_rpkm_normInput_bs10_smth200_pseudo10.bdg'
cov_field = 3
out_name = 'peaks_joined_filtered.bed'
mergedist = 50
minlen = 200

covbed = pb.BedTool(datadir+covfile)
cov = np.array([float(feat.fields[cov_field]) for feat in covbed])
gm = GaussianMixture(n_components=2).fit(cov.reshape(-1, 1))

plot_gaussian_mixture(cov, gm, f'gm_fit.png')
plot_components(cov, gm, f'gm_components.png')


## All sklearn estimators use as input a 2D array,
## with samples as rows and features as columns.
## Our data-set consists in one sample per bed feature
## with only one feature (coverage), so we have to reshape it
## into a 2D array with just one column.
## -1 in reshape just means whatever it takes to make it work!
## So (-1, 1) means: any number of rows and 1 column

cov.reshape(-1,1).shape
probs = gm.predict_proba(cov.reshape(-1,1))

labels = gm.predict(cov.reshape(-1,1))
labels

np.mean(cov[labels == 0])
np.mean(cov[labels == 1])

np.unique(labels)
mask
mask = probs[:,0] < probs[:,1]
# labels = gm.predict(cov.reshape(-1,1)) produces same result as 0s and 1s

with open(f'pre_{out_name}', 'w+') as outfile:
    for feat, mask in zip(covbed, mask):
        if mask:
            outfile.write(str(feat))

outbed_raw = pb.BedTool('peaks.bed')
out_pre = outbed_raw.sort().merge(d=mergedist)
out = out_pre.filter(lambda f: f.stop - f.start > minlen).saveas(out_name)


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


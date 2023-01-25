import os
import pandas as pd
from chip_seq_processing import *
from bowtie2_align import *

## wd = '/mnt/Disc4T/Projects/Miniprojects/Chip_Seq_Pipe/'
## os.chdir(wd)

samples = pd.read_excel('./chip_seq_samples.xlsx', engine='openpyxl', skiprows=1)
bamdir = './Alignments/DeDuplicated/'

def get_align_name(r1, r2):
    name1 = r1.rsplit('.', 1)[0]
    name2 = r2.rsplit('.', 1)[0]
    name = common_start(name1, name2)+'_q5_sort_noDup.bam'
    return(name)

bams = sorted([b for b in os.listdir(bamdir) if b.endswith('_noDup.bam')])
IPs = [get_align_name(r1, r2) for r1, r2 in zip(samples['Read 1'], samples['Read 2'])]
inputs = [get_align_name(r1, r2) for r1, r2 in zip(samples['Read 1.1'], samples['Read 2.1'])]


#### Get RPKMs ####

outdir = './RPKMs_noDup_bs10_smth_200/'
os.makedirs(outdir, exist_ok = True)

bs = 10
smooth = 200
norm = 'RPKM'

for bam in bams:
    out = outdir+bam.replace('.bam', '_RPKMs.bdg')
    get_RPKMs(bamdir+bam, bs, smooth, norm, outdir)

#### Get normalized by input RPKMs ####

outdir = './RPKMs_NormInput_noDup_bs10_smth_200_pseudo_10/'
os.makedirs(outdir, exist_ok = True)

bs = 10
smooth = 200
norm = 'RPKM'
of = 'bedgraph'
outfld = outdir
pseudo = 10
num_process = 8

for ip, inpt in zip(IPs, inputs):
    print(ip, inpt)
    print('----')

    get_RPKMs_normInput(
        bamdir+ip,
        bamdir+inpt,
        bs,
        smooth,
        norm,
        of,
        outfld,
        pseudo
    )

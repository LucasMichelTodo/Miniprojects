import os
import subprocess as sp
import multiprocessing as mp

## Functions

def subset_fastq(fastq, nreads):
    ext = '.'+'.'.join(fastq.split('.')[1:])
    out_fq = fastq.replace(ext, f'_subset{nreads}.fq')

    cmd = ['seqtk', 'sample', '-s123', fastq, str(nreads)]
    print(f'Subseting {fastq} to {out_fq} ({str(nreads)} reads) ...')
    subset = sp.run(cmd, stdout=sp.PIPE).stdout.decode('utf-8')
    with open(out_fq, 'w+') as outfile:
        outfile.write(subset)

def subset_fastq_paired_parallel(fastq1, fastq2, nreads):
    pool = mp.Pool(2)
    pool.starmap(subset_fastq, [[fastq1, nreads], [fastq2, nreads]])

## Calls

wd = '/mnt/Disc4T/Projects/Nuria/ChIP_151021_NCV/Raw_Data/'
os.chdir(wd)

fastq1s = [
    'NCV26_lib_08246AAC_ATCACG_R1_001.fastq.gz',
    'NCV27_lib_08247AAC_CGATGT_R1_001.fastq.gz',
    'NCV28_lib_08248AAC_TTAGGC_R1_001.fastq.gz',
    'NCV29_lib_08249AAC_TGACCA_R1_001.fastq.gz',
]
fastq2s = [
    'NCV26_lib_08246AAC_ATCACG_R2_001.fastq.gz',
    'NCV27_lib_08247AAC_CGATGT_R2_001.fastq.gz',
    'NCV28_lib_08248AAC_TTAGGC_R2_001.fastq.gz',
    'NCV29_lib_08249AAC_TGACCA_R2_001.fastq.gz',
]

for f1, f2 in zip(fastq1s, fastq2s):
    subset_fastq_paired_parallel(f1, f2, 1000)

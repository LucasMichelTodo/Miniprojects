import os
import sys
import pandas as pd
from preprocessing import *
from bowtie2_align import *

## Input Data

## wd = '/mnt/Disc4T/Projects/Miniprojects/Chip_Seq_Pipe/'
## os.chdir(wd)

raw_data_dir = './Raw_Data/'
samples = pd.read_excel('./chip_seq_samples.xlsx', engine='openpyxl', skiprows=1)

read1s = samples['Read 1'].tolist() + samples['Read 1.1'].tolist()
read2s = samples['Read 2'].tolist() + samples['Read 2.1'].tolist()

print('Cleaning Reads... ')
for f in read1s: print(f)
print('-----')
for f in read2s: print(f)
print('-----')

## Clean Reads

params = "ktrim=r k=22 mink=6 overwrite=t "
outpath = './Clean_Reads/'
os.makedirs(outpath, exist_ok=True)

for pair in zip(read1s, read2s):
    in1, in2 = raw_data_dir+pair[0], raw_data_dir+pair[1]
    call_BBDUK(in1, in2, outpath, params)

## Call fastq_screen

all_reads = [raw_data_dir + f for f in read1s+read2s]
outdir = './Fastq_Screens/'
os.makedirs(outdir, exist_ok=True)
call_fastq_screen(all_reads, threads = '8', outdir=outdir)

## Call fastqc on clean reads
print('Calling FastQC on clean reads...\n\n')

clean_reads_dir = './Clean_Reads/'
all_reads = [clean_reads_dir+f for f in os.listdir(clean_reads_dir) if f.endswith('clean.fq')]
outdir = './FastQC_Clean/'
os.makedirs(outdir, exist_ok=True)
call_fastqc(all_reads, threads = '8', outdir = outdir)

## Call fastq_screen on clean reads
print('Calling Fastq_Screen on clean reads...\n\n')

outdir = './Fastq_Screens_Clean/'
os.makedirs(outdir, exist_ok=True)
call_fastq_screen(all_reads, threads = '8', outdir=outdir)


## Algin Clean Reads
print('Starting Alignments...\n\n')

indir = './Clean_Reads/'
outdir = './Alignments/'
os.makedirs(outdir, exist_ok=True)

all_reads = sorted([f for f in os.listdir(indir) if f.endswith('clean.fq')])
reads1 = [indir+r for r in all_reads if '_R1_' in r or '_read1' in r]
reads2 = [indir+r for r in all_reads if '_R2_' in r or '_read2' in r]

params = ("-p 8 --very-sensitive --local "
          "-5 4 -3 4 -I 50 -X 2000 "
          "-x /home/lucas/Programs/bowtie2-2.3.0-legacy/Pf3D7")

bowtie2_align_pipe(reads1, reads2, params, outdir)

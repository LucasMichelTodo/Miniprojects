from tqdm import tqdm
import sys
import subprocess as sp
import os

""" This is a preprocessing script for raw-read files

It takes as fisrt input raw read fastq files and sequentially:
- Cleans reads using BBDUK
- Aligns reads to a reference genome using Bowtie2
- Converts SAM files into BAM files and filters them to q >=5.

"""

############ Functions ############


def call_BBDUK(in1, in2, out1, out2, outm, ref, params):
    cmd = ("bbduk.sh in={} in2={} "
           "out={} out2={} outm={} "
           "ref={}") .format(in1, in2, out1, out2, outm, ref)

    cmd = cmd+" "+params
    sp.call(cmd, shell=True)


def call_Bowtie2(in1, in2, out, params):
    cmd = "bowtie2 -1 {} -2 {} -S {}" .format(in1, in2, out)
    cmd = cmd+" "+params
    print(cmd)
    sp.call(cmd, shell=True)


def from_sam_to_bam(samfile):
    name = samfile.rsplit(".")[0]
    cmd = "samtools view -bS {} > {}" .format(samfile, name+".bam")
    sp.call(cmd, shell=True)

    # Erase SAM after creating BAM
    cmd = "rm {}" .format(samfile)
    sp.call(cmd, shell=True)

    cmd = "samtools sort {} > {}" .format(name+".bam", name+"_sort.bam")
    sp.call(cmd, shell=True)

    # Erase bam after creating sortedBAM
    cmd = "rm {}" .format(name+".bam")
    sp.call(cmd, shell=True)

    cmd = "samtools index {} > {}" .format(
        name+"_sort.bam", name+"_sort.bam.bai")
    sp.call(cmd, shell=True)

    # Filter only >=q5 reads
    cmd = "samtools view -b -q 5 {} > {}" .format(
        name+"_sort.bam", name+"_q5_sort.bam")
    sp.call(cmd, shell=True)


############ Calls ############


# BBDUK
params = "ktrim=r k=22 mink=6 overwrite=t"
ref = "/home/lucas/Programs/bbmap/resources/adapters.fa"
inpath = ""

read1s = sorted([f for f in os.listdir(inpath) if "read1.fastq.gz" in f])
read2s = sorted([f for f in os.listdir(inpath) if "read2.fastq.gz" in f])

for pair in tqdm(zip(read1s, read2s)):
    name = pair[0].split("_")[0]

    in1, in2 = inpath+pair[0], inpath+pair[1]
    out1, out2 = inpath+name+"_read1_clean.fastq", inpath+name+"_read2_clean.fastq"
    outm = inpath+name+"_badreads.fastq"

    call_BBDUK(in1, in2, out1, out2, outm, ref, params)

# BOWTIE2
params = ("-p 4 --very-sensitive --local "
          "-5 4 -3 4 -I 50 -X 200 "
          "-x /home/lucas/Programs/bowtie2-2.3.0-legacy/Pf3D7")

inpath = ""
files = [f for f in os.listdir(inpath)if f.endswith("_clean.fastq")]

read1s = sorted([f for f in files if "read1" in f])
read2s = sorted([f for f in files if "read2" in f])

for pair in tqdm(zip(read1s, read2s)):
    name = pair[0].split("_")[0]

    in1, in2 = inpath+pair[0], inpath+pair[1]
    out = inpath+name+".sam"
    call_Bowtie2(in1, in2, out, params)


# From SAM to BAM
indir = ""
samfiles = [f for f in os.listdir(indir) if f.endswith(".sam")]

for f in tqdm(samfiles):
    sam = indir+f
    from_sam_to_bam(sam)

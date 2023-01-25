import os
import subprocess as sp

def remove_duplicates(indir, outdir, bam):

    gatk_path = '/home/lucas/Programs/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar'
    i = indir+bam
    o = outdir+bam.replace(".bam", "_noDup.bam")
    m = outdir+bam.replace(".bam", "_metrics.txt")

    cmd = (f"java -jar {gatk_path} MarkDuplicates "
           "REMOVE_DUPLICATES=true I={} O={} M={}") .format(i, o, m)

    sp.call(cmd, shell=True)

def samtools_index(bam, np=8):
    cmd = ['samtools', 'index', '-@', str(np), bam]
    print(' '.join(cmd))
    sp.run(cmd)

## Remove Duplicates
indir = './Alignments/'
bams = sorted([b for b in os.listdir(indir) if b.endswith("q5_sort.bam")])
outdir = indir+'DeDuplicated/'
os.makedirs(outdir, exist_ok=True)

for bam in bams:
    remove_duplicates(indir, outdir, bam)

## Index files
nodup_bams = sorted([f for f in os.listdir(outdir) if f.endswith('_noDup.bam')])

for b in nodup_bams:
    samtools_index(outdir+b)

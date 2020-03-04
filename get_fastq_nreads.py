import os
import subprocess as sp

# Count number of reads in all fastq. files specified in [indirs]


indirs = ['/mnt/Disc4T/Projects/Choline_DHA/Broadbent/Raw_data/']

for indir in indirs:

    ext = (".fastq", ".fastq.gz") 
    fastqs = [f for f in os.listdir(indir) if f.endswith(ext)]

    for f in fastqs:
        print(f'Counting {f}:')
        sp.call(f'expr $(cat {indir+f} | wc -l) / 4', shell=True)
        print("----")
    


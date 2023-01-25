import os
import subprocess as sp

wd = '/mnt/Disc4T/Projects/Miniprojects/Heterochromatin_Boundaries/Strong_Insulators/'
os.chdir(wd)

intervals = './Intervals/manual_intervals_subclones.bed'
ref_fasta = '../Data/PlasmoDB-61_Pfalciparum3D7_Genome.fasta'
outfile = './Seqs/manually_curated_subclones.fasta'

cmd = f'bedtools getfasta -fi {ref_fasta} -bed {intervals} > {outfile}'
sp.call(cmd, shell = True)

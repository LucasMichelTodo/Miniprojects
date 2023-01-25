import os
import subprocess as sp

wd = '/mnt/Disc4T/Projects/Miniprojects/Heterochromatin_Boundaries/Strong_Insulators/'
os.chdir(wd)

outdir = 'Motifs/'
os.makedirs(outdir, exist_ok=True)

infile = './Seqs/merged.fasta'
control = '../Data/PlasmoDB-61_Pfalciparum3D7_Genome.fasta'
seed = 123
nmotifs = 5
function = 'de'
numprocess = 8

meme_call = (
    f'meme {infile} '
    f'-o {outdir} '
    f'-objfun {function} '
    f'-neg {control} '
    '-dna -revcomp '
    f'-nmotifs {nmotifs} '
    f'-seed {seed} '
    f'-p {numprocess} '
)

meme_call

sp.call(meme_call, shell = True)

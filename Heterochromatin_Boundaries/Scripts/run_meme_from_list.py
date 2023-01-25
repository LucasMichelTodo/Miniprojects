import os
import subprocess as sp

indir = '/PROJECTES/MALARIA_EPIGENETICS/Projects/Het_Boundaries/Fastas/'
infiles = [
	'ASR_border_flanked.fasta',
	'GII_border_flanked.fasta',
	'manually_curated_subclones.fasta'
]

for f in infiles:

	infile = indir+f
	control = indir+'PlasmoDB-61_Pfalciparum3D7_Genome.fasta'
	seed = 123
	nmotifs = 10
	function = 'de'
	numprocess = 32

	suffix = f.replace('.fasta', '')
	outdir = f'/PROJECTES/MALARIA_EPIGENETICS/Projects/Het_Boundaries/Motifs_{suffix}/'
	os.makedirs(outdir, exist_ok=True)

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

	print(meme_call)
	#sp.call(meme_call, shell = True)
        


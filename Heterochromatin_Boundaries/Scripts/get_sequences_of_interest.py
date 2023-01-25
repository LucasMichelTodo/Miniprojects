import os
import pybedtools as pb
import numpy as np

wd = '/mnt/Disc4T/Projects/Miniprojects/Heterochromatin_Boundaries/'
os.chdir(wd)

#### Functions ####

def get_intervals_from_boundaries(bed, len_left, len_right, outdir):
    outname = os.path.basename(bed).replace('.bed', '_intervals.bed')
    with open(bed, 'r+') as infile, open(outdir+outname, 'w+') as outfile:
        linelist = [l.strip().split('\t') for l in infile.readlines()]
        chroms = sorted(np.unique([x[0] for x in linelist]))

        for chrom in chroms:
            chrom_list = [x for x in linelist if x[0] == chrom]
            for idx, line in enumerate(chrom_list):
                chrom, start, stop = line
                start = int(start)
                stop = int(stop)
                if idx == 0:
                    newstart = str(stop-len_left)
                    newstop = str(stop+len_right)
                    outfile.write('\t'.join([chrom, newstart, newstop])+'\n')
                elif idx == len(chrom_list)-1:
                    newstart = str(start-len_left)
                    newstop = str(start+len_right)
                    outfile.write('\t'.join([chrom, newstart, newstop])+'\n')
                else:
                    newstart1 = str(start-len_left)
                    newstop1 = str(start+len_right)
                    newstart2 = str(stop-len_left)
                    newstop2 = str(stop+len_right)
                    outfile.write('\t'.join([chrom, newstart1, newstop1])+'\n')
                    outfile.write('\t'.join([chrom, newstart2, newstop2])+'\n')

#### Calls ####

indir = './Sandra/'
bed_limits_files = [
    'FRASCHKA_fulllength_merged_60000bp_gamIV.bed',
    'BUNNIK_fulllengthtelomere_GamIV_rep2.bed',
    'BUNNIK_fulllengthtelomere_GamIV_rep1.bed',
    'US_fulltel_merged_60000bp_G2rep1.bed',
]

outdir = './Strong_Insulators/Seqs/'
os.makedirs(outdir, exist_ok=True)

for bed in bed_limits_files:
    print(indir+bed)
    get_intervals_from_boundaries(indir+bed, 5000, 5000, outdir)

import os
import pybedtools as pb
import numpy as np
from Bio import SeqIO

wd = '/mnt/Disc4T/Projects/Miniprojects/Heterochromatin_Boundaries/'
os.chdir(wd)

#### Functions ####

def fasta_seqs_from_bed(bed, ref_fasta, outname):

    parser = SeqIO.parse(open(ref_fasta), 'fasta')
    dict_fasta = dict([(seq.id, seq) for seq in parser])
    inbed = pb.BedTool(bed)

    with open(outname, 'w+') as outfile:
        for feat in inbed:
            if feat.chrom in dict_fasta:
                file_id = outname.split('/')[-1].replace('.fasta', '')
                outfile.write(f'>{feat.chrom}_{feat.start}:{feat.stop}_{file_id}\n')
                outfile.write(str(dict_fasta[feat.chrom][feat.start:feat.stop+1].seq)+'\n')

#### Calls ####

ref_fasta = './Data/PlasmoDB-61_Pfalciparum3D7_Genome.fasta'
indir = './Strong_Insulators/Intervals/'
bed_files = os.listdir(indir)
outdir = './Strong_Insulators/Seqs/'
os.makedirs(outdir, exist_ok=True)

for bed in bed_files:
    outname = bed.replace('.bed', '_SEQs.fasta')
    fasta_seqs_from_bed(indir+bed, ref_fasta, outdir+outname)

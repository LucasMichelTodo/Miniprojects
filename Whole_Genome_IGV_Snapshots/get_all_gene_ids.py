import os
import pybedtools as pb

wd = '/mnt/Disc4T/Projects/Miniprojects/PhD_Project_Gene_Plots/'
os.chdir(wd)
ref_gff = pb.BedTool('PlasmoDB-52_Pfalciparum3D7.gff')

entry_types = set([entry.fields[2] for entry in ref_gff])
gene_types = [
    'ncRNA_gene',
    'protein_coding_gene',
    'pseudogene',
]

with open('all_gene_ids.txt', 'w+') as output:
    for entry in ref_gff:
        if entry.fields[2] in gene_types:
            info = entry.fields[8].split(';')
            info_dict = {k:v for k,v in [x.split('=') for x in info]}
            gid = info_dict['ID']
            output.write(gid+'\n')
            

import os
import re
from collections import defaultdict

wd = '/mnt/Disc4T/Projects/Miniprojects/Gene_Rosetta/'
os.chdir(wd)

rosetta = defaultdict(list)
with open('./plasmoDB_57_gene_aliases.txt', 'r+') as infile:
    for line in infile:
        linelist = line.strip().split('\t')
        new_id = linelist[0]
        old_ids = linelist[1:]
        for i in old_ids:
            rosetta[i].append(new_id)

with open('gene_rosetta.tsv', 'w+') as outfile:
    outfile.write('Old_id\tGene_id\n')
    for k, v in rosetta.items():
        outfile.write(k+'\t'+','.join(v)+'\n')
 
one_old_multi_new = {}
for k,v in rosetta.items():
    if len(v) > 1:
        one_old_multi_new[k] = v

with open('one_old_id_multiple_new.txt', 'w+') as outfile:
    for k, v in one_old_multi_new.items():
        outfile.write(k+'\n')


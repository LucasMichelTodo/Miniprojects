#!/usr/bin/env python

# Import packages
import pybedtools as py
from pybedtools.featurefuncs import TSS
from pybedtools.featurefuncs import three_prime
import subprocess as sp

# Parameters to set
upstream = 1000
dwstream = 1000
ref = '/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-41_Pfalciparum3D7.gff'

out_name = "anotation_PlasmoDB-41"

# Create genome dict
genome_file = "/home/lucas/ISGlobal/Gen_Referencies/Pf3D7.genome"

genome={}
with open(genome_file, "r+") as infile:
    for line in infile:
        genome[line.split()[0]] = (0, int(line.split()[1]))

gff = py.BedTool(ref)
sel = ["gene"]#, "rRNA", "snoRNA", "snRNA", "tRNA", "ncRNA"]

gene_gff = gff.filter(lambda x: x.fields[2] in sel)

fivePrime = gene_gff.each(TSS,
                          upstream=upstream,
                          downstream=1,
                          add_to_name="_fivePrime",
                          genome=genome).saveas(out_name+"_fivePrime.gff")

gene_gff = gff.filter(lambda x: x.fields[2] in sel)

threePrime = gene_gff.each(three_prime,
                          upstream=1,
                          downstream=dwstream,
                           add_to_name="_threePrime",
                           genome=genome).saveas(out_name+"_threePrime.gff")

cmd = "cat {} {}_fivePrime.gff {}_threePrime.gff > {}.gff" .format(ref,
                                                                   out_name,
                                                                   out_name,
                                                                   out_name)
sp.call(cmd, shell = True)

cmd = "gff_sorter.py {}.gff > {}_sorted.gff" .format(out_name, out_name)

sp.call(cmd, shell = True)

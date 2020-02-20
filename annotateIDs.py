import os
import pybedtools as pb

gffile = "/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-46_Pfalciparum3D7.gff"
gff = pb.BedTool(gffile)
gff = gff.filter(lambda x: x[2] == "gene")

genesfile = "/mnt/Disc4T/Projects/PhD_Project/chip_seq_genes.txt"
genes = []
with open(genesfile, "r+") as infile:
    for line in infile:
        genes.append(line.strip())


def getAnnot(gffentry):

    info = gffentry.fields[8].split(";")
    gid = info[0].replace("ID=", "")
    anot = info[1].replace("description=", "")
    return([gid, anot])


annot = {}
for entry in gff:
    ga = getAnnot(entry)
    annot[ga[0]] = ga[1]

with open(genesfile.replace(".txt", "_annotated.csv"), "w+") as outfile:
    for gene in genes:
        outfile.write(gene+"\t"+annot[gene]+"\n")

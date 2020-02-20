import pybedtools as py

ref = '/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-46_Pfalciparum3D7.gff'

gff = py.BedTool(ref)
gene_gff = gff.filter(lambda x: x.fields[2] == "gene")

for gene in gene_gff:
    infoline = gene.fields[8]
    gid = infoline.split(";")[0].replace("ID=", "")
    print(gid)

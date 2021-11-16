import collections as cll
import pandas as pd

gff_file = ""
gaf_file = "/mnt/Disc4T/Projects/Anastasia/Arrays_Lipids/plasmoDB_52.gaf"


# def get_GO_from_gff(gff):
#     with open(gff, "r+") as gfile:
#         anot = []
#         for line in gfile:
#             if line.startswith("#"):
#                 pass
#             else:
#                 if line.split()[2] == "mRNA":
#                     info = line.split("\t")[8].split(";")
#                     go = 0
#                     for i in info:
#                         if i.startswith("Ontology_term="):
#                             go = i.replace("%3B", ",").replace(
#                                 "Ontology_term=", "").strip().split(",")
#                     if go == 0:
#                         go = ["Unknown function"]

#                     anot.append((info[0].replace("ID=", "").split(".")[0], go))

#     GOs = cll.defaultdict(list)
#     for gene in anot:
#         for g in gene[1]:
#             GOs[g].append(gene[0])

#     return(GOs)


def get_GO_from_gaf(gaf):
    with open(gaf, "r+") as gfile:
        GOs = cll.defaultdict(list)
        for line in gfile:
            if line.startswith(("#", "!")):
                pass
            else:
                go = line.split("\t")[4]
                # print(go)
                gene = line.split("\t")[1]
                # print(gene)
                GOs[go].append(gene)

    return(GOs)


#anot = get_GO_from_gff("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-41_Pfalciparum3D7.gff")

anot = get_GO_from_gaf(gaf_file)
tbl = pd.concat([pd.Series(v, name=k) for k, v in anot.items()], axis=1)

ncols = len(tbl.columns)
tbl.loc[-1] = ["na"]*ncols # Add a row of 'nas'
tbl.index = tbl.index + 1  # shifting index
tbl = tbl.sort_index()  # sorting by index
tbl.head()
tbl.to_csv(r'/mnt/Disc4T/Projects/Anastasia/Arrays_Lipids/gsea_from_PlsmoDB_52_gaf.csv',
           sep='\t', index=False)

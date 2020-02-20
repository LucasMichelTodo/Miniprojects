import os

os.chdir("/home/lucas/ISGlobal/Gen_Referencies/")

sizes = {}
with open("Pf3D7.sizes", "r+") as infile:
    for line in infile:
        chrom = line.strip().split()[0]
        n = line.strip().split()[1]
        sizes[chrom] = n

for chrom, n in sizes.items():
    i = 0
    while i < int(n):
        line = [chrom, str(i), str(i+1)]
        print("\t".join(line))
        i += 1

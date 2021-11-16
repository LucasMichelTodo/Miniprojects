import pybedtools as pb

ref = '/mnt/Disc4T/Projects/PhD_Project/Data/PlasmoDB-52_Pfalciparum3D7.gff'
outfile = ref.replace('.gff', '_parsed_annotation.tsv')

gff = pb.BedTool(ref)

types = set([entry.fields[2] for entry in gff])
gene_types = [
    'ncRNA_gene',
    'protein_coding_gene',
    'pseudogene'
]

gene_gff = gff.filter(lambda x: x.fields[2] in gene_types)

with open(outfile, 'w+') as out_file:

    outcols = [
        'Gene_id',
        'Chrom',
        'Start',
        'Stop',
        'Type',
        'Strand',
        'Name',
        'Annot'
    ]
    
    out_file.write('\t'.join(outcols)+'\n')
    
    for gene in gene_gff:
        line = [
            gene.attrs['ID'],
            gene.chrom,
            str(gene.start),
            str(gene.stop),
            gene.fields[2],
            gene.strand,
            gene.attrs.get('Name', ''),
            gene.attrs.get('description', '').replace('+', ' ').replace('%2C', ', ')
        ]
        
        out_file.write('\t'.join(line)+'\n')




    

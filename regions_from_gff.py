import pybedtools as pb
import pathlib as pth
from Bio import SeqIO

def get_extension(filename):
    '''
    Get extension from a filename. It can have up to 2 extensions,
    a filetype extension and optionally a compressing extension (e.g., '.gz')
    '''
    p = pth.Path(filename)
    s = p.suffixes[-2:]
    return(''.join(s))

def genome_dict(ref_fa):
    # Create genome dict
    genome={}
    records = list(SeqIO.parse(ref_fa, "fasta"))
    genome_file = ref_fa
    for record in records:
        genome[record.id] = (0, len(record.seq))
    return(genome)

def cut_TSS_and_CDS(gff, ref_fa,
                    up5=1000, cds='whole', up3=500,
                    allowoverlapp=False):

    genome = genome_dict(ref_fa)

    ref = pb.BedTool(gff)
    current_chrom = ''
    ngenes = len(ref)
    str_bed = ''
    first_in_chrom = False
    last_in_chrom = False

    extension = get_extension(gff)
    suffix = f'_{str(tss)}5p_{str(cds)}CDS.bed'
    outname = gff.replace(extension, suffix)
    with open(outname, 'w+') as outfile:

        for idx, gene in enumerate(ref):

            ## Check Orientation:
            strain = gene.fields[6]

            ## Check if first/last in chromosome
            chrom = gene.chrom

            if current_chrom != chrom:
                first_in_chrom = True
                # print('First in chrom!')
                # print(gene)
                # print('---------------')

            if idx == ngenes-1:
                ## First check if we are in the last gene!
                last_in_chrom = True
            else:
                if ref[idx+1].chrom != chrom:
                    last_in_chrom = True

            ## Set new start and new stop depending on strain:
            if strain == '+':
                newstart = gene.start-tss
                newstop = gene.start+cds
            else:
                newstart = gene.stop-cds
                newstop = gene.stop+tss

            # ## Check overlapp previous gene if +strain or next gene if -strain
            # if strain == '+':
            #     if first_in_chrom:
            #         pass
            #     else:
            #         if newstart < ref[idx-1].stop:
            #             newstart = ref[idx-1].stop+1
            #             # print('Previous gene:')
            #             # print(ref[idx-1])
            #             # print('Current gene:')
            #             # print(gene)
            #             # print('New start-stop:')
            #             # print(newstart, newstop)
            #             # print('--------------------------')
            # else:
            #     if last_in_chrom:
            #         pass
            #     else:
            #         if newstop > ref[idx+1].start:
            #             newstop = ref[idx+1].start-1
            #             # print('Next gene:')
            #             # print(ref[idx+1])
            #             # print('Current gene:')
            #             # print(gene)
            #             # print('New start-stop:')
            #             # print(newstart, newstop)
            #             # print('--------------------------')

            ## Check we dont go < 0
            if newstart < 0: newstart = 0

            ## Check we don't go > chrom length
            if newstop > genome[chrom][1]: newstop = genome[chrom][1]

            ## Check start always < stop
            if newstart >= newstop: newstop = newstart+1

            first_in_chrom = False
            last_in_chrom = False
            current_chrom = chrom

            newline = [gene.chrom, newstart, newstop, gene.fields[8].replace('\t', ' ')]
            newline = [str(x) for x in newline]
            outfile.write('\t'.join(newline)+'\n')
    print(f'Generated: {outname}')



ref_fa = '/home/lucas/ISGlobal/Projects/Alba/ChIP_Seqs_01_23_edited_genomes/Custom_Genomes/Genome_DD/genome_DD.fa'

gff = '/home/lucas/ISGlobal/Projects/Alba/ChIP_Seqs_01_23_edited_genomes/Custom_Genomes/Genome_DD/PlasmoDB-61_Pfalciparum3D7_edited_DD_final.gff'

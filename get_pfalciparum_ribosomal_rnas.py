from Bio import SeqIO

ref_fasta = '/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-54_Pfalciparum3D7_AnnotatedTranscripts.fasta'

out_fasta = '/home/lucas/ISGlobal/Gen_Referencies/pfalciparum_ribosomal_RNAs.fasta'

with open(out_fasta, 'w+') as outfile:
    for seq in SeqIO.parse(ref_fasta, "fasta"):
        info_list = seq.description.split(' | ')[1:]
        info = {x.split('=')[0]:x.split('=')[1] for x in info_list}
        if 'S ribosomal RNA' in info['gene_product']:
            outfile.write('>'+seq.description+'\n')
            outfile.write(str(seq.seq)+'\n')



            

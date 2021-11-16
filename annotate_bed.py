import pybedtools as pb
import sys

def extract_gff_info(feat, bdlen):
    """
    Extract information from the 9th field of the file
    generated intersecting a bed and a gff file.
    feat -> is the line of the file.
    bdlen -> number of fields of the original bedfile
    """
    if feat.fields[bdlen] == ".":
        gene, anot, strand = ".", ".", "."
    else:
        info = feat.fields[bdlen+8]
        gene = info.split(";")[0].replace("ID=", "")
        anot = info.split(";")[1].replace("description=", "")
        strand = feat.fields[bdlen+6]
    return(gene, anot, strand)


def annotate_bed(bed, gff, only_genes=True):
    """
    Annotate a valid bedfile.
    """

    # Make sure we are not appending to an existing file:
    outfile = bed.rsplit(".", 1)[0]+"_annotated.csv"
    out = open(outfile, "w+")

    anot, gene, dist, strand = ".", ".", ".", "."

    try:
        ref = pb.BedTool(gff)

        if only_genes:
            ref = ref.filter(lambda x: x[2] == "gene")

        ref = ref.sort()

        bedfile = pb.BedTool(bed)
        bdlen = len(bedfile[0].fields)

        # Write column titles
        bedcols = ['Bed_field_'+str(x) for x in range(1,bdlen-2)]
        out.write("\t".join(['Chrom', 'Start', 'Stop'] + bedcols +
                            ['Position', 'Gene_id', 'Annot', 'Dist', 'Strand'])+"\n")

        intersect = bedfile.intersect(ref, wao=True)

        for feat in intersect:
            # If there is a "perfect" hit report it and blank lines for nearest hits.
            if feat.fields[bdlen] != ".":
                leftline = feat.fields[0:bdlen]+['Nearest left', gene, anot, dist, strand]
                rightline = feat.fields[0:bdlen]+['Nearest right', gene, anot, dist, strand]

                gene, anot, strand = extract_gff_info(feat, bdlen)
                line = feat.fields[0:bdlen]+['Overlapping', gene, anot, dist, strand] 

                outfile = bed.rsplit(".", 1)[0]+"_annotated.csv"

                out.write("\t".join(line)+"\n")
                out.write("\t".join(leftline)+"\n")
                out.write("\t".join(rightline)+"\n")
                anot, gene, dist, strand = ".", ".", ".", "."
                
            else:
                # If no "perfect" hit is found, run closest.
                blen = 3 # Closest only returns the 3 first fields of the original bed.

                # First print empty line for exact match
                line = feat.fields[0:bdlen]+['Overlapping', gene, anot, dist, strand]
                out.write("\t".join(line)+"\n")

                ## Create a mini-bed for the gene and look for
                ## closest 2 genes before and after it (but not overlapping)
                linebed = "\t".join([feat.chrom, str(feat.start), str(feat.stop)])
                gene_bed = pb.BedTool(linebed, from_string=True)

                pregene = gene_bed.closest(ref, D = "ref", id = True, io = True)[0]
                postgene = gene_bed.closest(ref, D = "ref", iu = True, io = True)[0]

                gene, anot, strand = extract_gff_info(pregene, blen)
                dist = str(pregene.fields[blen+9])
                line = feat.fields[0:bdlen]+['Nearest left', gene, anot, dist, strand]

                out.write("\t".join(line)+"\n")
                anot, gene, dist, strand = ".", ".", ".", "."

                gene, anot, strand = extract_gff_info(postgene, blen)
                dist = str(postgene.fields[blen+9])
                line = feat.fields[0:bdlen]+['Nearest right', gene, anot, dist, strand]

                out.write("\t".join(line)+"\n")
                anot, gene, dist, strand = ".", ".", ".", "."


    finally:
        out.close()

bed = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/New_Coverage/New_Peaks/all_peaks_tests.bed"
gff = "/mnt/Disc4T/Projects/PhD_Project/Data/PlDB-46_Pf3D7_GDV1_ncRNAs_500bp_bothends.gff"    

annotate_bed(bed, gff)

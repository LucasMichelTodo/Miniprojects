#!/usr/bin/env nextflow

process Sam_to_bam {
    tag "Sam_to_bam on $samfile"
    publishDir params.outdir+'/Alignments', mode:'copy'
    
    input:
    path samfile

    output:
    path("*.bam")

    script:
    """
    samtools view -bS $samfile > ${samfile.baseName}.bam
    """
}

process Sort_bam {
    tag "Sorting $bamfile"
    publishDir params.outdir+'/Alignments', mode:'copy'
    
    input:
    path bamfile

    output:
    path("*_sort.bam")

    script:
    """
    samtools sort $bamfile > ${bamfile.baseName+'_sort.bam'}
    """
}

process Index_bam {
    tag "Indexing $bamfile"
    publishDir params.outdir+'/Alignments', mode:'copy'
    
    input:
    path bamfile

    output:
    path("*.bam.bai")

    script:
    """
    samtools index $bamfile > ${bamfile.baseName+'.bam.bai'}
    """
}

process Qfilter_bam {
    tag "Filtering by Qval on $bamfile"
    publishDir params.outdir+'/Alignments', mode:'copy'
    
    input:
    path bamfile

    output:
    path("*_q5.bam")

    script:
    """
    samtools view -b -q $params.qval_filter $bamfile > ${bamfile.baseName+'_q5.bam'}
    """
}
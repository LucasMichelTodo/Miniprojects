#!/usr/bin/env nextflow

/*
 * pipeline input parameters
 */
 
params.reads = "$projectDir/Raw_Data/*{R1,R2}*.fq"
params.bbduk_ref = "/home/lucas/Programs/bbmap/resources/adapters.fa"
params.align_index = "/home/lucas/Programs/bowtie2-2.3.0-legacy/Pf3D7"
params.align_params = "--very-sensitive --local -5 4 -3 4 -I 50 -X 2000"
params.qval_filter = 5 
params.outdir = "Results"

log.info """\


    C H I P - S E Q - N F   P I P E L I N E
    =======================================
    reads              : ${params.reads}
    bbduk reference    : ${params.bbduk_ref} 
    bowtie2 index      : ${params.align_index}
    bowtie2 params     : ${params.align_params}
    Qval filter (BAM)  : ${params.qval_filter}
    outdir             : ${params.outdir}


    """
    .stripIndent()


process FastQC_Raw {
    tag "FASTQC on $reads"
    publishDir params.outdir+'/FastQC_RawReads', mode:'copy'

    input:
    path reads

    output:
    path "*"

    script:
    """
    fastqc -f fastq -q ${reads}
    """
}

process FastQC_Clean {
    tag "FASTQC on $reads"
    publishDir params.outdir+'/FastQC_CleanReads', mode:'copy'

    input:
    path reads

    output:
    path "*"

    script:
    """
    fastqc -f fastq -q ${reads}
    """
}

process CleanReads {
	tag "BBDUK on $sample_id"
    publishDir params.outdir+'/CleanReads', mode:'copy'
	
	input:
    tuple val(sample_id), path(reads)

   	output:
    tuple val(sample_id), path("*_cleanReads.fq")

    script:
    """
    bbduk.sh in=${reads[0]} in2=${reads[1]}\
    out=${reads[0].baseName+"_cleanReads.fq"} out2=${reads[1].baseName+"_cleanReads.fq"}\
    outm=${sample_id+"_badReads.fq"}\
    ktrim=r k=22 mink=6 overwrite=t\
    ref=$params.bbduk_ref
    """
}

process Fastq_Screen {
    tag "FASTQ_SCREEN on $reads"
    publishDir params.outdir+'/Fastq_Screens', mode:'copy'

    input:
    path reads

    output:
    path "*"

    script:
    """
    fastq_screen ${reads}
    """
}

process Bowtie2_Align {
    tag "BOWTIE2 on $sample_id"
    publishDir params.outdir+'/Alignments', mode:'copy'
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.sam"), stdout

    script:
    """
    bowtie2 -1 ${reads[0]} -2 ${reads[1]}\
    -S ${sample_id+".sam"} $params.align_params\
    -x $params.align_index\
    --met-file ${sample_id+".metrics"}\
    2>&1
    """
}

process Bwt2_Report {
    tag "Collecting Bwt2 report for $sample_id"
    publishDir params.outdir+'/Alignments/Reports', mode:'copy'

    input:
    tuple val(sample_id), path(alignment), val(report)

    output:
    path "*.log"

    script:
    """
    echo \$'$sample_id\n\n$report' > ${sample_id+'_align.log'}
    """
}

process JoinBwt2reports {
    publishDir params.outdir+'/Alignments/Reports', mode:'copy'

    input:
    path reports

    output:
    path "join_report.log"

    script:
    """
    cat $reports > join_report.log
    """
}

process ParseBwt2reports {
    publishDir params.outdir+'/Alignments/Reports', mode:'copy'
    
    input:
    path join_report

    output:
    path '*_report.tsv'

    script:
    """
    parse_bwt2_reports.py -r $join_report -sl '^[0-9]+ reads; of these:' -o agregated_report.tsv
    """
}

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


workflow {
	Channel.fromFilePairs(params.reads, checkIfExists: true).set{read_ch}
	FastQC_Raw(read_ch.map{id, files -> files}.flatten())
    cl_read_ch = CleanReads(read_ch)
    FastQC_Clean(cl_read_ch.map{id, files -> files}.flatten())
    Fastq_Screen(cl_read_ch.map{id, files -> files}.flatten())
    align_ch = Bowtie2_Align(cl_read_ch)
    reports_ch = Bwt2_Report(align_ch)
    joinreport_chanel = JoinBwt2reports(reports_ch.collect())
    ParseBwt2reports(joinreport_chanel)
    bam_ch1 = Sam_to_bam(align_ch.map{id, files, report -> files})
    bam_ch2 = Sort_bam(bam_ch1)
    bam_ch3 = Index_bam(bam_ch2)
    bam_ch4 = Qfilter_bam(bam_ch2)
    bam_ch5 = Index_bam(bam_ch4)
}


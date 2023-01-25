#!/usr/bin/env nextflow

params.samples_csv = "$projectDir/samples_csv.csv"
params.outdir = "$projectDir/Results_2"

process FastQC_Raw {
    tag "FASTQC on $reads"
    publishDir params.outdir+'/FastQC_RawReads', mode:'copy'

    input:
    tuple val(sample), val(rep), val(ab), path(r1), path(r2)

    output:
    tupple val(sample_id), path("*")

    script:
    """
    fastqc -f fastq -q $r1
    """
}

workflow {
    Channel
    .fromPath( params.samples_csv )
    .splitCsv( header: true, sep: ',' )
    .map { row -> tuple(row.sample, row.replicate, row.antibody, file(row.reads1), file(row.reads2 )) }
    .set { sample_ch }
    //sample_ch.collect().view()
    sample_ch.flatMap { it -> [([it[0], it[1], it[2]].join('_'), it[3]), ([it[0], it[1], it[2]].join('_'):it[4]) ] }
    .view()

    //results = sample_ch.flatMap{ sample, rep, ab, r1, r2 -> [sample, rep, ab].join('_'), r1}   
    //results.subscribe onNext: { println it }, onComplete: { println 'Done' }
    //FastQC_Raw(sample_ch.map{ sample, rep, ab, r1, r2 -> ([sample, rep, ab].join('_'), r1), ([sample, rep, ab].join('_'), r2)})
}

#!/bin/bash -ue
mkdir -p FastQC_RawReads
fastqc -o FastQC_RawReads -f fastq -q NCV28_me_R2_subset1000.fq

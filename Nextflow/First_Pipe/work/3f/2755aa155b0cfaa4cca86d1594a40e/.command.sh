#!/bin/bash -ue
mkdir -p FastQC_CleanReads
fastqc -o FastQC_CleanReads -f fastq -q NCV29_in_R1_subset1000_cleanReads.fq

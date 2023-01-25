#!/bin/bash -ue
mkdir -p FastQC_RawReads
fastqc -o FastQC_RawReads -f fastq -q NCV27_in_R1_subset1000.fq

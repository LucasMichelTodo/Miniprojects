#!/bin/bash -ue
mkdir -p Fastq_Screens
fastqc --outdir Fastq_Screens NCV26_me_R2_subset1000_cleanReads.fq

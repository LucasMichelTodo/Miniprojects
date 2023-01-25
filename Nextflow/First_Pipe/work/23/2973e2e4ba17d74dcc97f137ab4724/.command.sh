#!/bin/bash -ue
mkdir -p Alignments/Metrics
bowtie2 -1 NCV26_me_R1_subset1000_cleanReads.fq -2 NCV26_me_R2_subset1000_cleanReads.fq    -S Alignments/NCV26_me.sam --very-sensitive --local -5 4 -3 4 -I 50 -X 2000    -x /home/lucas/Programs/bowtie2-2.3.0-legacy/Pf3D7    --met-file Alignments/Metrics/NCV26_me.log

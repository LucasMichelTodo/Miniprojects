#!/bin/bash -ue
mkdir -p Alignments/Metrics
mkdir -p Alignments/Reports
bowtie2 -1 NCV29_in_R1_subset1000_cleanReads.fq -2 NCV29_in_R2_subset1000_cleanReads.fq    -S Alignments/NCV29_in.sam --very-sensitive --local -5 4 -3 4 -I 50 -X 2000    -x /home/lucas/Programs/bowtie2-2.3.0-legacy/Pf3D7    --met-file Alignments/Metrics/NCV29_in.metrics    > Alignments/Reports/NCV29_in.log

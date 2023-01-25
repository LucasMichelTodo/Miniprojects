#!/bin/bash -ue
mkdir -p Alignments/Metrics
mkdir -p Alignments/Reports
bowtie2 -1 NCV27_in_R1_subset1000_cleanReads.fq -2 NCV27_in_R2_subset1000_cleanReads.fq    -S Alignments/NCV27_in.sam --very-sensitive --local -5 4 -3 4 -I 50 -X 2000    -x /home/lucas/Programs/bowtie2-2.3.0-legacy/Pf3D7    --met-file Alignments/Metrics/NCV27_in.metrics    2>&1 > Alignments/Reports/NCV27_in.log

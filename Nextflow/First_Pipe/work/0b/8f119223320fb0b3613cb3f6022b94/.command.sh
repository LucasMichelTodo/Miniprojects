#!/bin/bash -ue
mkdir -p Alignments
    bowtie2 -1 NCV27_in_R1_subset1000_cleanReads.fq -2 NCV27_in_R2_subset1000_cleanReads.fq -S Alignments/NCV27_in.sam 
--very-sensitive --local \n
-5 4 -3 4 -I 50 -X 2000 \n
-x /home/lucas/Programs/bowtie2-2.3.0-legacy/Pf3D7 \n

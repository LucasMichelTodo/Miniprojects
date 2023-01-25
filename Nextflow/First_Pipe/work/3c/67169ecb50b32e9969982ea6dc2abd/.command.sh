#!/bin/bash -ue
mkdir -p Alignments
samtools sort NCV28_me.bam > Alignments/NCV28_me_sort.bam

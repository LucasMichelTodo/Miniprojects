#!/bin/bash -ue
mkdir -p Alignments
samtools sort NCV29_in.bam > Alignments/NCV29_in_sort.bam

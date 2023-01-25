#!/bin/bash -ue
mkdir -p Alignments
samtools sort NCV27_in.bam > Alignments/NCV27_in_sort.bam

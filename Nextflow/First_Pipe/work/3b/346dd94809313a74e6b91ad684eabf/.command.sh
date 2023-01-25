#!/bin/bash -ue
mkdir -p Alignments
samtools sort NCV26_me.bam > Alignments/NCV26_me_sort.bam

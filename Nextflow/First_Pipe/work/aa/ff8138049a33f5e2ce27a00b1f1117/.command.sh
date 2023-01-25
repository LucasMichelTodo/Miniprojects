#!/bin/bash -ue
mkdir -p Alignments
samtools view -bS NCV29_in.sam > Alignments/NCV29_in.bam

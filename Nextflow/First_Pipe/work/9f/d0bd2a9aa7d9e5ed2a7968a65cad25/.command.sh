#!/bin/bash -ue
mkdir -p Alignments
samtools view -bS NCV27_in.sam > Alignments/NCV27_in.bam

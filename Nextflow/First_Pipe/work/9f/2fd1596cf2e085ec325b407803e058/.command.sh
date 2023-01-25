#!/bin/bash -ue
mkdir -p Alignments
samtools view -bS NCV28_me.sam > Alignments/NCV28_me.bam

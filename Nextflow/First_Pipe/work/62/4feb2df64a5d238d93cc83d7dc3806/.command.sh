#!/bin/bash -ue
mkdir -p Alignments
samtools view -bS NCV26_me.sam > Alignments/NCV26_me.bam

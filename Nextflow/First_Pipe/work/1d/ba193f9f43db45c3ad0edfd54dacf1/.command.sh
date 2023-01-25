#!/bin/bash -ue
mkdir -p CleanReads
bbduk.sh in=NCV29_in_R1_subset1000.fq in2=NCV29_in_R2_subset1000.fq    out=CleanReads/NCV29_in_R1_subset1000_cleanReads.fq out2=CleanReads/NCV29_in_R2_subset1000_cleanReads.fq    outm=CleanReads/NCV29_in_badReads.fq    ktrim=r k=22 mink=6 overwrite=t    ref=/home/lucas/Programs/bbmap/resources/adapters.fa

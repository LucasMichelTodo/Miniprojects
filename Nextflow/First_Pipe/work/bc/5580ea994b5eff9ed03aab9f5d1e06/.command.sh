#!/bin/bash -ue
bbduk.sh in=NCV26_me_R1_subset1000.fq in2=NCV26_me_R2_subset1000.fq    out=NCV26_me_R1_subset1000_cleanReads.fq out2=NCV26_me_R2_subset1000_cleanReads.fq    outm=NCV26_me_badReads.fq    ktrim=r k=22 mink=6 overwrite=t    ref=/home/lucas/Programs/bbmap/resources/adapters.fa

#!/bin/bash -ue
mkdir -p CleanReads
bbduck.sh in=NCV28_me_R1_subset1000.fq in2=NCV28_me_R2_subset1000.fq     out=CleanReads/NCV28_me_R1_subset1000_cleanReads.fq out2=CleanReads/NCV28_me_R2_subset1000_cleanReads.fq outm=CleanReads/NCV28_mebadReads.fq     ktrim=r k=22 mink=6 overwrite=t

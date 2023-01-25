#!/bin/bash -ue
bamCoverage -b NCV27_in_sort_q5_noDup.bam -of bedgraph --normalizeUsing RPKM -bs 10 --smoothLength 200 -o NCV27_in_sort_q5_noDup_rpkm.bdg

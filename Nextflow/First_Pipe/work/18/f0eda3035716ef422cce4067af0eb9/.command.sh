#!/bin/bash -ue
bamCoverage -b NCV26_me_sort_q5_noDup.bam -of bedgraph --normalizeUsing RPKM -bs 10 --smoothLength 200 -o NCV26_me_sort_q5_noDup_rpkm.bdg

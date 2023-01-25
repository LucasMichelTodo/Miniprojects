#!/bin/bash -ue
bamCoverage -b NCV29_in_sort_q5_noDup.bam --outFileFormat bedgraph --normalizeUsing RPKM -bs 10 --smoothLength 200 -o NCV29_in_sort_q5_noDup_rpkm.bdg

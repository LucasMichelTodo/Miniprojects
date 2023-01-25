#!/bin/bash -ue
mkdir -p Alignments/Reports
cat NCV26_me_align.log NCV28_me_align.log NCV29_in_align.log NCV27_in_align.log > Alignments/Reports/join_report.log

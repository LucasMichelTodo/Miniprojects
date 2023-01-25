#!/bin/bash -ue
parse_bwt2_reports.py -r join_report.log -sl '^[0-9]+ reads; of these:' -o /mnt/Disc4T/Projects/Miniprojects/Nextflow/First_Pipe/Results/Alignments/Reports/agregated_report.tsv

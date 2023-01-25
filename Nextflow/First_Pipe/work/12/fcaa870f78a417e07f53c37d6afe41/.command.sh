#!/bin/bash -ue
parse_bwt2_reports.py -r join_report.log -sl '^[0-9]+ reads; of these:' -o agregated_report.tsv

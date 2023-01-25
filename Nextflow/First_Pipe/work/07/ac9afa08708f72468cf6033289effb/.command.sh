#!/usr/bin/env python3

import os
from parse_bwt2_reports import parse_agregated_report

linestart = re.compile('^[0-9]+ reads; of these:')

parse_agregated_report('Alignments/Reports/join_report.log', linestart, 'Alignments/Reports/agregated_report.tsv')

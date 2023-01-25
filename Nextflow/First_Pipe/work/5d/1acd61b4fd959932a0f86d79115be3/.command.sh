#!/usr/bin/env python
from parse_bwt2_reports parse_agregated_report

linestart = re.compile('^[0-9]+ reads; of these:')

parse_agregated_report('Alignments/Reports/join_report.log', linestart, 'Alignments/Reports/agregated_report.tsv')

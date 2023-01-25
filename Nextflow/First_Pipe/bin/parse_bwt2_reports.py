#!/usr/bin/env python3

import re
import pandas as pd
import argparse

def parse_agregated_report(report, start_line, out_name):
    
    start_line = re.compile(start_line)

    with open(report, 'r+') as infile:

        infile = open(report, 'r+')
        lines = infile.readlines()
        lines = [line.rstrip() for line in lines]
        out = {}

        for idx, line in enumerate(lines):

            if re.match(start_line, line):

                name = lines[idx-2].strip()
                out[name] = {}
                out[name]['nreads'] = lines[idx].split()[0]
                paired = re.search('\((.*)\)', lines[idx+1])
                out[name]['paired'] = paired.group(1)
                conc_0 = re.search('\((.*)\)', lines[idx+2])
                out[name]['conc_0'] = conc_0.group(1)
                conc_1 = re.search('\((.*)\)', lines[idx+3])
                out[name]['conc_1'] = conc_1.group(1)
                conc_more = re.search('\((.*)\)', lines[idx+4])
                out[name]['conc_more'] = conc_more.group(1)

                conc_0_disc1 = re.search('\((.*)\)', lines[idx+7])
                out[name]['conc_0_disc1'] = conc_0_disc1.group(1)

                conc_0_disc0_al0 = re.search('\((.*)\)', lines[idx+11])
                out[name]['conc_0_disc0_al0'] = conc_0_disc0_al0.group(1)
                conc_0_disc0_al1 = re.search('\((.*)\)', lines[idx+12])
                out[name]['conc_0_disc0_al1'] = conc_0_disc0_al1.group(1)
                conc_0_disc0_almore = re.search('\((.*)\)', lines[idx+13])
                out[name]['conc_0_disc0_almore'] = conc_0_disc0_almore.group(1)

                out[name]['overall_al'] = lines[idx+14].replace(' overall alignment rate', '')
        

        out_df = pd.DataFrame.from_dict(out, orient='index')
        
        # Rename cols
        cols = out_df.columns.tolist()

        new_cols = [
            'Total Reads',
            'Paired Reads',
            'Aligned Concordantly 0 times',
            'Aligned Concordantly 1 time',
            'Aligned Concordantly >1 times',
            'Aligned Concordantly 0 times, Discordantly 1',
            'Non-concordant or discordant, Mate Aligned 0 times',
            'Non-concordant or discordant, Mate Aligned 1 time',
            'Non-concordant or discordant, Mate Aligned >1 times',
            'Overall Alignment Rate',
        ]

        re_cols = dict(zip(cols, new_cols))
        out_df.rename(columns=re_cols, inplace=True)
        out_df.to_csv(out_name, sep='\t')

parser = argparse.ArgumentParser()
parser.add_argument('-r','--report', required=True)
parser.add_argument('-sl','--startline', required=True)
parser.add_argument('-o','--outname', required=True)
args = parser.parse_args()

linestart = args.startline
report_file = args.report
out_file = args.outname

parse_agregated_report(report_file, linestart, out_file)
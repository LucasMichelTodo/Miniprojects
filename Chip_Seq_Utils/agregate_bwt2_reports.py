import re
import pandas as pd
import os

def parse_agregated_report(report, start_line, out_name):
    with open(report, 'r+') as infile:

        infile = open(report, 'r+')
        lines = infile.readlines()
        lines = [line.rstrip() for line in lines]
        out = {}

        for idx, line in enumerate(lines):

            if line.startswith(start_line):

                name = line.split('/')[-1].strip()
                out[name] = {}

                out[name]['call'] = lines[idx+2]
                out[name]['nreads'] = lines[idx+3].split()[0]

                paired = re.search('\((.*)\)', lines[idx+4])
                out[name]['paired'] = paired.group(1)
                conc_0 = re.search('\((.*)\)', lines[idx+5])
                out[name]['conc_0'] = conc_0.group(1)
                conc_1 = re.search('\((.*)\)', lines[idx+6])
                out[name]['conc_1'] = conc_1.group(1)
                conc_more = re.search('\((.*)\)', lines[idx+7])
                out[name]['conc_more'] = conc_more.group(1)

                conc_0_disc1 = re.search('\((.*)\)', lines[idx+10])
                out[name]['conc_0_disc1'] = conc_0_disc1.group(1)

                conc_0_disc0_al0 = re.search('\((.*)\)', lines[idx+14])
                out[name]['conc_0_disc0_al0'] = conc_0_disc0_al0.group(1)
                conc_0_disc0_al1 = re.search('\((.*)\)', lines[idx+15])
                out[name]['conc_0_disc0_al1'] = conc_0_disc0_al1.group(1)
                conc_0_disc0_almore = re.search('\((.*)\)', lines[idx+16])
                out[name]['conc_0_disc0_almore'] = conc_0_disc0_almore.group(1)

                out[name]['overall_al'] = lines[idx+17].replace(' overall alignment rate', '')

        out_df = pd.DataFrame.from_dict(out, orient='index')
        # Reorder cols
        cols = out_df.columns.tolist()
        cols = [cols[1]]+[cols[-1]]+cols[2:]+[cols[0]]
        out_df = out_df[cols]

        new_cols = [
            'Total Reads',
            'Overall Alignment Rate',
            'Paired Reads',
            'Aligned Concordantly 0 times',
            'Aligned Concordantly 1 time',
            'Aligned Concordantly >1 times',
            'Aligned Concordantly 0 times, Discordantly 1',
            'Non-concordant or discordant, Mate Aligned 0 times',
            'Non-concordant or discordant, Mate Aligned 1 time',
            'Non-concordant or discordant, Mate Aligned >1 times',
        ]

        re_cols = dict(zip(cols, new_cols))
        out_df.rename(columns=re_cols, inplace=True)
        out_df.to_csv(out_name, sep='\t')

wd = '/mnt/Disc4T/Projects/Checks_on_ChIPs/Test_Reports/'
os.chdir(wd)

report = 'join_report.txt'
start_line = '/mnt/Disc4T/Projects/Checks_on_ChIPs/Test_Alignments/NCV'

parse_agregated_report(report, start_line, 'join_report_stats.tsv')

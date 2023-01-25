import subprocess as sp
import re
import pandas as pd
from datetime import datetime
import os

## Funtions

def call_Bowtie2(in1, in2, out, params):
    cmd = ['bowtie2', '-1', in1, '-2', in2, '-S', out]
    params = params.split(' ')
    rep_name = out.rsplit('.', 1)[0]+'_bwt2_align_report.txt'

    cmd = cmd+params
    print(cmd)
    with open(rep_name, 'w+') as outfile:
        outfile.write(f'{out}\n\n'+' '.join(cmd)+'\n')

    result = sp.run(cmd, stdout=sp.PIPE, stderr=sp.STDOUT)
    str_result = result.stdout.decode('utf-8')

    with open(rep_name, 'a+') as outfile:
        outfile.write(str_result)


def from_sam_to_bam(samfile):
    
    name = samfile.rsplit(".", 1)[0]
    cmd = "samtools view -bS {} > {}" .format(samfile, name+".bam")
    print('\n'+cmd+'\n')
    sp.call(cmd, shell=True)

    ### Erase SAM after creating BAM
    cmd = "rm {}" .format(samfile)
    print('\n'+cmd+'\n')
    sp.call(cmd, shell=True)

    cmd = "samtools sort {} > {}" .format(name+".bam", name+"_sort.bam")
    print('\n'+cmd+'\n')
    sp.call(cmd, shell=True)

    ### Erase bam after creating sortedBAM
    cmd = "rm {}" .format(name+".bam")
    print('\n'+cmd+'\n')
    sp.call(cmd, shell=True)

    cmd = "samtools index {} > {}" .format(name+"_sort.bam", name+"_sort.bam.bai")
    print('\n'+cmd+'\n')
    sp.call(cmd, shell=True)

    ## Filter only >=q5 reads
    cmd = "samtools view -b -q 5 {} > {}" .format(name+"_sort.bam", name+"_q5_sort.bam")
    print('\n'+cmd+'\n')
    sp.call(cmd, shell=True)

    cmd = "samtools index {} > {}" .format(name+"_q5_sort.bam", name+"_q5_sort.bam.bai")
    print('\n'+cmd+'\n')
    sp.call(cmd, shell=True)

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
        
def common_start(sa, sb):
    """ returns the longest common substring from the beginning of sa and sb """
    def _iter():
        for a, b in zip(sa, sb):
            if a == b:
                yield a
            else:
                return

    return ''.join(_iter())       

def bowtie2_align_pipe(read1s_list, read2s_list, params, outpath = './'):
    """ 
    Take a sorted list of read1 and read 2 files (both lists must be in the same sort order),
    a params string with all wanted bowtie2 parameters and optionally a path for the outputs.
    Runs the alignments and generates an automated report.
    """
    startTime = datetime.now()
    os.makedirs(outpath, exist_ok=True)
    
    for r1, r2 in zip(read1s_list, read2s_list):
        
        name = common_start(os.path.basename(r1), os.path.basename(r2))
        out = outpath+name+'.sam'
        print('\n'+f'Aligning: {out}'+'\n')
        call_Bowtie2(r1, r2, out, params)
        from_sam_to_bam(out)

    finishTime = datetime.now()
    time_lines = (
        f'\nProcess started at {startTime}\n'
        f'Process ended at {finishTime}\n'
        f'Process started at {finishTime - startTime}'
    )

    ## Write report
    reports = [f for f in os.listdir(outpath) if f.endswith('bwt2_align_report.txt')]
    with open('full_report.txt', 'w+') as outfile:
        for fname in reports:
            with open(outpath+fname) as infile:
                outfile.write(infile.read())
        outfile.write(time_lines)

    ## Generate table report
    parse_agregated_report('full_report.txt', outpath, 'full_report_table.tsv')




import os
import sys
import subprocess as sp
import pandas as pd
import re
from datetime import datetime

#### FUNCTIONS ####

## Pre-Processing

def call_BBDUK(in1, in2, ref, outfolder='./', params=''):

    if not outfolder.endswith('/'): outfolder = outfolder+'/'

    extension = '.'+in1.rsplit('.', 1)[1]

    out1 = outfolder + in1.rsplit("/", 1)[1].replace(extension, "_clean.fq")
    out2 = outfolder + in2.rsplit("/", 1)[1].replace(extension, "_clean.fq")
    outm = outfolder + in1.rsplit("/", 1)[1].replace(extension, "_badreads.fq")

    cmd = [
        'bbduk.sh',
        'in='+in1, 'in2='+in2,
        'out='+out1, 'out2='+out2, 'outm='+outm,
        'ref='+ref,
    ]

    params = params.split(' ')
    cmd = cmd+params
    sp.run(cmd)

def call_fastq_screen(fastq_file_list, threads = '4', outdir = './'):
    ## input a list of fastq files to be screened
    cmd = ['fastq_screen', '--threads', threads, '--outdir', outdir] + fastq_file_list
    sp.run(cmd)

def call_fastqc(fastq_file_list, threads = '4', outdir = './'):
    ## input a list of fastq files to be screened
    cmd = ['fastqc', '--threads', threads, '--outdir', outdir] + fastq_file_list
    sp.run(cmd)

## Alignment
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
        name1 = r1.rsplit('/', 1)[1].rsplit('.', 1)[0]
        name2 = r2.rsplit('/', 1)[1].rsplit('.', 1)[0]
        name = common_start(name1, name2)
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


## Samtools
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

## Remove Duplicates

def remove_duplicates(indir, outdir, bam):

    gatk_path = '/home/lucas/Programs/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar'
    i = indir+bam
    o = outdir+bam.replace(".bam", "_noDup.bam")
    m = outdir+bam.replace(".bam", "_metrics.txt")

    cmd = (f"java -jar {gatk_path} MarkDuplicates "
           "REMOVE_DUPLICATES=true I={} O={} M={}") .format(i, o, m)

    sp.call(cmd, shell=True)

## Get coverage (Tracks)

def get_RPKMs(bam, bs, smooth, norm, outfld):
    
    outfld = outfld if outfld.endswith('/') else outfld + '/'
    name = bam.rsplit('/', 1)[1]
    outname = outfld+name.replace('.bam', f'_rpkm_raw_bs{bs}_smth{smooth}.bdg')

    cmd = (
        f'bamCoverage -b {bam} '
        '--outFileFormat bedgraph '
        f'--normalizeUsing {norm} '
        '-p 8 '
        f'-bs {bs} '
        f'--smoothLength {smooth} '
        f'-o {outname}'
    )
    print(cmd)
    subprocess.call(cmd, shell=True)


def get_RPKMs_normInput(
        bam_IP,
        bam_in,
        bs = 50,
        smooth = 150,
        norm = 'RPKM',
        of = 'bedgraph',
        outfld = './',
        pseudo = 1,
        num_process = 8,
):
    outfld = outfld if outfld.endswith('/') else outfld + '/'
    name = bam_IP.rsplit('/', 1)[1]
    outname = name.replace('.bam', f'_rpkm_normInput_bs{bs}_smth{smooth}_pseudo{pseudo}.bdg')
    cmd = (
        f'bamCompare -b1 {bam_IP} -b2 {bam_in} '
        '--outFileFormat bedgraph '
        '--scaleFactorsMethod None '
        f'--normalizeUsing {norm} '
        f'-p {num_process} '
        f'-bs {bs} '
        f'--smoothLength {smooth} '
        f'--pseudocount {pseudo} '
        f'-o {outfld+outname} '
        f'-of {of}'
    )

    print(cmd)
    subprocess.call(cmd, shell=True)

def params_to_name(params_dict):
    """
    Take a dicttionary with params and convert it to a string to add to name files.
    """
    suffix = ''
    for k,v in params_dict.items():
        k_str = k.replace('-', '').strip()
        v_str = v.strip()
        suffix += f'_{k_str}{v_str}'
    return(suffix)

def input_parser(input_file):
    with open(input_file, 'r+') as params:
        for line in params:
            


#### Parse Arguments and Run ####

bbduk_ref = '/home/lucas/Programs/bbmap/resources/adapters.fa'


## Input Data

## wd = '/mnt/Disc4T/Projects/Miniprojects/Chip_Seq_Pipe/'
## os.chdir(wd)

raw_data_dir = './Raw_Data/'
samples = pd.read_excel('./chip_seq_samples.xlsx', engine='openpyxl', skiprows=1)

read1s = samples['Read 1'].tolist() + samples['Read 1.1'].tolist()
read2s = samples['Read 2'].tolist() + samples['Read 2.1'].tolist()

print('Cleaning Reads... ')
for f in read1s: print(f)
print('-----')
for f in read2s: print(f)
print('-----')

## Clean Reads

params = "ktrim=r k=22 mink=6 overwrite=t "
outpath = './Clean_Reads/'
os.makedirs(outpath, exist_ok=True)

for pair in zip(read1s, read2s):
    in1, in2 = raw_data_dir+pair[0], raw_data_dir+pair[1]
    call_BBDUK(in1, in2, outpath, params)

## Call fastq_screen

all_reads = [raw_data_dir + f for f in read1s+read2s]
outdir = './Fastq_Screens/'
os.makedirs(outdir, exist_ok=True)
call_fastq_screen(all_reads, threads = '8', outdir=outdir)

## Call fastqc on clean reads
print('Calling FastQC on clean reads...\n\n')

clean_reads_dir = './Clean_Reads/'
all_reads = [clean_reads_dir+f for f in os.listdir(clean_reads_dir) if f.endswith('clean.fq')]
outdir = './FastQC_Clean/'
os.makedirs(outdir, exist_ok=True)
call_fastqc(all_reads, threads = '8', outdir = outdir)

## Call fastq_screen on clean reads
print('Calling Fastq_Screen on clean reads...\n\n')

outdir = './Fastq_Screens_Clean/'
os.makedirs(outdir, exist_ok=True)
call_fastq_screen(all_reads, threads = '8', outdir=outdir)


## Algin Clean Reads
print('Starting Alignments...\n\n')

indir = './Clean_Reads/'
outdir = './Alignments/'
os.makedirs(outdir, exist_ok=True)

all_reads = sorted([f for f in os.listdir(indir) if f.endswith('clean.fq')])
reads1 = [indir+r for r in all_reads if '_R1_' in r or '_read1' in r]
reads2 = [indir+r for r in all_reads if '_R2_' in r or '_read2' in r]

params = ("-p 8 --very-sensitive --local "
          "-5 4 -3 4 -I 50 -X 2000 "
          "-x /home/lucas/Programs/bowtie2-2.3.0-legacy/Pf3D7")

bowtie2_align_pipe(reads1, reads2, params, outdir)

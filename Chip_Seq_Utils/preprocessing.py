import subprocess as sp
import os

## Functions

def call_BBDUK(in1, in2, outfolder='./', params=''):

    ref = '/home/lucas/Programs/bbmap/resources/adapters.fa'
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

    
    

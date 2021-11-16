import subprocess as sp
import multiprocessing as mp

### Base functions

def get_nreads(bam):
    cmd = ['samtools', 'view', '-c', bam]
    result = sp.run(cmd, stdout=sp.PIPE)
    out = int(result.stdout.decode('utf-8').strip())
    return out

def get_nreads_mapped(bam):
    ## original cmd: samtools view -F 0x4 foo.sorted.bam | cut -f 1 | sort | uniq | wc -l
    cmd = ['samtools', 'view', '-F', '0x4', bam]
    result = sp.run(cmd, stdout=sp.PIPE).stdout.decode('utf-8')
    read_list = result.strip().split('\n')
    id_list = sorted([read.split('\t')[0] for read in read_list])
    n_unique = len(set(id_list))
    return n_unique

def calc_percentage(reads_mapped, reads_total, read_type = 'paired'):
    factor = {'single':1, 'paired':2}
    return (reads_mapped*factor[read_type]*100)/reads_total

### Multiprocessing

def parllalel_get_align_stats(bam_list, out_file):
    ## Get data
    pool = mp.Pool(mp.cpu_count())
    nreads = pool.map(get_nreads, bam_list)
    mappeds = pool.map(get_nreads_mapped, bam_list)
    percs = list(map(calc_percentage, mappeds, nreads))

    str_nreads = [str(x) for x in nreads]
    str_mappeds = [str(x) for x in mappeds]
    str_percs = [str(x) for x in percs]

    ## Generate output table
    colnames = ['Samples', 'Total Reads', 'Mapped Alignments', 'Alignment Rate (%)']
    with open(out_file, 'w+') as out:
        out.write('\t'.join(colnames)+'\n')
        for x in zip(bam_list, str_nreads, str_mappeds, str_percs):
            out.write('\t'.join(x)+'\n')

# if __name__ == "__main__":
#   parllalel_get_align_stats()

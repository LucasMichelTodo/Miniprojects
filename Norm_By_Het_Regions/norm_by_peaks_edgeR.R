library(edgeR)
library(Rsubread)
library(tidyverse)

wd <- '/mnt/Disc4T/Projects/Miniprojects/Norm_By_Het_Regions/'
setwd(wd)

## Load Common Peaks File 
#peaks_file <- '/mnt/Disc4T/Projects/Miniprojects/Norm_By_Het_Regions/common_peaks_ALL_me_files.bed'
peaks_file <- './join_peaks_sorted_merged.bed'
peaks_df <- read_tsv(peaks_file, col_names = c('Chr', 'Start', 'End'))

## Load whole genome binned for inputs
bins_file <- './whole_genome__binned_10000.bed'
bins_df <- read_tsv(bins_file, col_names = c('Chr', 'Start', 'End'))

## Convert to SAF (single annotation format)
peaks_saf <- peaks_df %>%
  select(Chr, Start, End) %>%
  mutate(Strand = '+') %>%
  mutate(GeneID = row_number())

bins_saf <- bins_df %>%
  select(Chr, Start, End) %>%
  mutate(Strand = '+') %>%
  mutate(GeneID = row_number())

## Load alignment data (bam files)
bamdir <- '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/Bams/'
bamfiles <- paste0(bamdir, list.files(bamdir, pattern = '.*_me_.*_q5.bam$'))
infiles <- paste0(bamdir, list.files(bamdir, pattern = '.*_in_.*_q5.bam$'))

## Call featureCounts to get counts matrix
feat_counts <- featureCounts(files = bamfiles, annot.ext = peaks_saf, isPairedEnd = T, requireBothEndsMapped = T, nthreads = 8)
raw.counts <- feat_counts$counts

feat_counts_in <- featureCounts(files = infiles, annot.ext = bins_saf, isPairedEnd = T, requireBothEndsMapped = T, nthreads = 8)
raw.counts.in <- feat_counts_in$counts

## Use edgeR to calculate normalizing factors
NormFactor <- calcNormFactors(object = raw.counts, method = "TMM")
NormFactor.rle <- calcNormFactors(object = raw.counts, method = "RLE")

NormFactor.in <- calcNormFactors(object = raw.counts.in, method = "TMM")
NormFactor.rle.in <- calcNormFactors(object = raw.counts.in, method = "RLE")

## Raw library size:
LibSize <- colSums(raw.counts)
LibSize.in <- colSums(raw.counts.in)

## Calculate size factors by which we shlould divide the raw coverage:
SizeFactors <- NormFactor * LibSize / 1000000
SizeFactors.in <- NormFactor.in * LibSize.in / 1000000

## Reciprocal (since we will be using deepTools bamCoverage --scaleFactors which multiply instead of dividing the raw coverage)
SizeFactors.Reciprocal <- 1/SizeFactors
SizeFactors.Reciprocal.in <- 1/SizeFactors.in

## Create output table
tibble(Files = bamfiles, Factors = SizeFactors.Reciprocal) %>%
  write_tsv('bamfiles_and_factors_join_peaks.tsv')

tibble(Files = infiles, Factors = SizeFactors.Reciprocal.in) %>%
  write_tsv('bamfiles_and_factors_inputs_binned10000.tsv')


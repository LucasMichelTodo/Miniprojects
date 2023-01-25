#!/bin/bash -euo pipefail
echo "blastn  -num_threads 8 -db testData/headtest -query headtest.3.fasta -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen frames salltitles qcovs' -evalue 1e-3 -out blastout" > blast.log
blastn  -num_threads 8 -db testData/headtest -query headtest.3.fasta -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen frames salltitles qcovs' -evalue 1e-3 -out blastout

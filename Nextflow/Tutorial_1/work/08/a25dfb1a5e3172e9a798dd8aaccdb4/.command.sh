#!/bin/bash -euo pipefail
makeblastdb -in /home/lucas/.nextflow/assets/isugifNF/blast/testData/headtest.fasta -dbtype 'nucl' -out testData/headtest
# makeblastdb -in /home/lucas/.nextflow/assets/isugifNF/blast/testData/headtest.fasta -dbtype 'prot' -out testData/headtest

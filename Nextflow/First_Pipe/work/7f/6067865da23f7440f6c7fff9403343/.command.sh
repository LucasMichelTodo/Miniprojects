#!/bin/bash -ue
java -jar /home/lucas/Programs/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar MarkDuplicates REMOVE_DUPLICATES=true    I=NCV27_in_sort_q5.bam    O=NCV27_in_sort_q5_noDup.bam    M=NCV27_in_sort_q5_metrics.txt

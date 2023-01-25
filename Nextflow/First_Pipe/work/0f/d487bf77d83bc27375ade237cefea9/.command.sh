#!/bin/bash -ue
mkdir -p Alignments/Reports
    echo $'NCV26_me
989 reads; of these:
  989 (100.00%) were paired; of these:
    694 (70.17%) aligned concordantly 0 times
    106 (10.72%) aligned concordantly exactly 1 time
    189 (19.11%) aligned concordantly >1 times
    ----
    694 pairs aligned concordantly 0 times; of these:
      1 (0.14%) aligned discordantly 1 time
    ----
    693 pairs aligned 0 times concordantly or discordantly; of these:
      1386 mates make up the pairs; of these:
        1340 (96.68%) aligned 0 times
        9 (0.65%) aligned exactly 1 time
        37 (2.67%) aligned >1 times
32.25% overall alignment rate
' > Alignments/Reports/NCV26_me_align.log

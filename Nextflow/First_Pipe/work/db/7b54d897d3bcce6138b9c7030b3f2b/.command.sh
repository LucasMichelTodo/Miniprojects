#!/bin/bash -ue
mkdir -p Alignments/Reports
    echo $'NCV28_me

915 reads; of these:
  915 (100.00%) were paired; of these:
    300 (32.79%) aligned concordantly 0 times
    383 (41.86%) aligned concordantly exactly 1 time
    232 (25.36%) aligned concordantly >1 times
    ----
    300 pairs aligned concordantly 0 times; of these:
      9 (3.00%) aligned discordantly 1 time
    ----
    291 pairs aligned 0 times concordantly or discordantly; of these:
      582 mates make up the pairs; of these:
        560 (96.22%) aligned 0 times
        2 (0.34%) aligned exactly 1 time
        20 (3.44%) aligned >1 times
69.40% overall alignment rate
' > Alignments/Reports/NCV28_me_align.log

#!/bin/bash -euo pipefail
echo "blastn -version" > software_check.txt
  blastn -version >> software_check.txt

  echo "
makeblastdb -version" >> software_check.txt
  makeblastdb -version >> software_check.txt

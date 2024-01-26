#!/usr/bin/bash
# USAGE: ./process_SRR_list.sh SRR_list.txt

if [ "$#" -ne 1 ]; then
  echo "Usage: ./get_SRR_data.sh SRR_list.txt"
  exit 1
fi


while IFS= read -r SRR_ID; do
  wget https://sra-pub-run-odp.s3.amazonaws.com/sra/$SRR_ID/$SRR_ID -O $SRR_ID.sra
  fastq-dump --gzip --split-3 $SRR_ID.sra
done < "$1"
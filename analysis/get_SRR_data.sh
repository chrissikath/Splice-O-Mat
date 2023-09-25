#!/usr/bin/bash

wget https://sra-pub-run-odp.s3.amazonaws.com/sra/$1/$1 -O $1.sra
fastq-dump --gzip --split-3 $1.sra
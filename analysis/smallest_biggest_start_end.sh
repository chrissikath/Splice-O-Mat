#!/bin/bash

gtf_file="CELSR1_exons.gtf"

smallest_start=$(awk '$3 == "exon" {print $4}' $gtf_file | sort -n | head -n 1)
greatest_end=$(awk '$3 == "exon" {print $5}' $gtf_file | sort -n | tail -n 1)

echo "Smallest Start Position: $smallest_start"
echo "Greatest End Position: $greatest_end"


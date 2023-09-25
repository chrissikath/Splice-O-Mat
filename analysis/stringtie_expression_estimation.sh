#!/bin/sh

#Author Christina Kuhn @github: https://github.com/chrissikath
# pre-used scripts from Alex 
# Shell Script, welches anhand von RNA-Seq Dateien Splice-Varianten findet Part 2
# STAR version 2.7.4a
# Stringtie version 2.1.3b


############################ 2.1 Stingtie MERGE #################################
#for i in $(cat $1); do
#    echo $i; 
#    OUTPUT="/home/christina/Transcript_Variants_aGPCR/analysis/$i";
#    ls ${OUTPUT}/stringtie/*.gtf >> new_mergefile.txt;
#done 

#inside the mergfile have to be the merge_old.gtf file from the last database/stringite run
#stringtie --merge -G $ANNO -o combo_version2.gtf new_mergefile.txt -l "NSTRG";

# ########################### 2.3 Stingtie expression estimation mode ##############
for i in $(cat $1); do
    RAW_SEQ="/home/christina/Transcript_Variants_aGPCR/data/$i/"
    SAMPLE_LIST="/home/christina/Transcript_Variants_aGPCR/data/$i/SRR_Acc_List.txt"
    OUTPUT="/home/christina/Transcript_Variants_aGPCR/analysis/$i"

    #Genome information
    STARINDEX_OUTPUT="/home/christina/Transcript_Variants_aGPCR/analysis/STARIndex"
    GENOME="hg38.fna"                                                              
    ANNO="/opt/genomes/human/hg38/ncbi_full_analysis_set/hg38.ncbiRefSeq.gtf"

    echo 'stringtie expression estimation mode';
    mkdir -p ${OUTPUT}/ballgown_redone_2/;
    cat ${SAMPLE_LIST} | cut -f2 | sed 's/\r$//' | while read sample; do 
    	echo $sample `date +%T`;
    	stringtie -e -B -p 32 -G combo_version2.gtf \
    	-o ${OUTPUT}/ballgown_redone_2/${sample}/${sample}.gff ${OUTPUT}/star/$sample/Aligned.out.bam_sorted.bam;
    done;
done;

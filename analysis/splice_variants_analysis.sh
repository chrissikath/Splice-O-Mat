#!/bin/sh

#Author Christina Kuhn @github: https://github.com/chrissikath
# pre-used scripts from Alexander Bernd Knierim
# Shell Script, which assembles possible splice-variants from RNA-Seq data Part 1
# STAR version 2.7.4a
# Stringtie version 2.1.3b

#used STAR index from "/opt/genomes/human/hg38/ncbi_full_analysis_set"                                                       
RAW_SEQ="$(pwd)/" # exmaple: RAW_SEQ="/home/christina/Transcript_Variants_aGPCR/data/GSE182321_OpiodUseDisorder/"
COHORT=$(basename $(pwd))
SAMPLE_LIST="$RAW_SEQ/SRR_Acc_List.txt"
OUTPUT="/home/christina/Transcript_Variants_aGPCR/analysis/$COHORT/"

#Genome information
STARINDEX_OUTPUT="/home/christina/Transcript_Variants_aGPCR/analysis/STARIndex"
GENOME="hg38.fna"                                                              
ANNO="/opt/genomes/human/hg38/ncbi_full_analysis_set/hg38.ncbiRefSeq.gtf"

############################# 1.1 STAR INDEX ################################
# creates STAR index into directory ${STARINDEX_OUTPUT}
# !need to run only one time!

# STAR --runThreadN 20 --runMode genomeGenerate \
# --genomeDir ${STARINDEX_OUTPUT} \
# -- genomeFastaFiles $GENOME \
# -- sjdbGTFfile $ANNO \
# -- sjdbOverhang 100 -- genomeSAindexNbases 12;

############################ 1.2 Mapping #####################################

# Mapping a list FASTA-files to the STARindex ${STARINDEX_OUTPUT} 

# raw data are zipped fasta files or non zipped possible  
# zipped = --readFilesCommand zcat
# and .gz 

# The sample_list file must be created beforehand (first column = name of the FASTA-file 
# without the "_1/2", second column = arbitrary identifier, for example "d1_adult_whole_1") 
# Second column can also be empty
# Example: 	G1-N2	d1_adult_gonad_wildtype_1
#			G2-N2	d1_adult_gonad_wildtype_2
#			G3-N2	d1_adult_gonad_wildtype_3
#			G4-N2	d1_adult_gonad_wildtype_4
#			G5-N2	d1_adult_gonad_wildtype_5

echo 'mapping: ' ${SAMPLE_LIST} 
cat ${SAMPLE_LIST} | cut -f2 | sed 's/\r$//' | while read sample; do 
	echo $sample `date +%T`;
	sra=`grep $sample ${SAMPLE_LIST} | cut -f1`;

	#check if paired end or single end read
	FILE=${RAW_SEQ}${sra}_1.fastq.gz
	if test -f "$FILE"; then #if paired end
		echo "$FILE exists."
		#location of the raw data files must be defined
		forw=`ls ${RAW_SEQ}${sra}_1.fastq.gz`; #.gz
		rev=`ls ${RAW_SEQ}${sra}_2.fastq.gz`; #.gz
		echo $forw;
		echo $rev;
			  
		#sub-directories for the resultinig files of each sample are created
		mkdir -p ${OUTPUT}/star/$sample/;
		# STAR itself
		STAR \
		--runThreadN 20 \
		--genomeDir $STARINDEX_OUTPUT \
		--readFilesIn $forw $rev \
		--outFileNamePrefix ${OUTPUT}/star/$sample/ \
		--outSAMtype BAM Unsorted \
		--quantMode GeneCounts \
		--readFilesCommand zcat \
		--outSAMstrandField intronMotif
	else #if signle end
		forw=`ls ${RAW_SEQ}${sra}.fastq.gz`; #.gz
		echo $forw;
		#sub-directories for the resultinig files of each sample are created
		mkdir -p ${OUTPUT}/star/$sample/;
		# STAR itself
		STAR \
		--runThreadN 20 \
		--genomeDir $STARINDEX_OUTPUT \
		--readFilesIn $forw\
		--outFileNamePrefix ${OUTPUT}/star/$sample/ \
		--outSAMtype BAM Unsorted \
		--quantMode GeneCounts \
		--readFilesCommand zcat \
		--outSAMstrandField intronMotif
		echo "$FILE does not exist."
	fi

done;

# ############################ 1.3 Sort bed files ######################
echo 'sort';
cat ${SAMPLE_LIST} | cut -f2 | sed 's/\r$//' | while read sample; do 
	echo $sample;
	samtools sort ${OUTPUT}/star/$sample/Aligned.out.bam \
		-o ${OUTPUT}/star/$sample/Aligned.out.bam_sorted.bam;
	done; 

############################ 2.1 Stingtie ######################################
# you might want to skip mitochondiral genes, this can be done with -x chrM
echo 'stringtie';
cat ${SAMPLE_LIST} | cut -f2 | sed 's/\r$//' | while read sample; do 
	echo $sample `date +%T`;
	mkdir -p ${OUTPUT}/stringtie/;
	stringtie -c 0.1 -f 0.0 -m 50 -a 1 -p 20 ${OUTPUT}/star/$sample/Aligned.out.bam_sorted.bam \
	  -o ${OUTPUT}/stringtie/$sample.gtf;
	# -c -> minimum read coverage 
	# -f -> disable trimming at the end of the transcripts ?! -t 
	# -m -> minimum length allowed for the predicted transcripts (normal=200)
	# -a -> filtered out junctions that dont have spliced reads align with at least 1
	# -p -> number of processing threads 
done; 

# Splice-O-Mat applied to adhesion GPCRs transcript variants

This repository contains all informations and scripts for the paper: <br>
[![DOI](https://zenodo.org/badge/696187476.svg)](https://zenodo.org/badge/latestdoi/696187476) <br>
**The repertoire and structure of adhesion GPCR transcripts assembled from deep-sequenced human samples.** <br>
**Christina Katharina Kuhn, Udo Stenzel, Sandra Berndt, Ines Liebscher, Torsten Schöneberg, Susanne Horn** <br>
Under Review at Nuclar Acid Research (NAC) <br> 

This includes preprocessing of the datasets (under analysis), the scripts for the webtool (under scripts), and additional analysis for the manuscript (under analysis).  <br>

## Table of Contents
- [Objective](#objective)
- [Data](#data)
  - [GSE182321_OpiodUseDisorder_braintissues](#gse182321_opiodusedisorder_braintissues)
  - [GSE173955_Alzeihmer_braintissue](#gse173955_alzeihmer_braintissue)
  - [GSE174478_Non_alcoholic_fatty_liver](#gse174478_non_alcoholic_fatty_liver)
  - [GSE217427_kidney](#gse217427_kidney)
  - [GSE165303_SRP302848_heart](#gse165303_srp302848_heart)
  - [SRP225193_many_tissues](#srp225193_many_tissues)
  - [melanoma](#melanoma)
- [Analysis](#analysis)
- [Results](#results)
- [Webtool Scripts](#webtool-scripts)
- [Additional](#additional)
- [Contact](#contact)
- [Todos](#todos)

## Objective

The main objective of this project was to analyze tissue-specific splicing of adhesion GPCRs (aGPCRs).
For this a webtool was created based on a database including over 900 samples and 48 different tissue types.
For investigation of tissue-specific splice pattern a gene of interest, please visit: https://tools.hornlab.org/Splice-O-Mat/

## Data

### GSE173955_Alzeihmer_braintissue
- GSE173955
- BioProject: PRJNA727602
- SRP ID: SRP318632
- Samples: 40
- Description: Analysis of postmortem human hippocampus brains from 8 subjects with Alzheimer's disease (AD) and 10 non-AD subjects using Illumina TruSeq stranded mRNA LT Sample Prep kit. Sequencing performed on HiSeq1500. One AD and one non-AD sample were applied independently twice to increase read coverage.
- Sequencing Depth: 9-25M reads per sample
- Participants: 8 AD subjects and 10 non-AD subjects
- PMID: 23595620

### GSE182321_OpiodUseDisorder_braintissues
- GSE182321
- Description: RNA sequencing analysis of postmortem human Brodmann Area 9 in the University of Texas Health Science Center at Houston Brain Collection for individuals with Opioid Use Disorder; 27 opioid users and 14 nonpsychiatric controls, as determined by postmortem consensus diagnosis by two trained psychiatrists.
- BioProject: PRJNA755746
- SRP ID: SRP332964
- Samples: 41
- Description: RNA sequencing analysis of postmortem human Brodmann Area 9 in the University of Texas Health Science Center at Houston Brain Collection for individuals with Opioid Use Disorder.
- Sequencing Depth: 27-34M  reads per sample
- PMID: 34385598

### GSE101521 Brain Dataset

- GSE101521
- Decription: on-psychiatric controls (CON, N=29), DSM-IV major depressive disorder suicides (MDD-S, N=21) and MDD non-suicides (MDD, N=9) in the dorsal lateral prefrontal cortex (Brodmann Area 9)
- Samples: 59
- Sequencing Depth: 7-61 million reads per sample
- Participants: 8 AD subjects and 10 non-AD subjects
- PMID: 27528462

### GSE174478_Non_alcoholic_fatty_liver
- GSE174478
- Description: a fatty liver diagnosed ultrasonically by an increase in hepatorenal contrast, a history of alcohol consumption of less than 30 g/d for men and less than 20 g/d for women, seronegativity for hepatitis B virus surface antigen and hepatitis C virus antibody, and the absence of autoimmune hepatitis, primary biliary cholangitis, primary sclerosing cholangitis, Budd-Chiari syndrome, Wilson disease, and drug-induced liver injury
- SRP ID: SRP319881
- BioProject: PRJNA730024
- Samples: 94
- Sequencing Depth: 31-44M reads per sample
- PMID: 35380992

### GSE217427_kidney
- GSE217427
- Description: medulla and coretex, with human kidney damage (KD) (n=22) and without KD (22),
- BioProject: 
- Samples: 44
- Sequencing Depth:37-51M reads per sample
- Not published yet
 
 ### GSE165303_SRP302848_heart
 - GSE165303
 - Description: with dilated cardiomyopathy, 50 non-failing, 2 Transfected with control adenovirus
and 2 transfected with HAND1 overexpressing adenovirus, only paired end was selected
- BioProject: SRP302848
- Samples: 101
- Sequencing Depth:57-80M reads per sample
- Not published yet

### SRP225193_many_tissues
- GSE138734
 - Description: 300 human samples, including 45 tissues, 162 cell types, and 93 cell lines, some paired some single end, total RNA (296 samples), only full RNA-seq (paired end) of 45 tissues was selected (cell line and cell types excluded)
- BioProject: SRP225193 
- Samples: 457
- Sequencing Depth:76-111M reads per sample
- PMID: 34140680

### melanoma 
- Description: RNA-seq of metastatic melanoma patients treated with anti-PD-1 alone or combined anti-PD-1 and anti-CTLA-4 immunotherapy
- BioProject: PRJEB23709 at ENA
- Samples: 91
- Sequencing Depth: around 50M reads per sample
- PMID: 30753825

## Analysis

The following steps were followed to perform the mapping/assembly and creation of the database:<br>
All files are within directory `analysis/`

1. SRR data retrieval:
   - Execute `../get_SRR_data.sh` to retrieve the SRR data.

2. STAR + StringTie:
   - Execute `../splice-variant-analysis.sh` to perform STAR mapping (sorting and indexing) and StringTie Assembly with hg38. Run this script in the `data/cohort/` directory, using "cohort" as the output name in the `analysis/cohort/` directory.

3. StringTie merged mode:
   - Execute `../stringtie_expression_estimation.sh` to generate a merged GTF file. This file combines multiple GTFs based on the `directories.txt` file. If a combo file already exists, include it in `new_mergefile.txt`. Requantification is performed using the combo GTF file.

4. Building StringTie.db:
   - Clone the repository: `git clone https://chrissi_kath@bitbucket.org/ustenzel/stringtiedb.git`
   - Build: `cabal build`
   - Update (if necessary): `cabal update`
   - Run StringTie-db:
     - `cabal run stringtie-db -- -d stringtie_2.db -m ../analysis/combo_new.gtf ../analysis/ballgown_version2//.gff -C ../analysis/pheno_data.csv`
     - `cabal run stringtie-db -- -d stringtie_2.db -C ../analysis/pheno_data.csv` (to update the samples table)
   - The resulting `stringtie_2.db` is the new database.

Other necessary files:
- ../analysis/hg38.fna

Additional analysis in the manuscript: 
- ../analysis/compare_CLSR1_with_genocode_and_diagostics/: Comparison with diagnostic exome sequencing of CELSR1/ADGRC1
- ../analysis/tissue_specific_splicing/: Spearman and violin plot of tissue-specific splicing

## Results 
Inside the `analysis/` directory, you will find the following directories for each dataset:

- ../cohort/star: Contains mappings generated by STAR.
- ../cohort/stringtie: Contains transcripts obtained from StringTie.
- ../cohort/ballgown: Contains outputs for database, including requantification with merged GTF.
- ..bzw ballgown_redone: Additional output for database, requantification with merged GTF.

## Webtool Scripts
The following scripts are available for the webtool: <br>
All files are within directory `scripts/`

- `../webtool.py`: DASH webtool named "Splice-o-mat".
- `../bootstrap.min`: CSS file for styling.
- `../assets/`: Images related to the webtool.
- `/var/tmp/process_id/`: Temporary storage directory for intermediate data, including SVG, TXT, and FASTA files.

### Dependencies
#### Python and Packages
- python version: Python 3.10.6
- used pyhton packages under: `../scripts/requirements.txt`

#### my_interproscan 
The is a need to include a local version of interproscan <br>
to search  domains in longest ORF of transcript <br>
version used: interproscan-5.60-92.0 <br>

- `../my_interproscan` includes a local Pfam database so search for domains in proteins

- with the script: `..analysis/insertDomainsinDb.py` domains can be added into the stringtie_2.db as additional domain table to save computation time of the webtool
 (has to be generated with the structure: domains (transcript, domain, start, end))

## Contact
Do you have any questions, suggestions about the webtool or the analysis please write me an [email](mailto:christina.kuhn@medizin.uni-leipzig.de)

---
## ToDos:
- [ ] 'gene_name' benutzt anstatt 'gene_id' benutzen
- [X] mitochondriale Gene löschen
- [ ] Reinzoomen
- [ ] alle Domänen anzeigen
- [ ] responsive
- [ ] an eine SNP Datenbank anknüpfen
- [ ] Svg flippen
- [ ] Eigene Domänen eintragen
- [ ]  Mutation labeln können (rs Code)
- [ ] Eigene Daten einspeißen (pipeline, datenbank etc)
  - müssen deep sequenced sein, paired-end, datenbank local speichern?
- [ ] gene_id entfernen, stattdessen gene_name

[Back to the top](#top)

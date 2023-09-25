import sqlite3
import pandas as pd
import pybedtools
from Bio.Seq import Seq
import regex as re
import subprocess
import os

def find_longest_ORF(cDNA, strandness, transcript):
    '''
    Get start and stop postition of longest ORF in cDNA
    Checks if gene comes from + or - strand
    :param DNA, strandness
    :return start_pos, stop_position <- on cDNA on original genomic direction, orf_sequence, protein_sequence
    '''
    # print("--Find ORFs of ",transcript,"--")
    originalcDNA = cDNA #make a copy of DNA
    start = 0
    ORFs = [] #list with all ORFs and start and end position

    # if strand is minus, take the reverse complement if the cDNA 
    # in order to search for the longest ORF
    if strandness == "-":
        cDNA = str(Seq(cDNA).reverse_complement())
        
    # print("cDNA:", cDNA)
    # print("originalcDNA:", originalcDNA)
    if 'ATG' in cDNA:
        for startMatch in re.finditer('ATG',cDNA): #check if substring starts with ATG
            remaining = cDNA[startMatch.start():]
            for stopMatch in re.finditer('TAA|TGA|TAG', remaining): #check for substring ending with stopcodons
                substring = remaining[:stopMatch.end()]
                if len(substring) % 3 == 0: #check if triplets
                    start = startMatch.start() #get start position in cDNA or if "-" on reverse complement cDNA
                    end = start+(len(substring)) #get end position start+len()-1
                    ORFs.append([substring,start,end]) #add ORF and start and end position to list
                    break #break when first stop codon is seen
    
    #sort list of lists [[ORF,start,stop],...] by len of first element
    ORFs = sorted(ORFs, key=lambda x: len(x[0]), reverse=True)
    # print(ORFs[0][0][0:10], len(ORFs[0][0]), ORFs[1][0][0:10], len(ORFs[1][0]), ORFs[2][0][0:10], len(ORFs[2][0]))
    
    if len(ORFs)==0:
        return (0, 0, "No ORF found", None)
    
    longest = ORFs[0] # save first item in list as longest ORF
    # match the start and stop position of the ORF on the original cDNA
    if strandness == "-": #if strand is minus, get original positions on + strand
        start_pos = (len(originalcDNA)-1)-longest[2]
        end_pos = (len(originalcDNA)-1)-longest[1]
    else:
        start_pos = longest[1]
        end_pos = longest[2]

    #write proteins to file
    with open("proteins.fasta","a") as f: 
        if len(ORFs)>=1:
            cDNA_Seq = Seq(ORFs[0][0])
            protein = cDNA_Seq.translate(to_stop=True, cds=True)
            f.write(">"+transcript+"\n"+str(protein)+"\n")
    f.close()
    
    return (start_pos, end_pos, longest[0], str(protein))    

con = sqlite3.connect("../analysis/stringtie_2.db")
cur = con.cursor()
# cur.execute("CREATE TABLE domains (id integer primary key, transcript references transcripts, domain text, start int, end int)")
# con.commit()

cur.execute("SELECT * FROM transcripts")
transcripts = cur.fetchall()
print(len(transcripts))

#if proteins.fasta exists, delete it
if os.path.exists("proteins.fasta"):
    os.remove("proteins.fasta")

for transcript in transcripts:
    statement = "SELECT e.chrom, e.start, e.end, e.strand, t.transcript_id, t.id from transcripts as t, exons as e WHERE e.transcript==t.id AND t.id==?"
    cur.execute(statement, (transcript[0],))
    #get cDNA sequence of transcript
    cDNA_transcript=""
    #write bed file with exons
    df = pd.DataFrame(cur.fetchall())
    df.columns=["chromosom", "start", "end", "strand", "transcript_id","id"]
    #take only first 3 columns
    df_bed = df.iloc[:,0:3]
    #safe tmp bed file
    df_bed.to_csv("gene.bed",sep="\t", index=False, header=False)
    gene_bed = pybedtools.BedTool("gene.bed")
    fasta = pybedtools.BedTool("hg38.fna")
    a = gene_bed.sequence(fi=fasta)
    with open(a.seqfn) as f:
        for line in f:
            if not(line.strip().startswith(">")):
                cDNA_transcript += line.strip()
    f.close()
    cDNA_transcript += cDNA_transcript
    #write cDNA sequence to file with ">transcript_id\n"
    with open("cDNA.fasta","a") as f:
        f.write(">"+df["transcript_id"][0]+"\n"+str(cDNA_transcript)+"\n")
    
    print("Transcript:", df["transcript_id"][0], " with id: ", df["id"][0])
    #get longest ORF of cDNA sequence
    start, stop, longest, protein = find_longest_ORF(cDNA_transcript,df["strand"][0],df["transcript_id"][0])
    #get protein sequence of longest ORF
    #search for domains in protein sequence
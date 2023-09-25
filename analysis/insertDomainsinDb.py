import sqlite3
import pandas as pd
import pybedtools
from Bio.Seq import Seq
import regex as re
import subprocess
import os
con = sqlite3.connect("../analysis/stringtie_2.db")
cur = con.cursor()

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
                    end = start+(len(substring)-1) #get end position start+len()-1
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

    #only the longest ORF to file
    # with open("longest_ORF.fasta","a") as f:
    #     if len(ORFs)>=1:
    #         f.write(">"+transcript+"\n"+str(ORFs[0][0])+"\n")
    # f.close()

    # #write proteins to file
    # with open("proteins.fasta","a") as f: 
    #     if len(ORFs)>=1:
    #         cDNA_Seq = Seq(ORFs[0][0])
    #         protein = cDNA_Seq.translate(to_stop=True, cds=True)
    #         f.write(">"+transcript+"\n"+str(protein)+"\n")
    # f.close()

    if len(ORFs)>=1:
        cDNA_Seq = Seq(ORFs[0][0])
        protein = cDNA_Seq.translate(to_stop=True, cds=True)
    
    return (start_pos, end_pos, longest[0], str(protein))

try:
    df_domains = pd.read_csv("all_domains.tsv", sep="\t", header=None) 
    df_domains.drop([1,2,3,4,8,9,10,11,12], axis=1, inplace=True)
    df_domains.columns=["transcript", "domain", "start", "end"]
except:
    print("no domain found")


cur.execute("SELECT * FROM transcripts")
transcripts = cur.fetchall()
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
    #safe tmp bed fileprint
    
    df_bed.to_csv("gene.bed",sep="\t", index=False, header=False)
    gene_bed = pybedtools.BedTool("gene.bed")
    fasta = pybedtools.BedTool("hg38.fna")
    a = gene_bed.sequence(fi=fasta)
    with open(a.seqfn) as f:
        for line in f:
            if not(line.strip().startswith(">")):
                cDNA_transcript += line.strip()
    f.close()
    #get longest ORF of cDNA sequence
    start_ORF_on_plus, stop_ORF_on_plus, ORF, protein = find_longest_ORF(cDNA_transcript,df["strand"][0],df["transcript_id"][0])

    if ORF == "No ORF found":
        continue

    #get all rows in df_domain where transcript_id == transcript
    df_domains_for_transcript = df_domains[df_domains["transcript"]==df["transcript_id"][0]]

    #if df_domains_for_transcript is empty, continue
    if df_domains_for_transcript.empty:
        continue
    
    print("Transcript:", df["transcript_id"][0], " with id: ", df["id"][0])
    #for each domain in database 
    for index, row in df_domains_for_transcript.iterrows():
        domain = row["domain"]
        #clean domain string from unwanted characters
        domain = re.sub(r'\(.*?\)', '', domain)
        start_domain_on_ORF = row["start"]*3
        end_domain_on_ORF = row["end"]*3

        if df["strand"][0] == "+":
            position_of_ORF_in_cDNA = cDNA_transcript.find(ORF)
            #jetzt muss ich ja eigentlich diese Position nehmen + der Start der Domäne in dem ORF;
            #und das gleiche für das ende der Domäen
            start_domain_on_cDNA = position_of_ORF_in_cDNA + start_domain_on_ORF
            end_domain_on_cDNA = position_of_ORF_in_cDNA + end_domain_on_ORF
        else:
            position_of_ORF_in_cDNA = start_ORF_on_plus
            # Jetzt muss ich diese Domänen die ja im ORF (length 3966) liegen wieder auf die DNA in plus strang richtung "mappen".
            # Hier muss die Position die ja vom ORF auf dem minus strang kamen erstmal auf das ORF in plus strang richtung "mappen".
            # Heißt: domäne_start_on_plus_ORF = len(ORF)-end_domäne(423)
            #  domäne_end_on_plus_ORF = len(ORF)-start_domäne(255)
            end_domain_on_plus_ORF = len(ORF)-start_domain_on_ORF
            start_domain_on_plus_ORF = len(ORF)-end_domain_on_ORF
            start_domain_on_cDNA = position_of_ORF_in_cDNA + start_domain_on_plus_ORF
            end_domain_on_cDNA = position_of_ORF_in_cDNA + end_domain_on_plus_ORF

        # domains.append((df["transcript_id"][0], domain, start_domain_on_cDNA, end_domain_on_cDNA))
        # save domains in database
        statement = "INSERT INTO domains (transcript, domain, start, end) VALUES (?,?,?,?)"
        try:
            cur.execute(statement, (int(df["id"][0]), domain, start_domain_on_cDNA, end_domain_on_cDNA))
            con.commit()
        except:
            print("Domain already exists")
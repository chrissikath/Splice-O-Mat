# Author: Christina Kuhn
# Date: 2022-03-10
# Description: webtool for splice-o-mat, which is a tool to visualize splice-variants
#               from RNA-seq data and compare them in different samples and tissues
#              the webtool is hosted at https://tools.hornlab.org/Splice-O-Mat/
# Use Case: GPCR analysis as in the paper from 10.1038/s41598-019-46265-x.

################## 0.0 Import Modules ############################
import datetime
import flask
import pandas as pd
from dash import Dash, dash_table, dcc, html, Input, Output, State, exceptions, no_update, ctx
import sqlite3
import dash_bootstrap_components as dbc
import svgwrite
import base64
from io import BytesIO
import twobitreader
import seaborn as sns
import plotly.express as px
import numpy as np
import matplotlib.pyplot as plt
import os
from Bio.Seq import Seq
import regex as re
import subprocess
import math
from scipy import stats
import json
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from config import *


# create a database connection to the SQLite database

# sqlite3 reports thread safety level "multithreaded", which means
# check_same_thread=False should be safe.  isolation_level=None means not
# implicit transactions are created.  We're only reading, so we don't need
# them and they can only cause spurious locking issues.
#
# Using a global database object means we can extend and configure it
# once, e.g. by defining functions or loading extensions.  Our extension
# implements access to compact storage of bulk data, and the stdev
# aggregate function.

con = sqlite3.connect(database_file,check_same_thread=False,isolation_level=None)
con.enable_load_extension(True)
con.load_extension("./mod_array")
con.enable_load_extension(False)                # to prevent arbitrary SQL from loading more extensions
r = con.execute("select attachArray(?)", (array_file,))
print(r.fetchone()[0])                          # prints "Attached..." if successful



def fill_dropdowns_in_dash_application():
    """Fill all dropdowns in dash application at the beginning
    Fill main dropdown for gene selection (by gene_name) and
    dropdowns for sample and tissue selection

    Returns:
        list: gene_names, tissues_and_samples, tissues
    """

    print("----load genes names----")
    gene_names = set()
    # pipe with gen localization in dropdown
    # res = con.execute("SELECT t.gene_name, e.chrom, MIN(e.start) AS smallest_start, MAX(e.end) AS biggest_end FROM transcripts t JOIN exons e ON t.id = e.transcript GROUP BY t.gene_name")
    # pipe without gen localization in dropdown
    res = con.execute("SELECT DISTINCT gene_name FROM transcripts WHERE gene_name IS NOT NULL ORDER BY gene_name")
    # gene_names = [item[0]+" "+item[1]+":"+str(item[2])+"-"+str(item[3]) for item in res.fetchall()]
    gene_names = [item[0] for item in res.fetchall()]
    gene_names = [{"label":i, "value":i} for i in gene_names]

    print("----load samples----")
    samples = set()
    res = con.execute("SELECT name FROM samples")
    samples = list(set([item[0] for item in res.fetchall()]))

    print("----load tissues----")
    tissues  = set()
    res = con.execute("SELECT DISTINCT tissue FROM samples")
    tissues = list(set([item[0] for item in res.fetchall()]))

    # tissues sort by given order in list was given by Torsten Schoenberg
    tissues = ["brain","brain opioid","brain alzheimer","brain fetal","brain frontal cortex",\
                "brain parietal cortex","brain occipital cortex","brain striatum","brain stem",\
                "brain cerebellum","esophagus","stomach","small intestine",\
                "duodenum","jejunum colon","ileum colon",\
                "colon","proximal colon","distal colon","liver","pancreas","kidney","bladder",\
                "breast","ovary","fallopian tube","uterus","cervix",\
                "placenta","prostate","testicle","heart","pericardium","right atrium heart",\
                "left atruim","right ventricle heart","left ventricle heart",\
                "vena cava","trachea","lung","thymus","spleen","lymph node","thyroid",\
                "adrenal","adipose","skeletal muscle","melanoma"]

    tissues_and_samples = tissues + samples
    tissues_and_samples = [{"label":i, "value":i} for i in tissues_and_samples]
    tissues = [{"label":i, "value":i} for i in tissues]
    return gene_names, tissues_and_samples, tissues

def b64_image(image):
    """encodes image for embedding in website

    Args:
        image (string): a png image

    Returns:
        bytes: rendered image encoded as png and base64
    """    """
    """
    return 'data:image/png;base64,' + base64.b64encode(image).decode('utf-8')

def b64_svg(drawing):
    """encodes svg for embedding in website

    Args:
        drawing (Drawing): an svgwrite.Drawing

    Returns:
        bytes: encoded base 64 svg
    """    """
    """
    return 'data:image/svg+xml;base64,' + base64.b64encode(drawing.tostring().encode('utf-8')).decode('utf-8')

def transform_to_columns_and_data(result_df, exon_number=True):
    """Adds a href link to the NCBI entry for each transcript_id which is not NSTRG

    Args:
        result_df (Dataframe): complete dataframe

    Returns:
        result_df (Dataframe): dataframe with hrefs
    """
    result_df = pd.DataFrame(result_df.apply(lambda row: {col: f"[{row[col]}](https://www.ncbi.nlm.nih.gov/nuccore/{row['transcript_id']})"
                                                      if (col == 'transcript_id' and not row[col].startswith('NSTRG'))
                                                      else row[col] for col in result_df.columns}, axis=1).tolist())

    if exon_number:
        number_inner_exons = calculate_inner_exons(result_df["gene_name"][0])
        result_df.loc[-1] = "" * len(result_df.columns)  # adding a row
        result_df.loc[-1, "# of exons"] = number_inner_exons

    columns = [
                {"name": i, "id": i, "presentation": "markdown"} if i == "transcript_id" else {"name": i, "id": i}
                for i in result_df.columns
            ]
    data = result_df.to_dict('records')

    return columns, data

def get_start_and_end_genomic_region(transcripts):
    """Get start coordinate and end coordinate of genomic region were
    first exon of first transcript and last exon of last transcript are located

    Args:
        transcripts (list): liste of transcript ids

    Returns:
        start_genomic_region (int): start and end coordinate of genomic region
        end_genomic_region (int): start and end coordinate of genomic region
        y (int): y coordinate for svg how much space is needed for transcript names
    """

    #initialize start and end of genomic region
    start_genomic_region = float("inf")
    end_genomic_region = 0
    y = 0
    #go through transcript and get start first exon = start_genomic_region and last exon = end_genomic_region
    for transcript in transcripts:
        y += 20
        statement = "SELECT e.start, e.end from transcripts as t, exons as e WHERE e.transcript==t.id and t.transcript_id==?"
        res = con.execute(statement, (transcript,))
        exons = res.fetchall()
        if (exons[0][0]) < start_genomic_region:
            start_genomic_region = exons[0][0]
        if (exons[-1][1]) > end_genomic_region:
            end_genomic_region = exons[-1][1]
    return start_genomic_region, end_genomic_region, y

def get_exons_from_transcript_id(transcript_id):
    """Get all exons from a transcript id

    Args:
        transcript_id (string): transcript_id

    Returns:
        exons (list): list of tuples with start, end, chrom, strand
    """    """
    """

    #database statement
    statement = "SELECT e.start, e.end, e.chrom, e.strand from transcripts as t, exons as e WHERE e.transcript==t.id and t.transcript_id==?"
    res = con.execute(statement, (transcript_id,))
    exons = res.fetchall()
    return exons

def split_fasta_string_to_dict(cDNAs):
    """Split a list of cDNAs into a dictionary with the transcript id as key and the sequence as value
    Since it was a fasta file string, the 1,3,5th.. element is the key and the 2,4,6th.. element is the value
    splitted by '\n'
    Args:
        cDNAs (string): content of a fasta file

    Returns:
        cDNAs_dict (dict): dictionary with transcript id as key and sequence as value
    """    """"""
    cDNAs = cDNAs.split("\n")
    keys = cDNAs[0::2]
    keys_new = [key[1:] for key in keys]
    values = cDNAs[1::2]
    cDNAs_dict = dict(zip(keys_new, values))
    return cDNAs_dict

#get the key:value of the len(value)==max
def get_max_flow(dict_input):
    """Get the key of the longest value in a dictionary
    meands the longest transcript

    Args:
        dict_input (dict): dictionary with transcript id as key and sequence as value

    Returns:
        value (string): longest protein sequence
        key (string): longest transcript id
    """    """
    """
    key_of_longest=max(dict_input, key=lambda k: len(dict_input[k]))
    return dict_input[key_of_longest],key_of_longest

def find_domains_from_database(longest_transcript):
    """Find domains in a protein sequence from a database
    and return the domains and their start and end on the DNA on plus strand
    #! if the gene lays on the minus strand, reverese complement is used, however,
    #! the return values start and end on the DNA is still on the plus strand for plotting on the svg
    # wirkt einfach, war es aber nicht
    Args:
        longest_transcript(string): a protein, using the standard single-letter code

    Returns:
        domains (list): list of domains (transcript_id, domain, start_on_cDNA, end_on_cDNA)
    """
    print("--Find domains from database--")

    #get domains from database
    statement = "SELECT d.domain, d.start, d.end FROM domains as d, transcripts as t WHERE t.id=d.transcript AND t.transcript_id=?"
    res = con.execute(statement,(longest_transcript,))
    df = pd.DataFrame(res.fetchall())
    if not df.empty:
        columns = ["domain","start","end"]
        df.columns = columns

    #for each entry in df append longest_transcript, domain, start, end to domains
    domains = []
    for _, row in df.iterrows():
        domains.append((longest_transcript, row["domain"], row["start"], row["end"]))
    return domains

def generate_heatmap(result_df, relatives=True):
    """
    Generate a heatmap of transcripts and their TPMs.

    Args:
        result_df (pd.DataFrame): Dataframe of transcripts and their TPMs.
        relatives (bool, optional): If True, the heatmap will be generated with relative values.If False heatmap \
            produces absolutes TPMS. Defaults to True.

    Returns:
        a png file in a buffer
    """

    print("--Generate heatmap--")
    try:
        transcripts = result_df["transcript_id"]
        if relatives==True:
            tpms_values = result_df.filter(regex='%')
            columns_tmp = tpms_values.columns
            columns_tmp = [col.replace('TPM(%)', '') for col in columns_tmp]
            tpms_values.columns = columns_tmp
            return plot_heatmap(tpms_values, transcripts, True)
        else:
            tpms_values = result_df.filter(regex='mean')
            columns_tmp = tpms_values.columns
            columns_tmp = [col.replace('TPM(mean)', '') for col in columns_tmp]
            tpms_values.columns = columns_tmp
            return plot_heatmap(tpms_values, transcripts, False)
    except Exception as e:
        print(f"An error occurred: {str(e)}")

def plot_heatmap(tpms, transcripts, relative=True):
    """
    Plot the heatmap of TPM values of transcripts.

    Args:
        tpms_relative_values (pd.DataFrame): Dataframe containing relative TPMs.
        transcripts (pd.Series): Series of transcript names.
    """

    x_labels = tpms.columns
    y_labels = transcripts

    if relative:
        cbar_legend = 'TPM (%, per tissue)'
        title_figure = "Relative expression of transcripts (% of total in each tissue)"
    else:
        cbar_legend = 'TPM (mean)'
        title_figure = "Absolute expression of transcripts (TPM)"

    fig = px.imshow(tpms.values,
                labels=dict(x="Tissues", y="Transcripts", color=cbar_legend),
                x=x_labels,
                y=y_labels,
                color_continuous_scale='Inferno',
                aspect="auto",
                width=num_entries_based_width(len(x_labels)),
                height=num_fields_based_height(len(y_labels)),
               )


    fig.update_layout(title_text=title_figure, title_x=0.5, title_font_size=16)
    fig.update_layout(xaxis_nticks=len(x_labels),yaxis_nticks=len(y_labels))
    fig.update_xaxes(tickangle=-45)

    fig.update_yaxes(tickmode='array', tickvals=list(range(len(y_labels))), ticktext=y_labels)


    try:
        return fig
    except Exception as e:
        print(f"An error occurred while saving the heatmap: {str(e)}")

def num_fields_based_height(num_fields: int) -> int:
    padding = 0 # arbitrary value depending on legends
    row_size = 40 # arbitrary value
    min = 400
    return max(min, padding + row_size * num_fields)

def num_entries_based_width(num_entries: int) -> int:
    padding = 0 # arbitrary value depending on legends
    entry_size = 28 # arbitrary value
    return padding + entry_size * num_entries

def generate_svg(transcripts, position_mut=None):
    """Generates a svg for a transcript variants of the list transcripts.
    Also plots the mutation position if position_mut is given.
    # TODO - x and y axis with more meanigful names

    Args:
        transcripts (list): list of transcript ids like '[NSTRG.8029.1, NSTRG.8029.2]'
        position_mut (string, optional): Mutations or list of mutations like '12344566,123456' Defaults to None.

    Returns:
        int, int, string, string, drawing: start_genomic_region_of_transcript, end_genomic_region_of_transcript, chrom, strand, drawing
    """

    print("--Generate svg--")

    start_genomic_region_of_transcript, end_genomic_region_of_transcript,y = get_start_and_end_genomic_region(transcripts)
    absolut = (end_genomic_region_of_transcript-start_genomic_region_of_transcript) #100% = genomic region in absolut number

    # [ ] TODO - x and y
    x_start = 0
    y_start = 20

    size = 2000 #size of whole svg
    d = svgwrite.Drawing(viewBox="0 0 "+ str(size+300) +" "+ str(y+450) +"")
    #explain colors in svg
    d.add(svgwrite.text.Text("red = displays the longest ORF for a potentially new transcript (NSTRG_), might not be the 'true' one.",\
                              insert=(x_start, y_start), style="font-size:24"))
    d.add(svgwrite.text.Text("orange = displays the longest ORF for the transcript, might not be the 'true' one.", \
                             insert=(x_start, y_start+20), style="font-size:24"))
    d.add(svgwrite.text.Text("yellow = domains presented are predicted for the longest protein, \
                             to get all domains for each transcript press button 'Get domains'", \
                                insert=(x_start, y_start+40), style="font-size:24"))

    #for all transcript_ids get gene_name and smallest start and biggest end position
    # creat dictionary with strucutre gene_name:(start,end)
    gene_names = {}
    for transcript_id in transcripts:
        #connect to database
        #get gene_name and smallest start and biggest end position of transcript
        statement = "SELECT t.gene_name, min(e.start), max(e.end) FROM transcripts as t, exons as e WHERE e.transcript==t.id and t.transcript_id==?"
        res = con.execute(statement, (transcript_id,))
        gene_name, start, end = res.fetchone()
        #add gene_name to dictionary and update start if smaller and end if bigger than before
        if gene_name in gene_names:
            if start < gene_names[gene_name][0]:
                gene_names[gene_name][0] = start
            if end > gene_names[gene_name][1]:
                gene_names[gene_name][1] = end
        else:
            gene_names[gene_name] = [start, end]

    formatted_entries = []
    for gene_name, values in sorted(gene_names.items()):
        formatted_entry = f"{gene_name} ({values[0]}, {values[-1]})"
        formatted_entries.append(formatted_entry)

    resulting_string = ', '.join(formatted_entries)

    d.add(svgwrite.text.Text("presented genes (start of first exon, end of latest exon): "+resulting_string, \
                                insert=(x_start, y_start+80), style="font-size:24"))

    x_start = 0
    y_start = y_start+100

    #draw genomic position
    genomic_region = svgwrite.shapes.Rect(insert=(x_start+40, y_start+20), size=(size, 5),fill='grey')
    start = svgwrite.text.Text("Start: "+str(start_genomic_region_of_transcript), insert=(x_start+0,y_start+15), style="font-size:24")
    end = svgwrite.text.Text("End: "+str(end_genomic_region_of_transcript), insert=(size,y_start+15), style="font-size:24")
    d.add(genomic_region)
    d.add(start)
    d.add(end)

    # for gene_name, values in sorted(gene_names.items()):
    #     print(gene_name, values)
    #     gene_position_start = ((absolut-(end_genomic_region_of_transcript-(values[0])))*size)//(absolut)
    #     gene_position_end = ((absolut-(end_genomic_region_of_transcript-(values[1])))*size)//(absolut)
    #     gene_name_on_genome =  svgwrite.text.Text(str(gene_name), insert=(40+((gene_position_end-gene_position_start)/2), y_start+10), style="font-size:24")
    #     gene_location_on_genome = svgwrite.shapes.Rect(insert=(40+((gene_position_start)), y_start+20), size=(2,200), fill='#bc5090')
    #     d.add(gene_name_on_genome)
    #     d.add(gene_location_on_genome)

    y_start = y_start+20

    longest_transcript=None
    max_length=0

    #draw each exon in svg
    y = y_start # y-axis of each exon in svg
    for transcript in transcripts: #for each transcript draw exon, intron, ORF
        y += 20
        exons = get_exons_from_transcript_id(transcript)
        # for each exon in list of exons, calculate position in svg
        for j in exons:
            start_exon = ((absolut-(end_genomic_region_of_transcript-j[0]))*size)//(absolut)
            end_exon = ((absolut-(end_genomic_region_of_transcript-j[1]))*size)//(absolut)
            exon = svgwrite.shapes.Rect(insert=(40+start_exon, y), size=((end_exon-start_exon), 5),fill='#58508d')
            d.add(exon)
        # and calculate first/last exon position in svg
        first_exon = ((absolut-(end_genomic_region_of_transcript-exons[0][0]))*size)//(absolut)
        last_exon = ((absolut-(end_genomic_region_of_transcript-exons[-1][1]))*size)//(absolut)
        #draw line from first exons to last exon
        line = svgwrite.shapes.Rect(insert=(40+first_exon, y+2), size=(last_exon-first_exon, 0.5), fill='#58508d')
        d.add(line)

        chrom = exons[0][2]
        strand = exons[0][3]
        start_pos, end_pos, orf_start_in_genome, orf_end_in_genome = find_longest_ORF(transcript) #find start stop in cDNA of longest ORF

        if start_pos!=None:
            if abs(end_pos-start_pos)+1>max_length:
                max_length=abs(end_pos-start_pos)+1
                longest_transcript=transcript

            # calculate start/end position of open reading frame
            start = ((absolut-(end_genomic_region_of_transcript-min(orf_start_in_genome,orf_end_in_genome)))*size)//(absolut)
            end = ((absolut-(end_genomic_region_of_transcript-max(orf_start_in_genome,orf_end_in_genome)))*size)//(absolut)
            # draw orf with no opacity
            d.add( svgwrite.shapes.Rect(insert=(40+start, y-2), size=((end-start) ,10), opacity='0.6',
                fill='#ff6361' if str(transcript).startswith("NSTRG") else '#ffa600' ))

        transcript_id = svgwrite.text.Text(str(transcript), insert=(size+50,y+5), style="font-size:24")
        d.add(transcript_id)


    #draw domains of longest protein in svg
    #find domains of each protein and return list of domains for longest protein
    domains = find_domains_from_database(longest_transcript)

    #for each domain in list of domains, calculate position in svg
    y_position = y+10
    for domain in domains:
        exons = get_exons_from_transcript_id(domain[0])
        domain_name = domain[1]
        start_domain_on_cDNA = int(domain[2])
        end_domain_on_cDNA = int(domain[3])
        #get cDNA position of each exon
        positions_of_exons_in_cDNA = get_cDNA_pos_of_each_exon(exons)
        # check in if domain pos lays in exon
        for i in positions_of_exons_in_cDNA:
            if start_domain_on_cDNA >= i[2] and start_domain_on_cDNA <= i[3]:
                domain_start_in_genome = i[0] + start_domain_on_cDNA-i[2]
            if end_domain_on_cDNA >= i[2] and end_domain_on_cDNA <= i[3]:
                domain_end_in_genome = i[0] + end_domain_on_cDNA-i[2]

        # calculate start/end position of domain
        start = ((absolut-(end_genomic_region_of_transcript-domain_start_in_genome))*size)//(absolut)
        end = ((absolut-(end_genomic_region_of_transcript-domain_end_in_genome))*size)//(absolut)
        # draw domain with no opacity

        #insert = (x,y) size = (width, height)
        domain = svgwrite.shapes.Rect(insert=(40+start, y_start), size=((end-start),y+20-y_start), \
                                      fill='yellow' ,opacity='0.2', stroke='black', stroke_width=0.3)
        d.add(domain)
        domain_name = svgwrite.text.Text(str(domain_name), insert=(start, y_position+30), style="font-size:24")
        y_position += 15
        d.add(domain_name)

    #draw mutation position in svg
    if position_mut!=None and position_mut!='':
        positions_mut = position_mut.split(",")
        for pos in positions_mut:
            mutation_position = ((absolut-(end_genomic_region_of_transcript-int(int(pos))))*size)//(absolut)
            mutation = svgwrite.shapes.Rect(insert=(40+mutation_position, y_start), size=(1,y+20-y_start), fill='#bc5090')
            d.add(mutation)

    print("--Generate svg done--")
    return start_genomic_region_of_transcript, end_genomic_region_of_transcript, chrom, strand, d

def get_cDNA_pos_of_each_exon(exons):
    """Generates a tuple for exon in exons list
    that contains the start and end position of the exon on genome and the start and end position of the exon on cDNA

    Args:
        exons (list): list of tuples with (start_exon, end_exon)

    Returns:
        list: (start_exon, end_exon, start_exon_on_cDNA, end_exon_on_cDNA)
    """
    cDNA_pos = []
    start_genom = 0
    for exon in exons:
        start_exon = exon[0] #20000
        end_exon = exon[1] #21000
        if start_genom == 0:
            cDNA_pos.append([start_exon, end_exon, 0, end_exon-start_exon-1])# zero indexin to fit to ORF position
            start_genom = start_genom+(end_exon-start_exon-1)
        else:
            cDNA_pos.append([start_exon, end_exon, start_genom+1, start_genom+1+(end_exon-start_exon-1)])
            start_genom = start_genom+1+(end_exon-start_exon-1)
    return cDNA_pos


def find_longest_ORF(transcript):
    """Find longest ORF of transcript_id "transcript" from database.
       Currently, the database layout works correctly only if all exons appear
       in syntheny.  start_pos and end_pos are the distances of the first and
       last base that are part of the longest ORF from the first base of the
       transcripts, counted on the forward strand of the exons; regardless of
       which strand in transcribed.  start_on_genome and end_on_genome are
       positions on the chromosome, and start_on_genome is smaller than
       end_on_genome if the reverse strand is transcribed.

    Args:
        transcript (string): transcript_id

    Returns:
        (int, int, int, int): start_pos on plus, end_pos on plus, start_on_genome, end_on_genome
    """    """"""

    s, e, gs, ge = con.execute( """select sum( min(e.end,t.cds_start,t.cds_end) - min(e.start,t.cds_start,t.cds_end) ),
                                            sum( min(e.end,max(t.cds_start,t.cds_end)) - min(e.start,max(t.cds_start,t.cds_end)) ) - 1,
                                            t.cds_start, t.cds_end-1
                                    from transcripts t, exons e
                                    where t.transcript_id=? and e.transcript=t.id""",
                        (transcript,) ).fetchone() ;

    if s==None or e==None or gs==None or ge==None:
        return (None,None,None,None) #no ORF Found
    else:
        return (s,e,gs,ge)


def get_cDNA(transcripts):
    """Get cDNA sequence for each transcript in list of transcripts.
    Uses hg38 genome 2bit.  (2bit is smaller and faster than Fasta.  A lot.)

    Args:
        transcripts (list): list of transcript ids like NM_001005484.1, NM_001005485.1

    Returns:
        (dict): maps each transcript_id to its cDNA sequence (a string)
    """
    cDNAs = {}
    revcom = str.maketrans("ACGT","TGCA")
    genome = twobitreader.TwoBitFile(genome_file)

    for current_transcript in transcripts:
        statement_get_exon_structure = "SELECT e.chrom, e.strand, e.start, e.end from transcripts as t, exons as e WHERE e.transcript==t.id and t.transcript_id==? ORDER BY e.sequence_number"

        cur_seq = []
        cur = con.execute(statement_get_exon_structure, [current_transcript])
        for chrom, strand, start, end in cur:
            if strand!='-':
                cur_seq.append(genome[chrom][start:end])
            else:
                cur_seq.insert(0, genome[chrom][start:end].translate(revcom)[::-1])
        cur.close()
        cDNAs[current_transcript] = ''.join(cur_seq)

    return cDNAs


def get_proteins(transcript_ids):
    """ Get proteins for a list of of transcripts

    Args:
        transcript_ids (list-ish): transcript ids

    Returns:
        (dict): maps transcript_ids to their protein sequences
    """
    cDNA_dict = get_cDNA(transcript_ids)

    proteins = {}
    for transcript in cDNA_dict:
        # print("transcript: %s" % transcript)
        cDNA_Seq = Seq(cDNA_dict.get(transcript)) #get cDNA of transcript
        # print("cDNA_string: %s" % cDNA_Seq)
        start_pos, end_pos, start_genome, end_genome = find_longest_ORF(transcript)
        if start_pos is None or end_pos is None:
            print("No ORF found")
            proteins[transcript]="No ORF found"
            continue
        # print("--ORF from %d to %d--" % (start_genome, end_genome))
        if start_genome<=end_genome:
            proteins[transcript] = cDNA_Seq[start_pos:end_pos+1].translate()
        else:
            s = len(cDNA_Seq) - (end_pos+1)
            e = len(cDNA_Seq) - start_pos
            proteins[transcript] = cDNA_Seq[s:e].translate()
    return proteins


def get_proteins_by_gene(gene_name):
    """ Get all proteins from all transcripts corresponding to gene_name

    Args:
        gene_name (string): gene name like "NSTRG.1" or "BRCA2"

    Raises:
        exceptions.PreventUpdate: if gene_name was not found

    Returns:
        (dict): maps transcript_ids to their protein sequences
    """
    print("---- get_proteins_by_gene from all transcript corresponding from gene_name ----")
    statement_get_all_transcripts = "SELECT transcript_id FROM transcripts WHERE gene_name=?"
    if gene_name is None:
        raise exceptions.PreventUpdate
    res = con.execute(statement_get_all_transcripts, [gene_name])
    tids = res.fetchall()

    if tids:
        return get_proteins([t for (t,) in tids])
    else:
        return (no_update, no_update, "Not available")


def get_mRNA_from_gene_name(gene_name):
    """Get all proteins from all transcripts corresponding to gene_name

    Args:
        gene_name (string): gene name like "NSTRG.1" or "BRCA2"

    Raises:
        exceptions.PreventUpdate: if gene_name was not foundd

    Returns:
        (dict): maps transcript_id to cDNA sequence
    """
    print("---- get cDNA from all transcript corresponding from gene_name ----")
    proteins = ""
    statement_get_all_transcripts = "SELECT transcript_id FROM transcripts WHERE gene_name=?"
    if gene_name is None:
        raise exceptions.PreventUpdate
    res = con.execute(statement_get_all_transcripts, [gene_name])
    tids = res.fetchall()
    if tids:
        return get_cDNA([t for (t,) in tids])
    else:
        return (no_update, no_update, "Not available")


def get_transcripts_by_gennomics_pos(chrom, start, stop, strand):
    """Get all transcripts from gennomics position
    # TODO - overlapping of genomic region
    Args:
        chrom (string): chromosome of gene
        start (string): start position on genome
        stop (string): end position on genome
        strand (string): strand of gene

    Returns:
        (DataFrame): DataFrame with all transcripts from gennomics position
    """
    if strand:
        statement = "SELECT t.gene_name, t.transcript_id FROM transcripts t where t.id in (select e.transcript from exons e where e.end>=? and e.start<=? and e.chrom=? and e.strand=?)"
        res = con.execute(statement,(start,stop,chrom,strand))
    else:
        statement = "SELECT t.gene_name, t.transcript_id FROM transcripts t where t.id in (select e.transcript from exons e where e.end>=? and e.start<=? and e.chrom=?)"
        res = con.execute(statement,(start,stop,chrom))

    df = pd.DataFrame(res.fetchall())
    df.columns = ["gene_name","transcript_id"]
    return df

def get_TPM_from_tissues(gene_name, tissues):
    """Get TPMs for all transcripts in given tissues for given gene name

    Args:
        gene_name (string): gene name like "NSTRG.1" or "BARD1"
        tissues (list): list of tissues (strings)

    Returns:
        (DataFrame): DataFrame with TPMs, mean and standard deviation for each transcript
    """
    print("--Get TPMs for all tissues--")
    print("tissues")
    #print current time
    print("Calculate TPMs for each tissue" ,datetime.datetime.now())

    # XXX  Shouldn't this whole loop be doable in a single SQL statement?
    df_result = pd.DataFrame()
    i = 0
    for tissue in tissues:
        statement = """SELECT t.gene_name, t.transcript_id,
                              AVG(getTpm(s.expresses_off,s.expresses_num,t.id)),
                              stdev(getTpm(s.expresses_off,s.expresses_num,t.id)) as tpm
                       FROM samples AS s, transcripts AS t
                       WHERE s.tissue=? AND t.gene_name=? GROUP BY t.transcript_id"""
        res = con.execute(statement, (tissue, gene_name))
        df_A = pd.DataFrame(res.fetchall())
        df_A.columns = ["gene_name","transcript_id","TPM(mean)"+tissue, "TPM(sd)"+tissue]

        #calculate the percentage of the TPM of the gene
        sum_gene_TPM = df_A["TPM(mean)"+tissue].sum()
        tpm_percentage = "TPM(%)"+tissue
        df_A[tpm_percentage] = (df_A["TPM(mean)"+tissue]/sum_gene_TPM)*100

        #if NAN in df_A then replace with 0
        df_A = df_A.fillna(0)
        df_A = df_A.round({"TPM(mean)"+tissue: 3, tpm_percentage:1, "TPM(sd)"+tissue: 3})
        df_A = df_A.reindex(columns = ["gene_name","transcript_id","TPM(mean)"+tissue, "TPM(sd)"+tissue, "TPM(%)"+tissue])

        if i==0:
            df_result = df_A
            df_result.columns = ["gene_name","transcript_id","TPM(mean)"+tissue, "TPM(sd)"+tissue, "TPM(%)"+tissue]
        else:
            df_result = df_result.merge(df_A)
        i += 1

    #sort df_result by transcript_id
    df_result = df_result.sort_values(by=["transcript_id"])

    # now get the length, number of exons, start and stop of each transcript

    length_transcript = []
    start = []
    end = []
    number_exons = []
    for transcript in  df_result["transcript_id"]: #for each transcript get length, number of exons, start, stop
        #reused code from get exon structure callback to get exon structure information about transcripts
        statement_get_exon_structure = "SELECT e.sequence_number, e.start, e.end from transcripts as t, exons as e WHERE e.transcript==t.id and t.transcript_id==?"
        res = con.execute(statement_get_exon_structure, (transcript,))
        df = pd.DataFrame(res.fetchall())
        columns = ["exon number", "start", "end"]
        df.columns = columns

        number_of_exon = df.iloc[len(df)-1][0] #number of exons
        first_exon_start = df.iloc[0][1] #get position of first exon of
        last_exon_stop = df.iloc[len(df)-1][2] #get position of last exon
        length = last_exon_stop - first_exon_start #get length

        length_transcript.append(length)
        start.append(first_exon_start)
        end.append(last_exon_stop)
        number_exons.append(number_of_exon)

    # TODO: check gene RBFOX2
    predicted_proteins = get_proteins(df_result["transcript_id"])
    # for each key in dict get the length of the proteins which is coded as the length of the Seq object
    length_protein = []
    for transcripts_id in predicted_proteins:
        length_protein.append(len(predicted_proteins[transcripts_id]))
    df_result.insert(2, "length transcript (bp)", length_transcript)
    df_result.insert(3, "length predicted protein (as)", length_protein)
    df_result.insert(4, "start", start)
    df_result.insert(5, "end", end)
    df_result.insert(6, "# of exons", number_exons)
    
    return df_result

def get_TPM_from_tissues_over_transcripts(transcripts, tissues):
    """Get TPMs for all transcripts in given tissues for given list of transcripts

    Args:
        transcripts (list): list of transcript_id
        tissues (list): list of tissues (strings)

    Returns:
        (DataFrame): DataFrame with TPMs, mean and standard deviation for each transcript
    """
    print("--Get TPMs for all tissues--")
    print("tissues")
    #print current time
    print("Calculate TPMs for each tissue" ,datetime.datetime.now())
    df_result = pd.DataFrame()
    i = 0
    for tissue in tissues:
        print(tissue)
        statement = """SELECT t.gene_name, t.transcript_id,
                              AVG(getTpm(s.expresses_off, s.expresses_num, t.id)) AS avg_tpm,
                              STDEV(getTpm(s.expresses_off, s.expresses_num, t.id)) AS stdev_tpm
                       FROM transcripts t, samples s
                       WHERE t.transcript_id IN (%s) AND s.tissue=?
                       GROUP BY t.gene_name, t.transcript_id""" % ",".join(["?"] * len(transcripts))
        res = con.execute(statement, [ t for t in transcripts ] + [tissue])
        df_A = pd.DataFrame(res.fetchall())
        df_A.columns = ["gene_name","transcript_id","TPM(mean)"+tissue, "TPM(sd)"+tissue]

        #calculate the percentage of the TPM of the gene
        sum_gene_TPM = df_A["TPM(mean)"+tissue].sum()
        tpm_percentage = "TPM(%)"+tissue
        df_A[tpm_percentage] = (df_A["TPM(mean)"+tissue]/sum_gene_TPM)*100

        #if NAN in df_A then replace with 0
        df_A = df_A.fillna(0)
        df_A = df_A.round({"TPM(mean)"+tissue: 3, tpm_percentage:1, "TPM(sd)"+tissue: 3})
        df_A = df_A.reindex(columns = ["gene_name","transcript_id","TPM(mean)"+tissue, "TPM(sd)"+tissue, "TPM(%)"+tissue])

        if i==0:
            df_result = df_A
            df_result.columns = ["gene_name","transcript_id","TPM(mean)"+tissue, "TPM(sd)"+tissue, "TPM(%)"+tissue]
        else:
            df_result = df_result.merge(df_A)
        i += 1

    #sort df_result by transcript_id
    df_result = df_result.sort_values(by=["transcript_id"])

    # now get the length, number of exons, start and stop of each transcript
    length_transcript = []
    start = []
    end = []
    number_exons = []
    for transcript in df_result["transcript_id"]: #for each transcript get length, number of exons, start, stop
        #reused code from get exon structure callback to get exon structure information about transcripts
        statement_get_exon_structure = "SELECT e.sequence_number, e.start, e.end from transcripts as t, exons as e WHERE e.transcript==t.id and t.transcript_id==?"
        res = con.execute(statement_get_exon_structure, (transcript,))
        df = pd.DataFrame(res.fetchall())
        columns = ["exon number", "start", "end"]
        df.columns = columns

        number_of_exon = df.iloc[len(df)-1][0] #number of exons
        first_exon_start = df.iloc[0][1] #get position of first exon of
        last_exon_stop = df.iloc[len(df)-1][2] #get position of last exon
        length = last_exon_stop - first_exon_start #get length

        length_transcript.append(length)
        start.append(first_exon_start)
        end.append(last_exon_stop)
        number_exons.append(number_of_exon)

    #insert value at column 1
    df_result.insert(2, "length (bp)", length_transcript)
    df_result.insert(3, "start", start)
    df_result.insert(4, "end", end)
    df_result.insert(5, "# of exons", number_exons)

    return df_result

def get_group_comparisons_from_gene(gene_name, groupA, groupB):
    """Get TPMs for all transcripts in given groupA and groupB for given gene name

    Args:
        gene_name (string): gene name like "NSTRG.1" or "BARD1"
        groupA (list): list of samples (strings) or tissues (strings)
        groupB (list): list of samples (strings) or tissues (strings)

    Returns:
        (DataFrame): DataFrame with TPMs, mean and standard deviation for each transcript for groupA and groupB
    """

    statement = """SELECT t.gene_name, t.transcript_id,
                            AVG(getTpm(s.expresses_off, s.expresses_num, t.id)),
                            STDEV(getTpm(s.expresses_off, s.expresses_num, t.id))
                    FROM samples AS s, transcripts AS t
                    WHERE s.%s IN (%s) AND t.gene_name=?
                    GROUP BY t.transcript_id"""

    if (groupA[0]).startswith("SRR") and (groupB[0]).startswith("SRR") :
        print("SRAs")
        print("Calculate TPMs over groups" ,datetime.datetime.now())
        statement_A = statement % ("name", ",".join(["?"] * len(groupA)))
        statement_B = statement % ("name", ",".join(["?"] * len(groupB)))
    else:
        print("tissues")
        #print current time
        print("Calculate TPMs over groups" ,datetime.datetime.now())
        statement_A = statement % ("tissue", ",".join(["?"] * len(groupA)))
        statement_B = statement % ("tissue", ",".join(["?"] * len(groupB)))

    groupA.append(gene_name)
    res = con.execute(statement_A, (groupA))
    print("group: ",groupA ,datetime.datetime.now())
    df_A = pd.DataFrame(res.fetchall())
    df_A.columns = ["gene_name","transcript_id","TPM(mean)A", "TPM(sd)A"]

    groupB.append(gene_name)
    res = con.execute(statement_B, (groupB))
    print("group: ",groupB ,datetime.datetime.now())
    df_B = pd.DataFrame(res.fetchall())
    df_B.columns = ["gene_name","transcript_id","TPM(mean)B", "TPM(sd)B"]
    df_result = pd.merge(df_A,df_B)

    return df_result

def get_group_comparisons_over_transcripts(transcripts, groupA, groupB):
    """Get TPMs for all transcripts in given groupA and groupB for given transcripts

    Args:
        transcripts (list): list of transcript_id
        groupA (list): list of samples (strings) or tissues (strings)
        groupB (list): list of samples (strings) or tissues (strings)

    Returns:
        (DataFrame): DataFrame with TPMs, mean and standard deviation for each transcript for groupA and groupB
    """

    statement = """SELECT t.gene_name, t.transcript_id,
                            AVG(getTpm(s.expresses_off,s.expresses_num,t.id)) AS avg_tpm,
                            STDEV(getTpm(s.expresses_off,s.expresses_num,t.id)) AS stdev_tpm
                    FROM transcripts AS t, samples AS s
                    WHERE s.%s IN (%s) AND t.transcript_id IN (%s)
                    GROUP BY t.gene_name, t.transcript_id"""

    # XX this is b0rk3n, isn't it?
    if (groupA[0]).startswith("SRR") and (groupB[0]).startswith("SRR") :
        print("SRAs")
        print("Calculate TPMs over groups" ,datetime.datetime.now())
        groupA = groupA.split(",")
        #for each element in list remove "['" and "']"
        for i in range(len(groupA)):
            groupA[i] = groupA[i].replace("['","")
            groupA[i] = groupA[i].replace("']","")
        groupB = groupB.split(",")
        #for each element in list remove "['" and "']"
        for i in range(len(groupB)):
            groupB[i] = groupB[i].replace("['","")
            groupB[i] = groupB[i].replace("']","")

        #  XXX  This code is untested, because I couldn't reach it.
        statement_A = statement % ("name", ",".join(["?"] * len(groupA)), ",".join(["?"] * len(transcripts)))
        statement_B = statement % ("name", ",".join(["?"] * len(groupB)), ",".join(["?"] * len(transcripts)))
    else:
        print("tissues")
        #print current time
        print("Calculate TPMs over groups" ,datetime.datetime.now())
        groupA = groupA.split(",")
        #for each element in list remove "['" and "']"
        for i in range(len(groupA)):
            groupA[i] = groupA[i].replace("['","")
            groupA[i] = groupA[i].replace("']","")
        groupB = groupB.split(",")
        #for each element in list remove "['" and "']"
        for i in range(len(groupB)):
            groupB[i] = groupB[i].replace("['","")
            groupB[i] = groupB[i].replace("']","")

        statement_A = statement % ("tissue", ",".join(["?"] * len(groupA)), ",".join(["?"] * len(transcripts)))
        statement_B = statement % ("tissue", ",".join(["?"] * len(groupB)), ",".join(["?"] * len(transcripts)))

    res = con.execute(statement_A, groupA + transcripts)
    print("group: ", groupA ,datetime.datetime.now())
    df_A = pd.DataFrame(res.fetchall())
    df_A.columns = ["gene_name","transcript_id","TPM(mean)A", "TPM(sd)A"]

    res = con.execute(statement_B, groupB + transcripts)
    print("group: ",groupB ,datetime.datetime.now())
    df_B = pd.DataFrame(res.fetchall())
    df_B.columns = ["gene_name","transcript_id","TPM(mean)B", "TPM(sd)B"]
    df_result = pd.merge(df_A,df_B)
    return df_result

def calculate_percentage_for_TPM(df_result):
    """Calculate percentage of TPMs for each transcript in groupA and groupB

    Args:
        df_result (DataFrame): DataFrame with TPMs, mean and standard deviation for each transcript for groupA and groupB

    Returns:
        (DataFrame): DataFrame extended for percentage of TPMs
    """
    length_transcript = []
    start = []
    end = []
    number_exons = []
    for transcript in df_result["transcript_id"]: #for each transcript get length, number of exons, start, stop
        #reused code from get exon structure callback to get exon structure information about transcripts
        statement_get_exon_structure = "SELECT e.sequence_number, e.start, e.end from transcripts as t, exons as e WHERE e.transcript==t.id and t.transcript_id==?"
        res = con.execute(statement_get_exon_structure, (transcript,))
        df = pd.DataFrame(res.fetchall())
        columns = ["exon number", "start", "end"]
        df.columns = columns

        number_of_exon = df.iloc[len(df)-1][0] #number of exons
        first_exon_start = df.iloc[0][1] #get position of first exon of
        last_exon_stop = df.iloc[len(df)-1][2] #get position of last exon
        length = last_exon_stop - first_exon_start #get length

        length_transcript.append(length)
        start.append(first_exon_start)
        end.append(last_exon_stop)
        number_exons.append(number_of_exon)

    #insert length of transcript
    df_result['length (bp)'] = length_transcript
    df_result['start'] = start
    df_result['end'] = end
    df_result['# of exons'] = number_exons

    #insert averages into df_A and df_B and calcualte percentages
    sum_gene_TPM_A = df_result['TPM(mean)A'].sum()
    df_result['TPM(%)A'] = (df_result['TPM(mean)A']/sum_gene_TPM_A)*100
    sum_gene_TPM_B = df_result['TPM(mean)B'].sum()
    df_result['TPM(%)B'] = (df_result['TPM(mean)B']/sum_gene_TPM_B)*100

    #round TPM values
    df_result = df_result.round({'TPM(mean)A': 3, 'TPM(%)A':1,'TPM(mean)B': 3, 'TPM(%)B':1, 'TPM(sd)A':3, 'TPM(sd)B':3})
    # sort columns of df_result by 'TPM(mean)A', 'TPM(%)A', 'TPM(mean)B', 'TPM(%)B'
    df_result = df_result.reindex(['gene_name','transcript_id','length (bp)','start','end','# of exons','TPM(mean)A','TPM(sd)A','TPM(%)A', 'TPM(mean)B','TPM(sd)B','TPM(%)B'], axis=1)
    return (df_result)

def calculate_inner_exons(gene_name):
    """Calculate the number of distict inner exons for a gene.
    Each transcript contributes all its exons, except the first and last.  Exons with identical coordinates count only once.
    this was a wish from Torsten Schoeneberg

    XXX Isn't this a job for SQL?

    Args:
        gene_name (string): name of gene
    """

    statement = "SELECT t.transcript_id, e.sequence_number, e.start, e.end FROM exons as e, transcripts as t WHERE e.transcript == t.id AND t.gene_name = ?"
    rows = con.execute(statement, (gene_name,))
    df = pd.DataFrame(rows, columns=['transcript_id','sequence_number', 'start', 'end'])
    #generate for each new transcript_id in df a new dataframe
    result_df = []
    for transcript_id in df['transcript_id'].unique():
        #get subset of df for transcript_id
        df_subset = df[df['transcript_id'] == transcript_id]
        #delete first and last entry in df_subset
        df_subset = df_subset.iloc[1:-1]
        #add df_subset to result_df
        result_df.append(df_subset)

    #sort result_df for start and end
    result_df = pd.concat(result_df)
    result_df = result_df.sort_values(by=['start', 'end'], ascending=True)
    #delete row if start and end numbers are the same as for another row in result_df
    result_df = result_df.drop_duplicates(subset=['start', 'end'], keep='first')
    number_inner_exons = result_df.shape[0]
    return number_inner_exons

################ 2.0 DASH APP #########################
# TODO - use dbc.Col and dbc.Row to make the app responsive
#APP Instance
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP,'assets/styles.css'], url_base_pathname=url_base_pathname)

# initialize the dropdowns in the app
ref_gene_names, tissues_and_samples, tissues = fill_dropdowns_in_dash_application()

#APP Layout
# the app has 4 cards to far
card_explanation = dbc.Card([
    dbc.CardHeader("How-to-use", style={"font-size": "24px"}),
    dbc.CardBody([
        dbc.Row([
            dbc.Col([
                html.P("With this app the user can perform qualitative and quantitative analyses of human transcript variants in 48 different deeply sequenced tissue types."),
            ], width=12),
            dbc.Col([
                html.H6("How to use the app: "),
                html.P(["The app is divided into 2 steps",
                        html.Br(),
                        " Step 1: Select Tissue/Samples for group comparison or across all tissues",
                        html.Br(),
                        " Step 2: Select Gene for visualization and analysis of transcripts"]),
            ], width=5),
            dbc.Col([
                html.H6("Video Tutorial: "),
                html.A(
                    html.Img(src=app.get_asset_url('VideoExplanation.png'), alt='Video Tutorial', style={'width':'100%'}),
                    href='https://youtu.be/fMvIgQYDtxs',
                    target='_blank'  # Open the link in a new tab
                )
            ], width=2),
            dbc.Col([
                html.H6("More explanations on all functions (press link): "),
                html.Ul([
                    html.Li([
                        html.A(
                            "How to use the app",
                            href=app.get_asset_url("download/HowTo.png"),
                            target="_blank",
                        )
                    ]),
                    html.Li([
                        html.A(
                            "Exportable graphics and tables.",
                            href=app.get_asset_url("download/GraphicsAndTables.png"),
                            target="_blank",
                        )
                    ]),
                    html.Li([
                        html.A(
                            "Heat-maps for comparisons of tissue specific-expression",
                            href=app.get_asset_url("download/Heatmaps.png"),
                            target="_blank",
                        )
                    ]),
                    html.Li([
                        html.A(
                            "Projection of protein domains on exons",
                            href=app.get_asset_url("download/Domains.png"),
                            target="_blank",
                        )
                    ]),
                    html.Li([
                        html.A(
                            "Projection of disease- or phenotype-causing variants on transcript’s exon structure",
                            href=app.get_asset_url("download/Mutations.png"),
                            target="_blank",
                        )
                    ]),
                    html.Li([
                        html.A(
                            "All genomic positions of the exons contributing to a transcript",
                            href=app.get_asset_url("download/ExonStructure.png"),
                            target="_blank",
                        )
                    ]),
                ])

            ], width=5),
            dbc.Col([
                html.P([
                    "If you use Splice-O-Mat for a scientific publication, please cite \"The repertoire and structure of adhesion GPCR transcript variants assembled from publicly available deep-sequenced human samples\", Christina Katharina Kuhn, Udo Stenzel, Sandra Berndt, Ines Liebscher, Torsten Schöneberg, Susanne Horn, Nucleic Acids Research, gkae145, ",
                    html.A("https://doi.org/10.1093/nar/gkae145", href="https://doi.org/10.1093/nar/gkae145")
                    ]),
            ], width=12)
        ]),
    ])
])

popover_step1_option1_explanations = "The user can select from tissue types (e.g. liver, brain...) or all tissues at once 'Across all tissues'."
popover_step1_option2_explanations  = "The user can select one or more samples (SRA...) or tissue types (e.g. liver, brain) for group comparisons (group A vs. group B)."

tab1_content_sample_selection = dbc.Card(
    dbc.CardBody(
        [
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dbc.Button(
                                "Info",
                                id="hover-target-option1",
                                color="info",
                                className="me-1",
                                n_clicks=0,
                            ),
                            dbc.Popover(
                                popover_step1_option1_explanations,
                                target="hover-target-option1",
                                body=True,
                                trigger="hover",
                                style={"font-size": "18px", "min-width": "50%"}
                            ),
                            html.Br(),
                            html.Br(),
                            dbc.Row(
                                [
                                    dbc.Col(
                                        [
                                            dcc.Dropdown(
                                                options=tissues,
                                                id='all-tissues-dropdown',
                                                multi=True
                                            ),
                                        ], width=6
                                    ),
                                    dbc.Col(
                                        [
                                            dbc.Button('Across all tissues', color="primary", id='all-tissues-button',
                                                       n_clicks=0),
                                        ], width=4
                                    ),
                                ]
                            ),
                            dbc.Row(
                                [
                                    dbc.Col(
                                        [
                                            html.Br(),
                                            dbc.Button("Download TPM without mean", id="btn-download-without-groups",
                                                       n_clicks=0,
                                                       style={
                                                           'background-color': 'rgba(211, 211, 211, 0.6)',
                                                           'border-color': 'rgba(211, 211, 211, 0.6)',
                                                           'color': 'grey'}
                                                       ),
                                            dcc.Download(id="download-without-groups"),
                                        ], width=6
                                    )
                                ]
                            )
                        ]
                    ),
                ]
            ),
            dbc.Row([html.Br()]),
        ]
    ),
    className="mt-3",
)


tab2_content_sample_selection = dbc.Card(
    [
        dbc.CardBody(
            [
                dbc.Col(
                    [
                        dbc.Button(
                            "Info",
                            id="hover-target-option2",
                            color="info",
                            className="me-1",
                            n_clicks=0,
                        ),
                        dbc.Popover(
                            popover_step1_option2_explanations,
                            target="hover-target-option2",
                            body=True,
                            trigger="hover",
                            style={"font-size": "18px", "min-width": "50%"}
                        ),
                        html.Br(),
                        html.Br(),
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        dbc.Label('Group A'),
                                        dcc.Dropdown(
                                            options=tissues_and_samples,
                                            id='group-comparisonA',
                                            multi=True
                                        ),
                                        html.Br(),
                                        dbc.Label('Chosen group A'),
                                        html.P("none", id='groupA', style={'color': 'blue'}),
                                    ]
                                ),
                                dbc.Col(
                                    [
                                        dbc.Label('Group B'),
                                        dcc.Dropdown(
                                            options=tissues_and_samples,
                                            id='group-comparisonB',
                                            multi=True
                                        ),
                                        html.Br(),
                                        dbc.Label('Chosen group B'),
                                        html.P("none", id='groupB', style={'color': 'blue'}),
                                    ]
                                ),
                            ]
                        ),
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        dbc.Button('Select Groups', id='group-button', n_clicks=0),
                                        html.Br(),
                                    ]
                                ),
                            ]
                        ),
                    ]
                ),
            ]
        ),
        dbc.Row([html.Br()]),
    ],
    className="mt-3",
)

card0 = dbc.Card([
    dbc.CardHeader("Step 1: Sample/Tissue selection", style={"font-size": "24px"}),
    dbc.Tabs(
                [
                    dbc.Tab(tab1_content_sample_selection, label="Option 1: Across all tissues/samples"),
                    dbc.Tab(tab2_content_sample_selection, label="Option 2: Group comparisons"),
                ]
            ),
])

popover_step2_explanations = "The user can select a gene (based on the annotation of the NCBI hg38 human genome).\
    The plot and the output table show all transcripts (transcript_id) from this gene (assigned by stringtie over gene_name) \
    including known transcripts (RefSeq accessions NM, NR, XM, XR) and new transcripts (NSTRG, not already annotated in the NCBI annotation) of the selected gene. \
    NSTRG's were assigned to the gene derived from the annotated transcripts with which they overlapped the most. \
    This was done counting overlapping exonic bases."


# card1
card1 = dbc.Card([
    dbc.CardHeader("Step 2: Discover transcripts of a gene", style={"font-size": "24px"}),
    dbc.CardBody([
        dbc.Button(
            "Info",
            id="hover-target-step2",
            color="info",
            className="me-1",
            n_clicks=0,
        ),
        dbc.Popover(
            popover_step2_explanations,
            target="hover-target-step2",
            body=True,
            trigger="hover",
            style={"font-size": "18px", "min-width": "50%"}
            ),
        html.Br(),
        html.Br(),
        dbc.Row([
            dbc.Col([
                html.H6('Select a gene and search for transcripts:'),
            ]),
        ]),
        dbc.Row([
            dbc.Col([
                dcc.Dropdown(id='my-dynamic-dropdown'),
                ], width=6),
            dbc.Col([
                dbc.Button('Search transcripts', id='transcript-button', color="primary", n_clicks=0),
            ], width=4),
            dbc.Col([
            ], width=2),
        ]),
        html.Br(),
        dbc.Row([
            dbc.Col([html.P(id='no-groups-warning',style={'color':'red'}),])
        ]),
        #position part
        dbc.Row([
            dbc.Col([
                html.H5('Output:'),
                html.Br(),
            ], width=4),
        ]),
        dbc.Row([
            dbc.Col([
                html.Div([dbc.Label("chrom: "),
                dbc.Input(id='chrom', value='', type='text'),])
            ],width=2),
            dbc.Col([
                html.Div([dbc.Label("start: "),
                dbc.Input(id='start', value='', type='text'),])
            ],width=2),
            dbc.Col([
                html.Div([dbc.Label("stop: "),
                dbc.Input(id='stop', value='', type='text'),])
            ],width=2),
            dbc.Col([
                html.Div([dbc.Label("strand: "),
                dbc.Input(id='strand', value='', type='text'),])
            ],width=2),
            dbc.Col([
                html.Div([dbc.Label("genomic positions to depict (positions separated by comma):"),
                dbc.Input(id='mutation', placeholder='position', type='text'),])
            ],width=2),
            dbc.Col([
                html.Br(),
                dbc.Button('Update transcripts', id='update-button', color="secondary", n_clicks=0)
            ],width=2),
        ]),
        #search transcripts and get cDNA
        html.Br(),
        dcc.Loading( type="default",children=[
            dbc.Row([
            dbc.Col(dbc.Button('Get mRNA', id='mRNA-button', color="secondary", n_clicks=0), width=2) ,
            dbc.Col(dbc.Button('Get proteins', id='protein-button', color="secondary", n_clicks=0), width=2) ,
            dbc.Col(dbc.Button('Get domains', id='domains-button',  color="secondary",n_clicks=0), width=2) ,
            dbc.Col(dcc.Loading(
                type="default", children=[
                    dcc.Download(id='download-cDNA'),
                    dcc.Download(id='download-proteins'),
                    dcc.Download(id='download-domains'),
                ]), width=0),
            ]),
            dbc.Row([
                dbc.Col([html.Br()]),
                dbc.Col([html.Div(html.Img(id="gene-png", width="100%"))],width=12),
            ]),
            # dbc.Row([
            #     dbc.Col([
            #         html.Div(
            #             dbc.Alert("Attention: For small display size or for high number of transcript variants, heatmap might show only every second transcript and tissue. Please hover over the heatmap for more informtion.", color="secondary"),
            #             id="display-size-warning")
            #     ],width=12),
            # ]),
            dbc.Row([
                dbc.Col([
                    # html.Div(id="heatmap-relatives")
                    html.Div([dcc.Graph(id="heatmap-relatives", responsive=False)], id="heatmap-relatives-div", style={'display':'none'})
                    # html.Div(html.Img(id="heatmap-relatives", width="100%"))
                ],width=12),

            ]),
            dbc.Row([html.Br(),]),
            dbc.Row([
                dbc.Col([
                    # html.Div(id="heatmap-absolutes")
                    html.Div([dcc.Graph(id="heatmap-absolutes", responsive=False)], id="heatmap-absolutes-div", style={'display':'none'})
                    # html.Div(html.Img(id="heatmap-absolutes", width="100%"))
                ],width=12),
            ]),
            dbc.Row([html.Br(),]),
            dbc.Row([
                dbc.Col([
                    dbc.Button("Export", id="btn",  color="secondary", n_clicks=0),
                    dcc.Download(id="download"),
                ], width=2),
            ]),
            dbc.Row([
                dbc.Col([
                    html.Div(
                    dash_table.DataTable(
                        id='search-output-ref-geneA',
                        columns=[],
                        data=[],
                        sort_action='native',  # Enable native sorting
                        style_cell={'fontSize':14},

                        )
                    ),
                ], width=12),
            ]),
        ]),
    ])
])

popover_gene_name_explanations = "The user can enter a gene_name. The output table shows all transcripts assigned to this gene_name."
popover_transcripts_by_region= "The user can specify a region (chromosome, start [bp], stop [bp], strand; annotation used: NCBI hg38 human genome). The output table shows all transcripts that have an exon overlapping this region (see figure)."

tab1_content_transcripts = dbc.Card(
    dbc.CardBody([
            dbc.Row([
            html.H6('Show transcripts by gene name'),
            ]),
            dbc.Button(
                        "Info",
                        id="hover-target-gene-id",
                        color="info",
                        className="me-1",
                        n_clicks=0,
                    ),
                    dbc.Popover(
                        popover_gene_name_explanations,
                        target="hover-target-gene-id",
                        body=True,
                        trigger="hover",
                        style={"font-size": "18px", "min-width": "50%"}
                    ),
            html.Br(),
            html.Br(),
            dbc.Row([
                dbc.Col([dbc.Input(type="text", id='gene-id',  placeholder='NSTRG.8029'),], width=6),
                dbc.Col([html.P(id='gene-id-err',style={'color':'red'}),], width=6),
            ]),
            dbc.Row([
                dbc.Col([html.Br(),]),
            ]),
            dbc.Row([
                dbc.Col([
                    dash_table.DataTable(id='search-output-gene-id',columns=[],data=[], sort_action='native',export_format="csv"),
                ]),
            ]),
            html.Br(),
         ]
    ),
    className="mt-3",
)

tab2_content_transcripts = dbc.Card(
    dbc.CardBody(
        [
        dbc.Row([
                html.Br(),
                html.H6('Show transcripts that overlap genomic region'),
            ]),
        dbc.Button(
                        "Info",
                        id="hover-target-region",
                        color="info",
                        className="me-1",
                        n_clicks=0,
                    ),
                    dbc.Popover(
                        popover_transcripts_by_region,
                        target="hover-target-region",
                        body=True,
                        trigger="hover",
                        style={"font-size": "18px", "min-width": "50%"}
                        ),
        html.Br(),
        html.Br(),
        dbc.Row([
                html.Div(html.Img(src=app.get_asset_url('AbbildungExon.png'), alt='image', style={'width':'20%'})),
                html.Br(),
            ]),
            dbc.Row([
                dbc.Col([
                    html.Div([dbc.Label("chromosom: "),
                    dbc.Input(id='text-chro', value='chr10', type='text'),]),
                ], width=6),
            ]),
            dbc.Row([
                dbc.Col([
                    html.Div([dbc.Label("start: "),
                    dbc.Input(id='text-start', value='133086592', type='text'),]),
                ], width=6),
            ]),
            dbc.Row([
                dbc.Col([
                    html.Div([dbc.Label("stop: "),
                    dbc.Input(id='text-stop', value='133131675', type='text'),]),
                ], width=6),
            ]),
            dbc.Row([
                dbc.Col([
                    html.Div([dbc.Label("strand: "),
                    dbc.Input(id='text-strand', value='+', type='text'),]),
                ], width=6),
            ]),
            dbc.Row([
                dbc.Col([
                    html.Br(),
                    dbc.Button('Find transcripts', id='my-button', n_clicks=0),
                    html.P(id='output-exon-err',style={'color':'red'}),
                    html.Br(),
                ], width=6),
            ]),
            dbc.Row([
                dbc.Col([
                    dbc.Button("Export", id="btn-exon",  color="secondary", n_clicks=0),
                    dcc.Download(id="download-exon"),
                ], width=2)
            ]),
            dbc.Row([
                dash_table.DataTable(id='search-output-exon',columns=[],data=[],sort_action='native'),
                html.Br(),
            ]),
        ]),
    className="mt-3",
)

card_1a = dbc.Card([
    dbc.CardHeader("Search Transcripts by Gene ID or Genomic Region", style={"font-size": "24px"}),
    dbc.Tabs(
                [
                    dbc.Tab(tab1_content_transcripts, label="Transcripts by Gene ID"),
                    dbc.Tab(tab2_content_transcripts, label="Transcripts that overlap genomic region"),
                ]
            ),
])

popover_exon_structure_explanations = "The user can select a specific transcript based on the transcript_id (e.g. NSTRG.8029.1, NTRG.8029.2, NM_001291085.1). The output table will show the exons of this transcript with start and end coordinates to be used e.g. in the genome browser (hg38)."

#card 2
card2 = dbc.Card([
    dbc.CardHeader("Exon structure", style={"font-size": "24px"}),
    dbc.CardBody([
        dbc.Row([
            html.H6('Get exon structure of specific transcript'),
        ]),
        dbc.Button(
                    "Info",
                    id="hover-target-exon-structure",
                    color="info",
                    className="me-1",
                    n_clicks=0,
                ),
                dbc.Popover(
                    popover_exon_structure_explanations,
                    target="hover-target-exon-structure",
                    body=True,
                    trigger="hover",
                    style={"font-size": "18px", "min-width": "50%"}
                    ),
        html.Br(),
        html.Br(),
        dbc.Row([
            dbc.Col([
                dbc.Input(type="text", id='transcript-id', placeholder='NSTRG.8029.1',),
                html.P(id='transcript-id-err',style={'color':'red'}),
                html.Br(),
            ], width=6),
        ]),
        dbc.Row([
                dbc.Col([
                    dbc.Button("Export", id="btn-transcript-id", color="secondary", n_clicks=0),
                    dcc.Download(id="download-transcript-id"),
                ], width=2)
            ]),
        dbc.Row([
            dash_table.DataTable(id='search-output-transcript-id',columns=[],data=[],sort_action='native'),
            html.Br(),
        ]),
      ]),
])

popover_sql_statement = "The user can make custom SQL queries to the database according to the structure of the database in the given entity-relationship model."

#card 3
card3 = dbc.Card([
    dbc.CardHeader("SQL statement:", style={"font-size": "24px"}),
    dbc.CardBody([
        dbc.Row([
            html.H6('Execute SQL statement'),
        ]),
        dbc.Button(
                    "Info",
                    id="hover-target-sql",
                    color="info",
                    className="me-1",
                    n_clicks=0,
                ),
                dbc.Popover(
                    popover_sql_statement,
                    target="hover-target-sql",
                    body=True,
                    trigger="hover",
                    style={"font-size": "18px", "min-width": "50%"}
                    ),
        html.Br(),
        html.Br(),
        dbc.Row([
            dbc.Col([html.Div(html.Img(src=app.get_asset_url('Stringtie_2.db.png'), alt='image', style={'width':'100%'}))],width=6),
        ]),
        dbc.Row([
            dbc.Col([
                html.Div([
                    "SQL statement: ",dbc.Input(id='sql-input', placeholder='SELECT * FROM ...', type='text')
                    ]),
                html.Br(),
            ],width=6),
        ]),
        dbc.Row([
            dbc.Col([
                dbc.Button('Execute SQL statement', id='sql-button', n_clicks=0),
                html.Br(),
            ],width=2),
            dbc.Col([
                html.P(id='sql-warning',style={'color':'red'}),
            ],width=6)
        ]),
        dbc.Row([
            dbc.Col([
                html.Br(),
                dcc.Loading(
                    type="default",
                    children = dash_table.DataTable(id='table-sql-output',columns=[],data=[],sort_action='native', export_format="csv"),
                )
            ]),
        ]),

    ]),
])

#card 4
card4 = dbc.Card([
    dbc.CardHeader("Info", style={"font-size": "24px"}),
    dbc.CardBody([
        html.P("We obtained three brain tissue dataset (GSE173955, GSE182321, GSE101521), one liver tissue dataset (GSE174478), one heart tissue dataset (GSE165303), one kidney tissue dataset (GSE217427) and one dataset including 45 different tissues types (GSE138734) from the Gene Expression Omnibus (GEO) public dataset. We further obtained one melanoma cancer dataset (PRJEB23709) from the European Nucleotide Archive. Only paired-end RNA-seq samples were included and datasets generated without random primers were excluded. The raw data was mapped against the hg38 human genome using STAR (version 2.7.6a) with default parameters. After indexing with samtools (version 1.9) the mapped reads were assembled to transcripts and quantified by StringTie (version v2.1.3b). StringTie parameters ‘read coverage’ (-c), ‘transcript length’ (-m) and ‘bases on both sides of a junction a spliced read has to cover’ (-a) were set to minimal values in order to avoid missing transcripts and generating a bias. The parameter ‘fraction of most abundant transcript at one locus’ (-f) was lowered from default (0.01) to 0. For all other StringTie parameters default values were used. To generate a global, unified set of transcripts across all three RNA-Seq samples, StringTie merge mode with providing the reference annotation (-G) was used. Quantification of abundance of the input transcripts was then performed using parameters ‘expression estimation mode’ (-e) with parameters ‘ballgown output’ (-B) and the beforehand generated ‘reference annotation transcripts’ (-G)."),
        html.P("Annotated transcripts are labeled with their RefSeq accessions: NM = Protein-coding transcripts (usually curated), NR = Non-protein-coding transcripts, XM = Predicted model protein-coding transcript, XR = Predicted model non-protein-coding transcript. Potential new transcripts assigned by Stringtie (not in the NCBI hg38 genome annotation) are annotated with 'NSTRG'"),
        html.P("The webtool was created by Christina Kuhn. For questions and suggestions please contact christina.kuhn@medizin.uni-leipzig.de")
    ]),
])

#put cards together
app.layout = html.Div([
    #Top
    html.Br(),
    html.H1('Splice-O-Mat'),
    html.H2('A Tool for Identifying Tissue-Specific Alternative Splicing from High-Throughput Data'),
    #Tabs
    dbc.Row([dbc.Col(card_explanation, width=12)]),
    html.Br(),
    dbc.Row([dbc.Col(card0, width=12)]),
    html.Br(),
    dbc.Row([dbc.Col(card1, width=12)]),
    html.Br(),
    dbc.Row([dbc.Col(card_1a, width=12)]),
    html.Br(),
    dbc.Row([dbc.Col(card2, width=12)]),
    html.Br(),
    dbc.Row([dbc.Col(card3, width=12)]),
    html.Br(),
    dbc.Row([dbc.Col(card4, width=12)]),
], style={'margin':50})

############# 3.0 Callbacks ############################
#TODO - Callbacks
# [ ] important - describe callbacks in comments
# [ ] important - explain input and output of callbacks
# [ ] future - sort callback by cards
# [ ] future - make each line not longer than 80 characters
# * Highlight comments
# ! alert comment
# ? question comment

@app.callback(
   [Output('groupA','children'),
    Output('groupB','children')],
    Input('group-button', 'n_clicks'),
    [State('group-comparisonA', 'value'),
     State('group-comparisonB','value')]
)
def selectGroups(n_clicks, group_comparisonA, group_comparisonB):
    """Select groups (tissue or samples) for comparison.

    Args:
        n_clicks (int): buton click
        group_comparisonA (string): dropdown value
        group_comparisonB (string): dropdown value

    Returns:
        string, string: groupA, groupB as like '["tissue1", "tissue2"]'
    """
    if (n_clicks is None) or (n_clicks == 0):
        return (no_update, no_update)
    if n_clicks is not None and n_clicks>0:
        print("----select groups----")
        if isinstance(group_comparisonA, str):
            groupA = []
            tmp_A = groupA + [group_comparisonA]
            group_comparisonA = tmp_A
        if isinstance(group_comparisonB, str):
            groupB = []
            tmp_B = groupB + [group_comparisonB]
            group_comparisonB = tmp_B
        print("Groups chosen:")
        print("GroupA: ", str(group_comparisonA))
        print("GroupB: ", str(group_comparisonB))
        return (str(group_comparisonA), str(group_comparisonB))
    else:
        return (no_update, no_update)

@app.callback(
        Output("download-without-groups", "data"),
        [Input("btn-download-without-groups", "n_clicks")],
        [State('all-tissues-dropdown','value'),
         State('my-dynamic-dropdown', 'value')],
         prevent_initial_call=True
)
def download_table_TPMS_without_means(n_clicks, all_tissues, input_value):
    """Download table with TPMs for each sample.

    Args:
        n_clicks (int): number of clicks from button btn-download-without-groups
        all_tissues (list): list of tissues from all-tissues-dropdown
        input_value (string): gene name from my-dynamic-dropdown

    Raises:
        exceptions.PreventUpdate: prevent update when no gene is found

    Returns:
        dict: dict(content=df_final.to_csv(sep="\t", index=False), filename=input_value+".tsv")
    """
    print("----download table TPMs without means----")
    df = pd.DataFrame() #dataframe for each group TPMs for each sample

    statement = "SELECT DISTINCT transcript_id FROM transcripts WHERE gene_name=?"
    transcript_ids = [t for (t,) in con.execute(statement, [input_value])]
    if not transcript_ids:
        print("no Update")
        raise exceptions.PreventUpdate

    df_final = []
    columns_list = []
    for tissue in tissues:
        # print(tissue['value'])
        df_tissue = pd.DataFrame(columns=columns_list)
        #get value at key "label" from dictionary
        for transcript in transcript_ids:
            statement = "SELECT s.tissue, s.name, getTpm(s.expresses_off, s.expresses_num, t.id) \
                         FROM samples as s, transcripts as t WHERE s.tissue=? AND t.transcript_id=?"
            res = con.execute(statement, (tissue['value'], transcript))
            df_result = pd.DataFrame(res.fetchall())
            df_result.columns= ['tissue', 'sample', 'TPM']

            #append a column transcript_tissue to df_result
            df_tissue[transcript+'_tissue'] = df_result['tissue']
            #append a column transcript_sample to df_result
            df_tissue[transcript+'_sample'] = df_result['sample']
            #append a column transcript_TPM to df_result
            df_tissue[transcript+'_TPM'] = df_result['TPM']

        #append df_tissue to df_final
        df_final.append(df_tissue)
    #drop every column including "tissue" in the name except the first one
    df_final = pd.concat(df_final)
    df_final.drop(df_final.filter(regex='tissue').columns[1:], axis=1, inplace=True)
    #drop every column including "sample" in the name except the first one
    df_final.drop(df_final.filter(regex='sample').columns[1:], axis=1, inplace=True)
    #rename first column to "tissue"
    df_final.rename(columns={df_final.columns[0]: "tissue"}, inplace=True)
    #rename second column to "sample"
    df_final.rename(columns={df_final.columns[1]: "sample"}, inplace=True)
    # TODO - calculate for each transcript kruskal-wallis test between all tissues and save p-value in a new column
    # ? hier aber eh an falscher stelle
    #calculate for each transcript kruskal-wallis test between all tissues and save p-value in a new column
    # for transcript in transcript_ids:
    #     df_final[transcript+'_p-value'] = df_final.groupby('tissue')[transcript+'_TPM'].apply(lambda x: stats.kruskal(*[group for name, group in x.groupby('tissue')])[1])

    return dict(content=df_final.to_csv(sep="\t", index=False), filename=input_value+".tsv")

@app.callback(
    Output("all-tissues-dropdown", "value"),
    Input("all-tissues-button", "n_clicks"),
    prevent_initial_call=True,
)
def selectAllTissues(n_clicks):
    """Select all tissues when clicking button all-tissues-button.

    Args:
        n_clicks (int): number of clicks from button all-tissues-button

    Returns:
        list: alle tissues in the tissues list (global variable)
    """
    if (n_clicks is None) or (n_clicks == 0):
        return (no_update, no_update)
    if n_clicks is not None and n_clicks>0:
        print("----select all tissues----")
        return [option["value"] for option in tissues]
    else:
        return (no_update)


@app.callback(
        Output("download", "data"),
        [Input("btn", "n_clicks")],
        [State("search-output-ref-geneA", "data")],
        prevent_initial_call=True
)
def download_table_transcripts(n_clicks: int, data: list[dict[str, any]]) -> dict[str, str]:
    """Downloads the data in search-output-ref-geneA as a tsv file
    by clicking on the button btn.

    Args:
        n_clicks (int): number of clicks of button btn
        data (list[dict]): table in search-output-ref-geneA

    Returns:
        dict: dict(content=df.to_csv(sep="\t", index=False),filename="data.tsv")

    Raises:
        ValueError: If data is not in the expected format or 'gene_name' column is missing.
    """
    try:
        df = pd.DataFrame.from_records(data)
    except ValueError:
        raise ValueError("Data is not in the expected format (expected a list of dictionaries).")

    if 'gene_name' not in df.columns:
        raise ValueError("'gene_name' column is missing in the data.")

    df = pd.DataFrame.from_records(data)
    gene_names = df["gene_name"].unique()
    gene_names = "_".join(gene_names)
    return dict(content=df.to_csv(sep="\t", index=False),filename=gene_names+"_data.tsv")

@app.callback(
        Output("download-exon", "data"),
        [Input("btn-exon", "n_clicks")],
        [State("search-output-exon", "data")],
        prevent_initial_call=True
)
def download_table_exons(n_clicks, data):
    """Download the data in search-output-exon as a tsv file
    by clicking btn-exon

    Args:
        n_clicks (int): number of clicks of button btn-exon
        data (table): table in search-output-exon

    Returns:
        dict: dict(content=df.to_csv(sep="\t", index=False),filename="data.tsv")
    """
    df = pd.DataFrame.from_records(data)
    return dict(content=df.to_csv(sep="\t", index=False),filename="data.tsv")

@app.callback(
        Output("download-transcript-id", "data"),
        [Input("btn-transcript-id", "n_clicks")],
        [State("search-output-transcript-id", "data")],
        prevent_initial_call=True
)
def download_table_transcripts_from_transcript_id(n_clicks, data):
    """Downloads the data in search-output-transcript-id as a tsv file
    by clicking btn-transcript-id

    Args:
        n_clicks (int): number of clicks of button btn-transcript-id
        data (table): table in search-output-transcript-id

    Returns:
        dict: dict(content=df.to_csv(sep="\t", index=False),filename="data.tsv")
    """
    df = pd.DataFrame.from_records(data)
    return dict(content=df.to_csv(sep="\t", index=False),filename="data.tsv")

@app.callback(
    Output("my-dynamic-dropdown", "options"),
    Input("my-dynamic-dropdown", "search_value")
)
def update_options(search_value):
    """Update the options dropdown dynamically based on the search_value.

    Args:
        search_value (string): input will change dropdown dynamically

    Raises:
        exceptions.PreventUpdate: if no input is given

    Returns:
        list: list of genes
    """

    if not search_value:
        raise exceptions.PreventUpdate

    search_value_upper = search_value.upper()
    # filter, then limit to a reasonable number to keep the UI responsive
    return [o for o in ref_gene_names if search_value_upper in o["label"]][0:max_dropdown_options-1:]


### download cDNAs
@app.callback(
   Output('download-cDNA','data'),
   [Input('mRNA-button', 'n_clicks')],
   [State('my-dynamic-dropdown', 'value')],
   prevent_initial_call=True,
)
def downloadmRNA(n_clicks, value):
    """Download a file comprising the mRNA sequences of selected transcripts.

    Args:
        n_clicks (int): number of clicks from mRNA-button
        value (string): gene name from dynamic dropdown

    Returns:
        dict: mRNA sequences in a fasta file
    """
    if (n_clicks is None) or (n_clicks == 0):
        return (no_update)
    if n_clicks is not None and n_clicks>0:
        mRNA = get_mRNA_from_gene_name(value)
        fasta = ''.join([ ">"+key+"\n"+mRNA[key]+"\n" for key in mRNA ])
        return dict(content=fasta,filename="mRNA_"+value+".fasta")

@app.callback(
   Output('download-proteins','data'),
   [Input('protein-button', 'n_clicks')],
   State('my-dynamic-dropdown', 'value'),
   prevent_initial_call=True,
)
def downloadProteins(n_clicks,value):
    """  Download the proteins without using tmp file

    Args:
        n_clicks (int): number of clicks from protein-button
        value (string): gene name from dynamic dropdown

    Returns:
        dict: proteins in a fasta file
    """
    print("----download proteins----")
    if (n_clicks is None) or (n_clicks == 0):
        return (no_update)
    if n_clicks is not None and n_clicks>0:
        proteins = get_proteins_by_gene(value)
        fasta = ''.join([ ">"+key+"\n"+str(proteins[key])+"\n" for key in proteins ])
        return dict(content=fasta,filename="proteins_"+value+".fasta")

@app.callback(
   Output('download-domains','data'),
   [Input('domains-button', 'n_clicks')],
   State('my-dynamic-dropdown', 'value'),
   prevent_initial_call=True,
)
def downloadDomains(n_clicks,value):
    """ Run InterProScan for all transcripts in proteins.fasta
    And return the file domains.fasta

    Args:
        n_clicks (int): number of clicks from domains-button
        value (string): genes name from dynamic dropdown

    Returns:
        dict: protein_domains as tsv file
    """
    if (n_clicks is None) or (n_clicks == 0):
        return (no_update)
    else:
        res = con.execute( """select t.transcript_id, d.domain, d.start, d.end,
                                     (select sum(min(e.end,t.cds_start) - min(e.start,t.cds_start))
                                      from exons e where e.transcript=t.id),
                                     t.cds_start, t.cds_end
                              from transcripts t, domains d
                              where t.id=d.transcript and t.gene_name=?""", (value,) )

        tsv = "transcript\tdomain\tstart\tend\n"

        # At this point, we have coordinates for the domains in terms of
        # offsets on the mRNA sequence, but on the forward strand of the
        # genome.  To get (back) to protein coordinates, we need to subtract
        # the position of the start codon in the same coordinate system, then
        # divide by three.  And then flip it around on the reverse strand.

        for tid, dom, start, end, cds_start, gs, ge in res:
            if gs<=ge:
                eff_start = start-cds_start
                eff_end = end-cds_start
            else:
                eff_start = cds_start-end
                eff_end = cds_start-start

            tsv += "%s\t%s\t%d\t%d\n" % (tid, dom, eff_start/3, eff_end/3)

        return dict(content=tsv,filename="domains_"+value+".tsv")

@app.callback(
    [Output('search-output-gene-id','data'),
     Output('search-output-gene-id','columns'),
     Output('gene-id-err','children')],
     Input('gene-id', 'value')
)
def get_transcripts_from_gene(input_value):
    """ Callback function to get all transcripts from a gene name

    Args:
        input_value (string): gene name from input field

    Raises:
        exceptions.PreventUpdate: _description_

    Returns:
        (data, columns, ""): returns the data and columns for the table and an empty string for the error message
    """
    statement_get_all_transcripts = "SELECT gene_name,transcript_id FROM transcripts WHERE gene_name=?"
    if input_value is None:
        raise exceptions.PreventUpdate
    print("----Get transcripts from gene_name----")
    res = con.execute(statement_get_all_transcripts, [input_value])
    df = pd.DataFrame(res.fetchall())
    if df.empty:
        return (no_update, no_update, "Not available")
    else:
        df.columns = ["gene_name","transcript_id"]
        columns = [{"name": i, "id": i} for i in df.columns]
        data = df.to_dict('records')
        return (data, columns, '')

@app.callback(
    [Output('search-output-transcript-id','data'),
     Output('search-output-transcript-id','columns'),
     Output('transcript-id-err','children')],
     Input('transcript-id', 'value')
)
def get_exon_infomration_from_transcript_id(input_value):
    """  Callback function to get exon informations from a transcript_id

    Args:
        input_value (string): gene name from input field

    Raises:
        exceptions.PreventUpdate: _description_

    Returns:
        (data, columns, ''): data and columns for the table and an empty string for the error message
    """
    statement_get_exon_structure = "SELECT e.sequence_number, e.chrom, e.start, e.end, e.strand from transcripts as t, exons as e WHERE e.transcript==t.id and t.transcript_id==?"

    if input_value is None:
        raise exceptions.PreventUpdate
    print("----Get exon structure from transcript_id----")

    res = con.execute(statement_get_exon_structure, [input_value])
    df = pd.DataFrame(res.fetchall())
    if df.empty:
        return (no_update, no_update, "Not available")
    else:
        columns = ["exon number", "chromosom", "start", "end", "strand"]
        df.columns = columns
        columns = [{"name": i, "id": i} for i in df.columns]
        data = df.to_dict('records')
        return (data, columns, '')

#### search by exon chr:start:end
@app.callback(
   [Output('search-output-exon','data'),
    Output('search-output-exon','columns'),
    Output('output-exon-err','children')],
    Input('my-button', 'n_clicks'),
    [State('text-chro', 'value'),
     State('text-start','value'),
     State('text-stop','value'),
     State('text-strand','value')]
)
def updateExonSearch(n_clicks, chro, start, stop, strand):
    """  Callback function to get all transcripts that have exon which lays within or overlaps genomic region


    Args:
        n_clicks (int): number of clicks from exon-search-button
        chro (string): chromosome
        start (string): start position on genome
        stop (string): stop position on genome
        strand (string): strand where the exons lays on

    Returns:
        (data, columns, ''): data and columns for the table and an empty string for the error message
    """
    statement = "SELECT DISTINCT t.gene_name, t.transcript_id FROM transcripts as t, exons as e WHERE t.id=e.transcript AND e.start>=? AND e.end<=? AND e.start<? AND e.chrom=? AND e.strand=?"
    if (n_clicks is None) or (n_clicks == 0):
        return (no_update, no_update, '')
    if n_clicks is not None:
        if n_clicks>0:
            print("----Search for transcripts in genomic region----")
            res = con.execute(statement,(start,stop,stop, chro,strand))
            df = pd.DataFrame(res.fetchall())
            if df.empty:
                return (no_update, no_update, "Not available")
            else:
                columns = ["gene_name","transcript_id"]
                df.columns = columns
                columns = [{"name": i, "id":i } for i in df.columns]
                data = df.to_dict('records')
                return (data, columns, '')
        else:
            return (no_update, no_update, '')

@app.callback(
        [Output('table-sql-output','data'),
        Output('table-sql-output','columns'),
        Output('sql-warning','children')],
        [Input('sql-button','n_clicks')],
        State('sql-input', 'value'),
        prevent_initial_call=True
)
def run_sql_statement(n_clicks, value):
    """Take SQL statement from input and run it

    Args:
        n_clicks (int): button click
        value (string): SQL statement

    Returns:
        (list, list, string): data, columns, warning
    """
    if n_clicks != 0:
        print("----Run SQL statement----")
        print(value)
        #if value includes "drop" or "delete" return error
        list = ["drop", "delete", "insert","select into", "update","create","alter","grant","revoke","attach","detach","vacuum","pragma"]
        # if any of the words in list is in value.lower() return error
        if any(x in value.lower() for x in list):
            return (no_update, no_update, "SQL statement not valid: drop/delete or update/insert not allowed")

        cur = con.cursor()
        #try execute sql statement else return error
        try:
            res = cur.execute(value)
        except sqlite3.OperationalError as e:
            return (no_update, no_update, "SQL statement not valid:"+ str(e))
        #if sql statement is valid, return data
        df = pd.DataFrame(res.fetchall())
        if df.empty:
            #return empty data table
            cur.close()
            return  ([], [], "No data found")
        else:
            print(cur.description)
            df.columns = [i[0] for i in cur.description]
            columns = [{"name": i, "id": i} for i in df.columns]
            data = df.to_dict('records')
            cur.close()
            return (data, columns, "")


# ! Main and most important callback in the app
@app.callback(
   [Output('search-output-ref-geneA','data'),
    Output('search-output-ref-geneA','columns'),
    Output('gene-png','src'),
    Output('no-groups-warning','children'),
    Output('start','value'),
    Output('stop','value'),
    Output('chrom','value'),
    Output('strand','value'),
    Output('mutation','value'),
    Output('heatmap-relatives','figure'),
    Output('heatmap-relatives-div','style'),
    Output('heatmap-absolutes','figure'),
    Output('heatmap-absolutes-div','style')],
    [Input('transcript-button', 'n_clicks'),
     Input('update-button', 'n_clicks')],
    [State('my-dynamic-dropdown', 'value'),
    State('groupA', 'children'),
    State('groupB','children'),
    State('start','value'),
    State('stop','value'),
    State('chrom','value'),
    State('strand','value'),
    State('mutation','value'),
    State('all-tissues-dropdown','value')],
    prevent_initial_call=True
)
def get_transcripts_from_gene(transcript_button_clicks, update_button_clicks, input_value, groupA, groupB, start, stop, chrom, strand, mutation, tissue_dropdown):
    """Get all transcripts from gene name, then update svg.
    If no groups are selected: without TPMs
    If groups are selected: with TPMs
    If all-tissues-dropdown is not empty: with TPMs for all tissues

    If Transcript-button is clicked: get all transcripts from gene name and calculate new SVG
    If Update-button is clicked: insert mutation or new start, end postition and update SVG

    Args:
        transcript_button_clicks (int): number of clicks on transcript-button
        update_button_clicks (int): number of clicks on update-button
        input_value (string): gene name from dynamic dropdown
        groupA (list): list of samples or sra ids for groupA getting from GroupA children component
        groupB (list): list of samples or sra ids for groupB getting from GroupB children component
        start (string): start position on genome
        stop (string): stop position on genome
        chrom (string): chromosome
        strand (string): strand
        mutation (string): position of mutation in genome
        tissue_dropdown (list): list of tissues from all-tissues-dropdown

    Raises:
        exceptions.PreventUpdate: _description_
        exceptions.PreventUpdate: _description_
        exceptions.PreventUpdate: _description_

    Returns:
        data, columns, svg, warning, start, stop, chrom, strand, mutation, heatmap_relatives, heatmap_absolutes
    """
    # either dropdown or start, stop, chrom als input
    triggered_id = ctx.triggered_id #check if update-button or search-button is triggered

    if triggered_id == 'transcript-button': #if search-search button is triggered
        #split inpute value for _
        if input_value is None: #if input gene is none
            print("--input gene None--")
            raise exceptions.PreventUpdate
        elif tissue_dropdown is not None:
            print("--across tissues--")
            result_df = get_TPM_from_tissues(input_value, tissue_dropdown)

            #generate heatmap with all tissues and transcripts
            heatmap_rel = generate_heatmap(result_df, True)
            heatmap_abs = generate_heatmap(result_df, False)
            print("Average TPM calcuated")
            start, stop, chrom, strand, drawing = generate_svg(result_df["transcript_id"], mutation)

            #add nhref for each transcript_id in NCBI
            columns, data = transform_to_columns_and_data(result_df)

            return (data, columns, b64_svg(drawing),'', start, stop, chrom, strand, None, heatmap_rel, {'display': 'block'}, heatmap_abs, {'display': 'block'})

        else:
            #if user did not select group A or B then present only transcripts in search-output-ref-geneA
            if ((groupA or groupB) == "none") or ((groupA or groupB) == []) or ((groupA or groupB) == None) or ((groupA or groupB) == "") or ((groupA or groupB) == "[]") or ((groupA or groupB) == "None"):
                print("--no groups selected--")
                # input_value = input_value.split(" ")[0]
                #get all transcripts with gene id and gene name over ref gene id
                statement_get_all_transcripts_from_gene = "SELECT DISTINCT gene_name, transcript_id FROM transcripts WHERE gene_name=? ORDER BY transcript_id"
                res = con.execute(statement_get_all_transcripts_from_gene, [input_value])
                df = pd.DataFrame(res.fetchall())
                transcripts = df[1]
                df.columns = ["gene_name","transcript_id"]
                #sort df by transcript_id
                columns, data = transform_to_columns_and_data(df, False)

                start, stop, chrom, strand, drawing = generate_svg(transcripts,mutation)


                return (data, columns, b64_svg(drawing), "No groups selected! The transcripts \
                        shown here are possible variants resulting from all analyzed datasets (a list of all datasets see 'info').\
                         This means the user sees a global, unified set of transcripts across multiple RNA-Seq samples.\
                         To compare the quantification (TPM) of transcript-variants per tissue or sample, \
                        select the corresponding groups.", start, stop, chrom, strand, None, None, {'display': 'none'}, None, {'display': 'none'})

            else: #if user did select group then present calculate avg tpm and fpkm from groupA and groupB
                print("--groups selected--")
                # input_value = input_value.split(" ")[0]
                groupA = groupA.split(",")
                groupB = groupB.split(",")
                #for each element in groupA and groupB remove "[" and "]" and "'"
                groupA = [x.strip('[').strip(']').strip("'").strip(' ') for x in groupA]
                groupB = [x.strip('[').strip(']').strip("'").strip(' ') for x in groupB]
                print(groupA)
                print(groupB)

                merged_df = get_group_comparisons_from_gene(input_value, groupA, groupB)
                df_result = calculate_percentage_for_TPM(merged_df)
                df_result = df_result.sort_values(by=['transcript_id'])
                heatmap_rel = generate_heatmap(df_result, True)
                heatmap_abs = generate_heatmap(df_result, False)

                columns, data = transform_to_columns_and_data(df_result)

                print("Average TPM calcuated")
                start, stop, chrom, strand, drawing = generate_svg(df_result["transcript_id"],mutation)
                return (data, columns, b64_svg(drawing),'', start, stop, chrom, strand, None, heatmap_rel, {'display': 'block'}, heatmap_abs, {'display': 'block'})


    elif triggered_id == 'update-button': #if update-button is triggered
        # if tissue_dropdown is None:
        #     print("no tissues selected but update view")

        if tissue_dropdown is not None:
            print("--across tissues--")

            # get all transcripts by genomic position and get TPM values from them
            transcripts = get_transcripts_by_gennomics_pos(chrom, start, stop, strand)
            print(transcripts)
            if transcripts.empty:
                return ([], [], None, "No transcript found for genomic region", start, stop, chrom, strand, mutation, None, {'display': 'none'}, None, {'display': 'none'})
            df_result = get_TPM_from_tissues_over_transcripts(transcripts["transcript_id"], tissue_dropdown)
            # generate heatmap with all tissues and transcripts
            heatmap_rel = generate_heatmap(df_result, True)
            heatmap_abs = generate_heatmap(df_result, False)

            columns, data = transform_to_columns_and_data(df_result)

            print("Average TPM calcuated")
            start, stop, chrom, _strand, drawing = generate_svg(df_result["transcript_id"],mutation)
            return (data, columns, b64_svg(drawing),'', start, stop, chrom, strand, mutation, heatmap_rel, {'display': 'block'}, heatmap_abs, {'display': 'block'})

        elif ((groupA or groupB) == "none") or ((groupA or groupB) == []) or ((groupA or groupB) == None):
            print("--no groups selected--")
            df_result = get_transcripts_by_gennomics_pos(chrom, start, stop, strand)
            if df_result.empty:
                return ([], [], None, "No transcript found for genomic region", start, stop, chrom, strand, mutation, None,{'display': 'none'}, None, {'display': 'none'}, None, None)
            #sort df_result by transcript_id
            df_result = df_result.sort_values(by=['transcript_id'])
            columns, data = transform_to_columns_and_data(df_result)


            columns = [{"name": i, "id":i } for i in df_result.columns]
            data = df_result.to_dict('records')
            start_result, stop_results, chrom_results, strand_results, drawing = generate_svg(df_result["transcript_id"], mutation)
            return (data, columns, b64_svg(drawing), "Updated", start, stop, chrom, strand, mutation, None, {'display': 'none'}, None, {'display': 'none'})
        else:
            print("--groups selected--")
            transcripts = get_transcripts_by_gennomics_pos(chrom, start, stop, strand)
            if transcripts.empty:
                return ([], [], None, "No transcript found for genomic region", start, stop, chrom, strand, mutation, None,{'display': 'none'}, None, {'display': 'none'})
            print(transcripts)
            df_result = get_group_comparisons_over_transcripts([t for t in transcripts["transcript_id"]], groupA, groupB)
            df_result = calculate_percentage_for_TPM(df_result)
            df_result = df_result.sort_values(by=['transcript_id'])
            columns, data = transform_to_columns_and_data(df_result)

            heatmap_rel = generate_heatmap(df_result, True)
            heatmap_abs = generate_heatmap(df_result, False)
            start_result, stop_results, chrom_results, strand_results, drawing = generate_svg(df_result["transcript_id"], mutation)
            return (data, columns, b64_svg(drawing), "Updated", start, stop, chrom, strand, mutation, None, {'display': 'none'}, None, {'display': 'none'})

# start the server of the  app
server = app.server

# run the app
if __name__ == '__main__':
    app.run_server(debug=True)

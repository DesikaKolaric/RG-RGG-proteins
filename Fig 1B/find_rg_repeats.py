import os,sys
from Bio import SeqIO,Seq,SeqFeature
import os,sys
import re
import numpy as np
import matplotlib.pyplot as plt
import copy
import pandas as pd
import re

"""
Given the fasta files with canonical sequences from the human proteome
    >> extract IDR sequences
    >> find RG repeats by iteratively increasing the number of RG motifs in the sequence 
    from 1 to 12 and count the number of such RG regions per each IDR 
"""

# reads in multi-entry fasta file using SeqIO methods
#@param: fasta_file - path to protein seq file in fasta format
#@retval: record - SeqIO sequence data storage structure
def read_fasta(fasta_file):
    if not os.path.isfile(fasta_file):
        raise Exception("FASTA file path is incorrect or the file does not exist! Exiting...")
        sys.exit(1)
    fasta_data = list(SeqIO.parse(fasta_file, "fasta"))
    return fasta_data

#reads SeqIO sequence data storage structure and creates a dictionary
#@param: fasta_data - SeqIO data storage variable
#@retval: fasta dictionary - UniprotIDs_IDRNumber as keys, sequences as values
def make_fasta_dict(fasta_data):
    fasta_dict={}
    for i in range(len(fasta_data)):
        fasta_id=str(fasta_data[i].id.split("|")[1]) + '_' + str(fasta_data[i].id.split(":")[1])
        sequence=fasta_data[i].seq
    return fasta_dict

if __name__=="__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <fpath> <rg>")
        sys.exit(1)
    fpath = sys.argv[1]   # path\to\fasta\file\UP000005640_9606_SPOTD_MIN_30AA.fasta
    rg = sys.argv[2]
    if rg not in ['RG', 'RGG']:
        print("Invalid value for rg. Allowed values are 'RG' or 'RGG'.")
        sys.exit(1)

    dfasta=read_fasta(fpath) # returns SeqOD the fasta file data storage
    dfasta=make_fasta_dict(dfasta) # returns fasta dictionary keys: UniprotIds_IDRNumber, values: protein sequence
    N_IDR = len(dfasta)
    final_matrix = np.zeros((13,N_IDR), dtype = object) #output matrix with UniprotIDs_IDRNumber as columns and the number of occurrences of each RG motif (1-12 RG) in each of them
    check_cnt = 0
    cnt_IDR = -1
    # in each sequence from fasta file look for the number of RG.{0.5}RG repeats and iteratively increase by one RG.{0,5} to max 12 RG.{0,5}
    for key in dfasta:
        cnt_IDR += 1
        query_seq = dfasta[key]
        final_matrix[0][cnt_IDR] = key
        if rg == 'RG':
            string_rg='(?=(RG.{0,5}RG'  # the first RG motif we look for in the sequences using Regex - later we iteratively add .{0,5}RG to it
            add_rg = '.{0,5}RG'
        if rg == 'RGG':
            string_rg='(?=(RGG.{0,5}RGG'  # the first RG motif we look for in the sequences using Regex - later we iteratively add .{0,5}RG to it
            add_rg = '.{0,5}RGG'
        for tg_number in range(1,13):
            final_rg = first_rg + '))'
            s = str(query_seq)
            next_rg = string_rg + '))'
            first_rg = re.compile(next_rg)
            matches = re.finditer(first_rg, s)
            results = [(match.start(0)+1,int(match.start(0))+len(match.group(1))) for match in matches]
            final_matrix[tg_number][cnt_IDR] = len(results)
            string_rg = string_rg + add_rg

    fn_mtx = final_matrix.transpose()
    df = pd.DataFrame(fn_mtx)
    df.columns = df.iloc[0]
    df = df[1:]
    if rg == 'RG':
        df.to_csv('occurence_matrix_RG.txt', sep='\t', encoding='utf-8', header=True, index = None)
    if rg == 'RGG':
        df.to_csv('occurence_matrix_RGG.txt', sep='\t', encoding='utf-8', header=True, index = None)
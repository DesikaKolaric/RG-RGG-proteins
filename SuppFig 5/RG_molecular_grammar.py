import os,sys
import re
from Bio import SeqIO,Seq,SeqFeature
import os,sys
import re
import numpy as np
import matplotlib.pyplot as plt
import copy
import pandas as pd

"""
Given the IDR bounds and the fasta files with canonical sequences from the human proteome
    >> extract IDR sequences
    >>
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
def make_fasta_dict(fasta_data, infile):
    fasta_dict={}
    for i in range(len(fasta_data)):
        if infile == 'UP000005640_9606_SPOTD_MIN_30AA.fasta':
            fasta_id=str(fasta_data[i].id.split("|")[1]) + '_' + str(fasta_data[i].id.split(":")[1])
        else:
            fasta_id = str(fasta_data[i].id)
        sequence=fasta_data[i].seq
        fasta_dict.setdefault(fasta_id,sequence)
    return fasta_dict

# the ptah to sequences of human IDRome, the sequences of proteins that localize to stress granules and that bind to TNPO1 need to be given as input separately
if __name__=="__main__":
    fpath=sys.argv[1] #path to the fasta file
    filename = fpath.split('\\')[-1]
    dfasta=read_fasta(fpath) # all the fasta file data storage
    dfasta=make_fasta_dict(dfasta, filename) # fasta dictionary keys: uniprot ids, values: protein sequence
    string_rg='(?=(RG.{0,5}RG.{0,5}RG'
    add_rg = '.{0,5}RG'
    global_rg_cnt = 0
    global_dict = {'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0, 'I':0, 'K':0, 'L':0, 'M':0, 'N':0, 'P':0, 'Q':0,
                    'R':0, 'S':0, 'T':0, 'V':0, 'W':0, 'Y':0}
    check_cnt = 0
    for rg in range(3,6):
        next_rg = string_rg + '))'
        for key in dfasta:
            query_seq = dfasta[key]
            s = str(query_seq)
            first_rg = re.compile(next_rg)
            matches = re.finditer(first_rg, s)
            results = [(match.start(0)+1,int(match.start(0))+len(match.group(1))-1) for match in matches]
            if (len(results) != 0):
                used_indices = []
                for res in range(len(results)):
                    indices = results[res]
                    if indices[0] >=11:
                        strt = indices[0] - 11
                    else:
                        strt = 0
                    if (indices[1]+10) <= (len(s)-1):
                        ed = indices[1] + 10
                    else:
                        ed = len(s)-1
                    for aa in range(strt,ed+1):
                        if aa not in used_indices:
                            used_indices.append(aa)
                            residue = query_seq[aa]
                            global_dict[residue] += 1
                            global_rg_cnt += 1
        global_occ = open(f'RG_{rg}_{filename}_composition.txt', "w")
        for aa_name in global_dict.keys():
            global_dict[aa_name] = global_dict[aa_name]/(global_rg_cnt-global_dict[aa_name])
            global_occ.write('%s\t%.5f\n'%(aa_name, global_dict[aa_name]))
        global_occ.close()
        string_rg = string_rg + add_rg
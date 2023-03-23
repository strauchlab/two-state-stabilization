#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 12:34:51 2019

@author: Karen J. Gonzalez
"""

import os
import pandas as pd
import subprocess
import argparse

def contractSeq(sequence): #sequence has to be a list with each aa as an independent element
    seqExtended = []
    for i in sequence:
        if i == "ARG":
            seqExtended.append("R")
        elif i == "HIS":
            seqExtended.append("H")    
        elif i == "LYS":
            seqExtended.append("K")            
        elif i == "ASP":
            seqExtended.append("D")            
        elif i == "GLU":
            seqExtended.append("E")            
        elif i == "SER":
            seqExtended.append("S")            
        elif i == "THR":
            seqExtended.append("T")            
        elif i == "ASN":
            seqExtended.append("N")    
        elif i == "GLN":
            seqExtended.append("Q")    
        elif i == "CYS":
            seqExtended.append("C")    
        elif i == "GLY":
            seqExtended.append("G")            
        elif i == "PRO":
            seqExtended.append("P")            
        elif i == "ALA":
            seqExtended.append("A")            
        elif i == "VAL":
            seqExtended.append("V")            
        elif i == "ILE":
            seqExtended.append("I")            
        elif i == "LEU":
            seqExtended.append("L")    
        elif i == "MET":
            seqExtended.append("M")    
        elif i == "PHE":
            seqExtended.append("F")    
        elif i == "TYR":
            seqExtended.append("Y")
        elif i == "TRP":
            seqExtended.append("W") 
    return seqExtended 

def resfile_postfusion(out_dir, positions_post, positions_pre, seq, positions_to_exclude = [], name_resfile="resfile"):
    os.chdir(out_dir)
    resfile = open(name_resfile,'w')
    resfile.write("nataa\nstart\n") 
    counter = -1    
    for i in range(len(positions_pre)):
        if positions_pre[i] in positions_to_exclude:
            continue
        else:
            counter+=1
            resfile.write(positions_post[counter][0]+'\t'+positions_post[counter][1]+'\tPIKAA\t')
            resfile.write(seq[i])
            resfile.write("\n")           
    resfile.close()  

def unique_sequences(list_folders, results_dir, hotspots):               
    seqs_hotspots = {}
    repeated_seqs = {}       
    for folder in list_folders:
        new_dir = os.path.join(results_dir, folder)
        files = [f for f in os.listdir(new_dir) if f.endswith(".pdb")]
        position_old = []        
        for f in range(len(files)):
            name = os.path.join(new_dir, files[f])
            sequence = ''
            with open(name,"r") as fopen:
                for line in fopen:
                    if line.startswith("ATOM"):
                        new_position = "".join([line[22:27].strip(),line[21]])                       
                        if new_position in hotspots:
                            if position_old != new_position:
                                aa = contractSeq([line[17:20]])
                                sequence+=aa[0]
                                position_old = new_position                                
            #check if the sequence is identical to another one in the set
            keys = [key  for (key, value) in seqs_hotspots.items() if value == sequence]
            if keys == []: #if the key is empty, that means that the sequence is different from the other ones and we need to add this new seq
                seqs_hotspots[name] = sequence
            else: #if the seq already exist in the dictionary, then add this repeated seq to the repeated seq dictionary, being the key the first pdb that appears with that sequence
                if keys[0] in repeated_seqs.keys():
                    repeated_seqs[keys[0]].append(name)
                else:
                    repeated_seqs[keys[0]]= [name]
    return seqs_hotspots, repeated_seqs
   
def read_alignment(alignment):  #alignment in fasta format. Make sure the sequence alignment corresponds with the structural alignment
    with open(alignment, "r") as fopen:
        sequences = {}
        for line in fopen:
            line = line.strip()
            if line.startswith(">"):
                seq_name=line[1:]
                sequences[seq_name] = []
            else:
                for aa in line:
                    sequences[seq_name].append(aa)
    sequences = pd.DataFrame.from_dict(sequences)
    return sequences

def add_data_according_to_alignment(protein, alignment_dataframe, data_values, column_name="PDB_#" ):
    new = []
    counter=0
    for index, row in alignment_dataframe.iterrows():
        if alignment_dataframe.loc[index,protein] == "-":
            new.append("-") 
        else:
            new.append(data_values[counter])
            counter+=1                
    new_column = protein + "_"+column_name       
    alignment_dataframe[new_column] = new
    return alignment_dataframe

def sequence_from_pdb(pdb, chain="all", contracted=True):
    positions = []
    seq = []
    chains =[] 
    fopen = open(pdb,"r")
    old = ''
    for line in fopen:
        if line.startswith("ATOM"):
            new = line[22:27].strip()
            if new != old:
                ch = line[21] #chain in the pdb
                if chain =="all":
                    chains.append(ch)
                else:
                    if not ch in chain:
                        continue
                    else:
                        chains.append(ch)
                aa = line[17:20]
                seq.append(aa)
                positions.append(new)
                old = new
    if contracted:
        seq = contractSeq(seq)
    return seq, positions, chains

def gather_total_score(pdbs):
    total_score = []
    idxes = []
    for pdb in pdbs:
        with open(pdb,'r') as fopen:
            for line in fopen:
                if line.startswith("score_MUT_total_energy"):
                    line = line.split()
                    total_score.append([float(line[1]), pdb])
                    idx = pdb.split("/")
                    idx = idx[-2]+"_"+idx[-1][-8:-4]
                    idxes.append(idx) 
                    break
    total_score = pd.DataFrame(total_score, columns=["total_score","description"], index=idxes)    
    return total_score

############################################################################### 
root = os.getcwd()   
### Command line arguments 
parser = argparse.ArgumentParser(description='Identification of designs improving the prefusion energy')
parser.add_argument('arg_file', help = "File containing all input arguments")
args = parser.parse_args()
# Input arguments from arg_file
arg_dict = {}
with open(args.arg_file, "r") as fopen:
    for line in fopen:
        if line and line[0].isalpha():
            line = line.split("=")
            arg_dict[line[0].strip()] = line[1].strip()    

### Identify positions where substitutions were made
hotspots = []
with open(arg_dict["resfile_pre"], 'r') as fopen:
    lines = fopen.readlines()[2:]
    for l in lines:
        position = "".join(l.split()[:2])
        hotspots.append(position)  
 
### Gather unique sequences from prefusion designs                  
results_dir = arg_dict["results_dir"]
list_folders = arg_dict["results_folder"].split(",")
seqs_hotspots, repeated_seqs = unique_sequences(list_folders, arg_dict["results_dir"], hotspots)

### Compare total energy of designed sequences with prefusion control 
control_energy = gather_total_score([arg_dict["control_dir"]])
designs_energy = gather_total_score(seqs_hotspots.keys())
# Select sequences where the prefusion design improved energy compared to control
candidates = designs_energy[designs_energy['total_score']<= control_energy["total_score"].tolist()[0]]

## Transfer candidate sequences to a separate folder
selection_path = os.path.join(results_dir, "selection_pre_E")
subprocess.call(["mkdir", selection_path])  
for index,row in candidates.iterrows():
    new_name = os.path.join(selection_path,index+".pdb")
    subprocess.call(["cp", candidates.loc[index,"description"],new_name])  
    
## Make resfiles for postfusion design based on sequences improving prefusion energy
# Open sequence alignment to identify correspondance between residue position
seq_pre, pdb_positions_pre, ch_pre = sequence_from_pdb(arg_dict["pdb_pre"], chain=arg_dict["pre_monomer_ch"]) 
seq_post, pdb_positions_post, ch_post = sequence_from_pdb(arg_dict["pdb_post"], chain=arg_dict["post_monomer_ch"]) 

alignment = read_alignment(arg_dict["alignment"])
alignment = alignment.rename(columns={arg_dict["pre_name_alignment"]: "prefusion", arg_dict["post_name_alignment"]: "postfusion"})
alignment = add_data_according_to_alignment("postfusion", alignment, pdb_positions_post, column_name="PDB_#" )
alignment = add_data_according_to_alignment("prefusion", alignment, pdb_positions_pre, column_name="PDB_#" )
alignment = add_data_according_to_alignment("postfusion", alignment, ch_post, column_name="chain")
alignment = add_data_according_to_alignment("prefusion", alignment, ch_pre, column_name="chain")

positions_post = []
positions_to_exclude = [] # for when the residue is not present in the postfusion structure (e.g. fusion peptide)
control_seq = []
for residue in hotspots:
    chain = residue[-1]
    position = int(residue[:-1])
    row = alignment[(alignment['prefusion_PDB_#'] == position) & (alignment['prefusion_chain'] == chain)]
    if row["postfusion_PDB_#"].tolist()[0] != "-":
           positions_post.append([str(row["postfusion_PDB_#"].tolist()[0]), row["postfusion_chain"].tolist()[0]])  
           control_seq.append(row["postfusion"].tolist()[0])                              
    else:
        positions_to_exclude.append(residue)
        
## Make postfusion resfiles    
for index, row in candidates.iterrows():
    output_dir = os.path.join(root,"postfusion_designs")
    subprocess.call(["mkdir", output_dir])
    folder = os.path.join(output_dir, index)
    subprocess.call(["mkdir", folder])
    resfile_postfusion(folder, positions_post, hotspots, seqs_hotspots[candidates.loc[index,"description"]], positions_to_exclude = positions_to_exclude, name_resfile="resfile_post")
    
# Make control postfusion resfile
control_out = os.path.join(output_dir, "control")
subprocess.call(["mkdir", control_out])  
resfile_postfusion(control_out, positions_post, positions_post, control_seq,  name_resfile="resfile_post_control")




























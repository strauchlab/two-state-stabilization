#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: karen
"""
import os
import pandas as pd
import subprocess
import argparse
from pyrosetta import *
from rosetta import *
from rosetta.core.scoring import *


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
                    idx = pdb.split("/")[-1][:-4]
                    idxes.append(idx) 
                    break
    total_score = pd.DataFrame(total_score, columns=["total_score","description"], index=idxes)    
    return total_score

def mutated_residues(seq_wt, seq_mut):
    positions = []
    mutations = []    
    for i in range(len(seq_wt)):
        if seq_wt[i] != seq_mut[i]:
            mut = [seq_wt[i], i+1 ,seq_mut[i]]
            mutations.append(mut)
            positions.append(i+1)
    return positions, mutations

def energy_terms(pose, weight_scores):  
    energies_array = pose.energies().residue_total_energies_array() #With exception of total_score, all score terms are unweighted energies
    columns = list(energies_array.dtype.names)
    unweighted_data = energies_array.tolist()
    unweighted_data = pd.DataFrame(unweighted_data, columns= columns, index=list(range(1,len(unweighted_data)+1))) #rosetta_numbering
   
    weighted_data = unweighted_data.copy(deep=True)
    for key, value in weighted_data.iteritems():
        if key in weight_scores:
            weighted_data[key] = weighted_data[key] * weight_scores[key]
    return weighted_data, unweighted_data


############################################################################### 
root = os.getcwd()   
### Command line arguments 
parser = argparse.ArgumentParser(description='Identification of designs improving prefusion over postfusion energy')
parser.add_argument('arg_file', help = "File containing all input arguments")
args = parser.parse_args()
# Input arguments from arg_file
arg_dict = {}
with open(args.arg_file, "r") as fopen:
    for line in fopen:
        if line and line[0].isalpha():
            line = line.split("=")
            arg_dict[line[0].strip()] = line[1].strip()    

### Gather prefusion designs energies               
pre_results_dir = arg_dict["pre_results_dir"]
prefusion_pdbs = [os.path.join(pre_results_dir,i) for i in os.listdir(pre_results_dir) if i.endswith("pdb")]
pre_energies = gather_total_score(prefusion_pdbs)
pre_control_energy = gather_total_score([arg_dict["pre_control_pdb"]])
pre_control_energy.index = ["control"]
pre_energies = pd.concat([pre_control_energy, pre_energies])

### Gather postfusion designs energies
post_results_dir = arg_dict["post_results_dir"]
# Transfer postfusion designs to one folder
post_folders = [i for i in os.listdir(post_results_dir) if os.path.isdir(post_results_dir+"/"+i)]
post_out_dir = os.path.join(post_results_dir, "selection_pre_E")
subprocess.call(["mkdir", post_out_dir])
for folder in post_folders:
    path = os.path.join(post_results_dir,folder)
    pdb = [i for i in os.listdir(path) if i.endswith("1.pdb")][0]
    old_name = os.path.join(path,pdb)
    new_name = path.split("/")[-1]+".pdb"
    new_name = os.path.join(post_out_dir,new_name)
    subprocess.call(["cp", old_name, new_name])  
# Gather energies
postfusion_pdbs = [os.path.join(post_out_dir,i) for i in os.listdir(post_out_dir) if i.endswith("pdb")]
post_energies = gather_total_score(postfusion_pdbs)

### Select sequences where the postfusion design increased energy compared to control
candidates = post_energies[post_energies['total_score']> post_energies.loc["control","total_score"]]

### Write output file with mutations present on each design
# Open sequence alignment to identify correspondance between residue position
pre_pdb_wt = pre_control_energy.loc["control", "description"]
post_pdb_wt = post_energies.loc["control","description"]
control_seq_pre, pdb_positions_pre, ch_pre = sequence_from_pdb(pre_pdb_wt, chain=arg_dict["pre_monomer_ch"]) 
control_seq_post, pdb_positions_post, ch_post = sequence_from_pdb(post_pdb_wt, chain=arg_dict["post_monomer_ch"]) 

alignment = read_alignment(arg_dict["alignment"])
alignment = alignment.rename(columns={arg_dict["pre_name_alignment"]: "prefusion", arg_dict["post_name_alignment"]: "postfusion"})
alignment = add_data_according_to_alignment("postfusion", alignment, list(range(1,len(control_seq_post)+1)), column_name="rosetta_#")
alignment = add_data_according_to_alignment("prefusion", alignment, list(range(1,len(control_seq_pre)+1)), column_name="rosetta_#")

### Analyze per-residue energies of each substitution on each design
init()
sfxn = core.scoring.get_score_function()
weight_scores = {"fa_rep":0.550, "fa_intra_rep":0.005, "pro_close":1.250, "dslf_fa13":1.250, "omega":0.400, "fa_dun":0.700, "p_aa_pp":0.600, "yhh_planarity":0.625,"rama_prepro":0.450  }

#data controls
#prefusion
pre_pose_wt = pose_from_pdb(pre_pdb_wt)
sfxn(pre_pose_wt)
energies_wt_pre, _ = energy_terms(pre_pose_wt, weight_scores)
energies_wt_pre = energies_wt_pre[:len(control_seq_pre)]
energies_wt_pre["seq"] = control_seq_pre
#postfusion
post_pose_wt = pose_from_pdb(post_pdb_wt)
sfxn(post_pose_wt)
energies_wt_post, _ = energy_terms(post_pose_wt, weight_scores)

all_energies = {}
mutations = {}
for idx in candidates.index:
    #prefusion
    pre_pdb_mut = pre_energies.loc[idx,"description"]
    pre_pose_mut = pose_from_pdb(pre_pdb_mut)
    sfxn(pre_pose_mut)
    energies_mut_pre, _ = energy_terms(pre_pose_mut, weight_scores)
    #postfusion
    post_pdb_mut = post_energies.loc[idx,"description"]
    post_pose_mut = pose_from_pdb(post_pdb_mut)
    sfxn(post_pose_mut)
    energies_mut_post, _ = energy_terms(post_pose_mut, weight_scores)
    #Identify positions with mutations (prefusion numbering)  
    seq_mut, _, _ = sequence_from_pdb(pre_pdb_mut, chain=arg_dict["pre_monomer_ch"]) 
    mutated_positions_pre, mutation = mutated_residues(control_seq_pre, seq_mut) # mutated positions are in rosetta numbering
    mutations[idx] = mutation 
    #calculate the energetic differences between wild-type and mutated residue
    for i in range(len(mutated_positions_pre)):
        difference_pre = energies_mut_pre.loc[mutated_positions_pre[i],"total_score"] - energies_wt_pre.loc[mutated_positions_pre[i],"total_score"]
        position_post = alignment[alignment["prefusion_rosetta_#"] == mutated_positions_pre[i]]["postfusion_rosetta_#"].tolist()[0]
        if position_post != "-":
            difference_post = energies_mut_post.loc[position_post,"total_score"] - energies_wt_post.loc[position_post,"total_score"]
        else:
            difference_post = "N/A"        
        native_res = mutation[i][0]+str(mutation[i][1])
        mut = mutation[i][2]
        data = [mutated_positions_pre[i],position_post ,mut,difference_pre,difference_post,idx]
        if not native_res in all_energies.keys():
            all_energies[native_res] = pd.DataFrame([data], columns=["position_pre","position_post","mutation","energy_difference_pre","energy_difference_post","pdb"])
        else:
            all_energies[native_res] = all_energies[native_res].append(pd.DataFrame([data], columns=["position_pre","position_post","mutation","energy_difference_pre","energy_difference_post","pdb"]), ignore_index=True)
    
#write output files containing energetic differences per-residue  
out_path_results = os.path.join(root, "candidates"); subprocess.call(["mkdir", out_path_results])
out_path_energies = os.path.join(out_path_results, "energies"); subprocess.call(["mkdir", out_path_energies])
out_path_sequences = os.path.join(out_path_results, "sequences"); subprocess.call(["mkdir", out_path_sequences])

for key in all_energies:
    all_energies[key].to_csv(out_path_energies+"/"+ key+".csv")
    
#write output file containing all mutations (prefusion numbering)  
description =  ["wild_type"] + list(mutations.keys())
values = list(mutations.items())
columns = [j[1] for x in values for j in x[1]]
columns = sorted(list(set(columns)))
positions_post = []
for i in columns:
    positions_post.append(alignment[alignment["prefusion_rosetta_#"] == i]["postfusion_rosetta_#"].tolist()[0])

positions_post = pd.DataFrame([positions_post],columns=columns, index=["postfusion_#"])
mutations_dataframe = pd.DataFrame(columns=columns, index=description)
for value in values:
    idx = value[0]
    for i in value[1]:
        pos = i[1]
        mut = i[2]
        wt = i[0]
        mutations_dataframe.loc[idx,pos] = mut
        mutations_dataframe.loc["wild_type",pos] = wt
        
        
mutations_dataframe.fillna(value="0",inplace=True)  
mutations_dataframe = pd.concat([positions_post,mutations_dataframe])
mutations_dataframe.to_csv(out_path_sequences+"/mutations_per_design.csv")

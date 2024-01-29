#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Karen J. Gonzalez
"""
import os
import pandas as pd
import subprocess
import argparse
from pyrosetta import *
from rosetta import *
from rosetta.core.scoring import *
from statistics import mean
import numpy as np
import json

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
            new = int(line[22:26])
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

def energy_terms(pose):  
    weight_scores = {"fa_rep":0.550, "fa_intra_rep":0.005, "pro_close":1.250, "dslf_fa13":1.250, "omega":0.400, "fa_dun":0.700, "p_aa_pp":0.600, "yhh_planarity":0.625,"rama_prepro":0.450  }
    energies_array = pose.energies().residue_total_energies_array() #With exception of total_score, all score terms are unweighted energies
    columns = list(energies_array.dtype.names)
    unweighted_data = energies_array.tolist()
    unweighted_data = pd.DataFrame(unweighted_data, columns= columns, index=list(range(1,len(unweighted_data)+1))) #rosetta_numbering
   
    weighted_data = unweighted_data.copy(deep=True)
    for key, value in weighted_data.iteritems():
        if key in weight_scores:
            weighted_data[key] = weighted_data[key] * weight_scores[key]
    return weighted_data, unweighted_data

def custom_round(value):
    return round(value, 3) if not pd.isna(value) else value

def per_resi_sidechain_hbonds_and_energy(pose, total_resi, positions=""):
    hbond_set = bindings.pose.get_hbonds(pose,True,True,True,False, False) #as an object, this has the same properties as hbonds.HBondSet()
    '''(pose, calculate_derivative=True, exclude_bb=True, exclude_bsc=True, exclude_scb=False, exclude_sc=False) '''
    sc_hbonds_per_resi = [0 for x in range(total_resi)]
    hbonds_energies = [0 for x in range(total_resi)]
    if positions == "":
        positions = list(range(1, total_resi+1))        
    for i in positions:
        trace = basic.PyTracer() #this is to  be able to put the tracer in a variable
        rosetta.basic.Tracer.set_ios_hook(trace, rosetta.basic.Tracer.get_all_channels_string(), False) #this is also to manipulate the tracer
        hbond_set.show(pose,i, True, trace)
        out = str(trace.buf()).split("#")
        numb_sc_hbonds =len(out)-2
        if numb_sc_hbonds > 0:
            sc_hbonds_per_resi[i-1]= numb_sc_hbonds  
            hbond_energy = 0
            for j in range(2, len(out)):
                energy = out[j].split()
                hbond_energy+= float(energy[-1]) #the total hbond energy will be the sum of the individual hbonds and not the average as is calculated with the energy term hbond_sc
            hbonds_energies[i-1] = hbond_energy       
    return sc_hbonds_per_resi, hbonds_energies

def per_resi_packing_score(pose, repetitions=3):
    packing = protocols.pose_metric_calculators.PackstatCalculator(repetitions,True)
    packing_per_res = json.loads(packing.get("residue_packstat", pose))  
    return packing_per_res

############################################################################### 
root = os.getcwd()   
### Command line arguments 
parser = argparse.ArgumentParser(description= 'Selection of candidate prefusion designs')
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
post_folders = [i for i in os.listdir(post_results_dir) if os.path.isdir(post_results_dir+"/"+i) and not i == "selection_pre_E" ]
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
post_control_energy = gather_total_score([arg_dict["post_control_pdb"]])
post_control_energy.index = ["control"]
post_energies = pd.concat([post_control_energy, post_energies])


### Select sequences where the postfusion design increased energy compared to control
candidates = post_energies[post_energies['total_score']> post_energies.loc["control","total_score"]]

### Write output file with mutations present on each design and per-residue energy differences
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

#data controls
#prefusion
pre_pose_wt = pose_from_pdb(pre_pdb_wt)
sfxn(pre_pose_wt)
energies_wt_pre, _ = energy_terms(pre_pose_wt)
energies_wt_pre = energies_wt_pre[:len(control_seq_pre)]
energies_wt_pre["seq"] = control_seq_pre
#postfusion
post_pose_wt = pose_from_pdb(post_pdb_wt)
sfxn(post_pose_wt)
energies_wt_post, _ = energy_terms(post_pose_wt)

all_energies = {}
mutations = {}
for idx in candidates.index:
    #prefusion
    pre_pdb_mut = pre_energies.loc[idx,"description"]
    pre_pose_mut = pose_from_pdb(pre_pdb_mut)
    sfxn(pre_pose_mut)
    energies_mut_pre, _ = energy_terms(pre_pose_mut)
    #postfusion
    post_pdb_mut = post_energies.loc[idx,"description"]
    post_pose_mut = pose_from_pdb(post_pdb_mut)
    sfxn(post_pose_mut)
    energies_mut_post, _ = energy_terms(post_pose_mut)
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
            difference_post = np.nan       
        native_res = mutation[i][0]+str(mutation[i][1])
        mut = mutation[i][2]
        data = [mutated_positions_pre[i],position_post ,mut,difference_pre,difference_post,idx]
        if not native_res in all_energies.keys():
            all_energies[native_res] = pd.DataFrame([data], columns=["position_pre","position_post","mutation","energy_difference_pre","energy_difference_post","pdb"])
        else:
            all_energies[native_res] = all_energies[native_res].append(pd.DataFrame([data], columns=["position_pre","position_post","mutation","energy_difference_pre","energy_difference_post","pdb"]), ignore_index=True)
 

### write output file containing all mutations and energy differences
out_path_results = os.path.join(root, "candidates"); subprocess.call(["mkdir", out_path_results])
#Sequences and total energy    
description =  ["control"] + list(mutations.keys())
values = list(mutations.items())
columns = [j[1] for x in values for j in x[1]]
columns = sorted(list(set(columns)))
positions_post = []
for i in columns:
    positions_post.append(alignment[alignment["prefusion_rosetta_#"] == i]["postfusion_rosetta_#"].tolist()[0])

columns.extend(["raw_score_prefusion","raw_score_postfusion"])
positions_post.extend(["-","-"])
positions_post = pd.DataFrame([positions_post], columns=columns, index=["postfusion_rosetta_#"])
mutations_dataframe = pd.DataFrame(columns=columns, index=description)

for value in values:
    idx = value[0]
    for i in value[1]:
        pos = i[1]
        mut = i[2]
        wt = i[0]
        mutations_dataframe.loc[idx,pos] = mut
        mutations_dataframe.loc["control",pos] = wt
        
mutations_dataframe["raw_score_prefusion"] =  pre_energies.loc[mutations_dataframe.index,"total_score"]
mutations_dataframe["raw_score_postfusion"] =  post_energies.loc[mutations_dataframe.index,"total_score"]            
mutations_dataframe.fillna(value="-",inplace=True)  
mutations_dataframe = pd.concat([positions_post,mutations_dataframe])

#per-residue energy differences
per_res_energy_pre = pd.DataFrame(index=mutations_dataframe.iloc[2:,:].index, columns=mutations_dataframe.columns)
per_res_energy_post = pd.DataFrame(index=mutations_dataframe.iloc[2:,:].index, columns=mutations_dataframe.columns)
for key in all_energies:
    position = all_energies[key].loc[0,"position_pre"]
    for idx in all_energies[key].index:
        pdb = all_energies[key].loc[idx, "pdb" ]
        per_res_energy_pre.loc[pdb,position] = all_energies[key].loc[idx,"energy_difference_pre"]
        per_res_energy_post.loc[pdb,position] = all_energies[key].loc[idx,"energy_difference_post"]

per_res_energy_pre.loc[:,["raw_score_prefusion", "raw_score_postfusion"]] = mutations_dataframe.loc[per_res_energy_pre.index,["raw_score_prefusion", "raw_score_postfusion"]]
per_res_energy_post.loc[:,["raw_score_prefusion", "raw_score_postfusion"]] = mutations_dataframe.loc[per_res_energy_post.index,["raw_score_prefusion", "raw_score_postfusion"]]
# Round values
per_res_energy_pre = per_res_energy_pre.applymap(custom_round)
per_res_energy_post = per_res_energy_post.applymap(custom_round)

per_res_energy_pre.fillna(value="-",inplace=True)      
per_res_energy_post.fillna(value="-",inplace=True)    

#write output
out_path_summary = os.path.join(out_path_results,"summary_sequences_and_per-residue_energy.xlsx")
output = open(out_path_summary, "w")
output.write("prefusion_rosetta_#,"+",".join(mutations_dataframe.columns.astype(str).tolist()))
output.write("\npostfusion_rosetta_#,"+",".join(mutations_dataframe.loc["postfusion_rosetta_#",:].astype(str).tolist()))
output.write("\nControl,"+",".join(mutations_dataframe.loc["control",:].astype(str).tolist()))

for index, row in mutations_dataframe.iloc[2:,:].iterrows():
    output.write("\n"+index+","+",".join(list(row.astype(str))))
    output.write("\n"+"Per-res_E_difference_Prefusion,"+",".join(per_res_energy_pre.loc[index,:].astype(str)))
    output.write("\n"+"Per-res_E_difference_Postfusion,"+",".join(per_res_energy_post.loc[index,:].astype(str)))
output.close()  

### take average energy difference per mutation             
aa = ['R','H','K','D','E','S','T','N','Q','C','G','P','A','V','I','L','M','F','Y','W']
average_energy_diff = {"prefusion":pd.DataFrame(columns=aa),"postfusion": pd.DataFrame(columns=aa) }

for key2 in sorted(all_energies, key=lambda x: int(x[1:])): #sort dictionary keys based on residue position
    data = all_energies[key2]
    data_split = data.groupby("mutation") #this does not have output but what it does is to create a virtual dataframe where the data is saved in bins, being each bin the unique values in "mutation" column. The bins are sorted A-Z
    unique_mut= data.mutation.unique().tolist()
    unique_mut.sort() #this is to know the order of the bins on data_split
    substitutions_pre = [data['energy_difference_pre'].values.tolist() for i, data in data_split] #obtain the values stored on each bin 
    substitutions_post = [data['energy_difference_post'].values.tolist() for i, data in data_split]

    average_pre = pd.DataFrame(columns=aa, index=[key2])
    # change residue numbering in postfusion dataframe
    position_post  = data.loc[0, "position_post"]
    if position_post != "-":
        key_post = key2[0]+str(position_post)       
    else:
        key_post = "-"   
    average_post = pd.DataFrame(columns=aa, index=[key_post])  
    for i in range(len(substitutions_pre)): 
        mean_value_pre = mean(substitutions_pre[i])
        mean_value_post = mean(substitutions_post[i])        
        average_pre.loc[key2,unique_mut[i]] = mean_value_pre
        average_post.loc[key_post,unique_mut[i]] = mean_value_post

    average_energy_diff["prefusion"] = average_energy_diff["prefusion"].append(average_pre)    
    average_energy_diff["postfusion"] = average_energy_diff["postfusion"].append(average_post)
    
average_energy_diff["prefusion"].to_excel(os.path.join(out_path_results,"average_per_resi_data_prefusion.xlsx"))
average_energy_diff["postfusion"].to_excel(os.path.join(out_path_results,"average_per_resi_data_postfusion.xlsx"))

###################  filtering
if arg_dict["filter"].lower() == "true": 
    average_energy_diff_filter_pre = average_energy_diff["prefusion"].copy()
    average_energy_diff_filter_post = average_energy_diff["postfusion"].copy()
    
    #based on energy changes in pre
    idx_to_discard_energy_pre = []
    idx_to_discard_energy_post = []
    i = -1 #counter for postfusion
    for idx, rows in average_energy_diff_filter_pre.iterrows():
        i+=1
        stabilizing = (rows <= -0.5).any() #check if at least one mutation is prefusion stabilizing
        if not stabilizing: #If no mutations are stabilizing, check if the ones with changes in energy less than 1 are postfusion destabilizing by at least 2.5 units. If so, do not discard
            discard = True
            for index, value in rows.items():
                if value < 1:
                    idx_post = average_energy_diff_filter_post.index[i]
                    if idx_post != "-":
                        value_post = average_energy_diff_filter_post.loc[idx_post,index]
                        if value_post >= 2.5:
                            discard = False
            if discard == True:
                idx_to_discard_energy_pre.append(idx)
                idx_to_discard_energy_post.append(i) #postfusion is based on index positions because it has missing values

    average_energy_diff_filter_pre.drop(index=idx_to_discard_energy_pre, inplace=True)
    tokeep_post = [i for i in range(len(average_energy_diff_filter_post)) if i not in idx_to_discard_energy_post]
    average_energy_diff_filter_post = average_energy_diff_filter_post.iloc[tokeep_post, :] 
    
    #output file
    average_energy_diff_filter_pre.to_excel(os.path.join(out_path_results,"average_per_resi_data_prefusion_FILTERED.xlsx"))
    average_energy_diff_filter_post.to_excel(os.path.join(out_path_results,"average_per_resi_data_postfusion_FILTERED.xlsx"))

############################ Examples on how to run per-residue metrics
#sc_hbonds_per_resi, hbonds_energies = per_resi_sidechain_hbonds_and_energy(pre_pose_wt, len(control_seq_pre), positions="")
#packing = per_resi_packing_score(pre_pose_wt, 1)





























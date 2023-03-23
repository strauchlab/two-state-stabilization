#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Karen J. Gonzalez
"""
import os
import pandas as pd
import numpy as np
from statistics import mean 
import re
import argparse
import subprocess

######################## Functions ########################

def expandSequence(sequence): #Sequence can be a list or a string
    seqExtended = []
    for i in sequence:
        if i == "R":
            seqExtended.append("ARG")
        elif i == "H":
            seqExtended.append("HIS")    
        elif i == "K":
            seqExtended.append("LYS")            
        elif i == "D":
            seqExtended.append("ASP")            
        elif i == "E":
            seqExtended.append("GLU")            
        elif i == "S":
            seqExtended.append("SER")            
        elif i == "T":
            seqExtended.append("THR")            
        elif i == "N":
            seqExtended.append("ASN")    
        elif i == "Q":
            seqExtended.append("GLN")    
        elif i == "C":
            seqExtended.append("CYS")    
        elif i == "G":
            seqExtended.append("GLY")            
        elif i == "P":
            seqExtended.append("PRO")            
        elif i == "A":
            seqExtended.append("ALA")            
        elif i == "V":
            seqExtended.append("VAL")            
        elif i == "I":
            seqExtended.append("ILE")            
        elif i == "L":
            seqExtended.append("LEU")    
        elif i == "M":
            seqExtended.append("MET")    
        elif i == "F":
            seqExtended.append("PHE")    
        elif i == "Y":
            seqExtended.append("TYR")
        elif i == "W":
            seqExtended.append("TRP") 
    return seqExtended

def contractSeq(sequence): #sequence has to be a list with each amino acid as an independent element
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

def mut_file_trimer(positions, seq, path_output, alanine=False, all_aa=True, include_original=False):
    if alanine==True and all_aa==True: #this is to calculate ddg for the 20 amino acids except the original amino acid in the sequence (unless include_original=True)
        total = 57
        aa = ['R','H','K','D','E','S','T','N','Q','C','G','P','A','V','I','L','M','F','Y','W']
    elif alanine==False and all_aa==True: #If you do not want to calculate alanine again
        total = 54
        aa = ['R','H','K','D','E','S','T','N','Q','C','G','P','V','I','L','M','F','Y','W']
    else:
        total = 3
        aa = ['A'] 

    if include_original == True: 
        total+= 3
        total = str(total)
    else:
        total = str(total)
        
    os.chdir(path_output)
    trimer_positions = list(range(1,len(seq)+1))
    trimer_positions = trimer_positions + trimer_positions + trimer_positions 
    trimer_seq = seq + seq + seq

    positions_idx = []
    for i in positions: #these are the positions we want to modify
        resi_idx = [j for j, x in enumerate(trimer_positions) if x == i] #this is to find all the indexes where i appears.
        positions_idx.append(resi_idx)
    
    for i in range(len(positions_idx)):
        doc = open("mut_file"+str(i+1),"w")
        if trimer_seq[positions_idx[i][0]] == "A" and alanine==False: #when we want to calculate all aa and do not repeat alanine, if the native aa is an ala, then the total number of mutations is 57 and not 54
            doc.write("total "+str(int(total)+3) +"\n3\n") 
        else:
            doc.write("total "+total+"\n3\n")
        if include_original == True:
            for j in range(len(aa)):
                if j == len(aa)-1 :
                    doc.write(trimer_seq[positions_idx[i][0]]+' '+str(positions_idx[i][0]+1)+' '+aa[j])
                    doc.write('\n'+trimer_seq[positions_idx[i][1]]+' '+str(positions_idx[i][1]+1)+' '+aa[j])
                    doc.write('\n'+trimer_seq[positions_idx[i][2]]+' '+str(positions_idx[i][2]+1)+' '+aa[j])
                else:
                    doc.write(trimer_seq[positions_idx[i][0]]+' '+str(positions_idx[i][0]+1)+' '+aa[j])
                    doc.write('\n'+trimer_seq[positions_idx[i][1]]+' '+str(positions_idx[i][1]+1)+' '+aa[j])
                    doc.write('\n'+trimer_seq[positions_idx[i][2]]+' '+str(positions_idx[i][2]+1)+' '+aa[j])
                    doc.write('\n3\n')
        else:
            for j in range(len(aa)):
                if aa[j] != trimer_seq[positions_idx[i][0]]:
                    if j == len(aa)-1 :
                        doc.write(trimer_seq[positions_idx[i][0]]+' '+str(positions_idx[i][0]+1)+' '+aa[j])
                        doc.write('\n'+trimer_seq[positions_idx[i][1]]+' '+str(positions_idx[i][1]+1)+' '+aa[j])
                        doc.write('\n'+trimer_seq[positions_idx[i][2]]+' '+str(positions_idx[i][2]+1)+' '+aa[j])
                    else:
                        doc.write(trimer_seq[positions_idx[i][0]]+' '+str(positions_idx[i][0]+1)+' '+aa[j])
                        doc.write('\n'+trimer_seq[positions_idx[i][1]]+' '+str(positions_idx[i][1]+1)+' '+aa[j])
                        doc.write('\n'+trimer_seq[positions_idx[i][2]]+' '+str(positions_idx[i][2]+1)+' '+aa[j])
                        doc.write('\n3\n')
        doc.close()
                  
def extractEnergy(files,seq):
    ddgs = []
    wild_type_energies = []
    mutation_energies = []
    for j in range(len(files)):
        with open(files[j],'r') as fopen:
            lines = fopen.readlines()
            lines = [x.split() for x in lines]
            limit = len(lines)
            if limit > 3: #ignore files where the ddg  was only calculated  for WT
                values_WT = [float(lines[0][3]),float(lines[1][3]),float(lines[2][3])]
                wild_type_energies.append(values_WT)
                averageWT = mean(values_WT)
                position = lines[3][2].split('_') #split the line where there is the mutation, e.g "MUT_25ARG_141ARG". result = [MUT,25ARG,141ARG]
                position = re.findall(r"\d+",position[1])#split the position 1 by the digits (the 'r"\d+"' means regular expression to find only digits).
                position = int(position[0])-1 #position[0] because position is a list

            for i in range(3,limit,3):
                mutation = lines[i][2].split('_')
                mutation = re.findall(r"\D+",mutation[1])#split the position 1 by non-digit characters (to obtain only the letters for example)
                mutation = re.findall(r"\w+",mutation[0]) #gives us the aa identity
                values_MUT = [float(lines[i][3]),float(lines[i+1][3]),float(lines[i+2][3])]
                mutation_energies.append([position+1, mutation[0], values_MUT] ) #+1 because position was indicating an index previously
                averageMUT = mean(values_MUT)
                ddg = averageMUT - averageWT
                ddgs.append([position+1, mutation[0], [ddg]] )
    return ddgs, wild_type_energies, mutation_energies
 
def statisticts(data, data_frame = False): 
    # this analysis was taken from https://towardsdatascience.com/ways-to-detect-and-remove-the-outliers-404d16608dba          
    # we are applying The interquartile range (IQR) to identify outliers

    #data has to be a dataframe
    #this calculates quartiles over all columns independently. However, for this case
    # we are just entering one column, that is, a Series object instead of a dataframe       
          
    Q1 = data.quantile(0.25)
    Q3 = data.quantile(0.75)
    IQR = Q3 - Q1                   
    outlier_lower = data < (Q1 - 1.5 * IQR) #this returns a dataframe of Trues and Falses. True for the values that are lower than the lowest limit
    outlier_upper = data > (Q3 + 1.5 * IQR)
    
    if data_frame == False:
        for index, value in outlier_lower.iteritems():
            if value == True or outlier_upper[index] == True:
                data[index] = np.nan
    else:
        #for a dataframe with several columns        
        for key, value in outlier_lower.iteritems(): 
            for index, row in outlier_lower.iterrows():
                if outlier_lower.loc[index, key] == True or outlier_upper.loc[index, key] == True:
                    data.loc[index, key] = np.nan          
    # average of each column
    average = data.mean(axis=0)
    return average

def collect_data(root, folder, seq):
    aa = ['ARG','HIS','LYS','ASP','GLU','SER','THR','ASN','GLN','CYS','GLY','PRO','ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP']
    limit = len(seq)+1 
    os.chdir(root)
    
    dictionary = {} 
    for i in range(1,limit):  #create keys of each position
        dictionary[i]= {}
        for amino in aa:
            dictionary[i][amino]= [] 

    for f in folder:
        directory = os.path.join(root, f)
        os.chdir(directory)
        files = [i for i in os.listdir(".") if i.endswith(".ddg")]
        ddgs,_,_ = extractEnergy(files,seq)
        aa_dict = {}
        
        for i in range(len(ddgs)):
            aa_dict[ddgs[i][1]] = ddgs[i][2] #[1] contains the aa identity
            try:
                if ddgs[i][0] != ddgs[i+1][0]: #check when the position changes
                    for key, values in aa_dict.items():
                        for value in values:
                            dictionary[ddgs[i][0]][key].append(value) #this might throw error for the first iteration, because dictionary_pre is empty
                    aa_dict = {}             
            except:
                for key, values in aa_dict.items():
                    for value in values:
                        dictionary[ddgs[i][0]][key].append(value) #this might throw error for the first iteration, because dictionary_pre is empty
                    aa_dict = {}
             #the last row will throw an error, because i+1 is out of range. That means we reached the last set of positions 
    return dictionary

def results(seq, dictionary, pdb_positions, scaling_factor = 1.0/2.94  ):
    aa = ['ARG','HIS','LYS','ASP','GLU','SER','THR','ASN','GLN','CYS','GLY','PRO','ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP']
    energies = [['' for i in range(20)] for n in range(len(seq))]    
    
    for j in range(1,len(dictionary)+1):
        for key in dictionary[j]:
            idxMutation = aa.index(key)
            if dictionary[j][key] != []:
                data_average = pd.Series(dictionary[j][key])
                averageMUT =  statisticts(data_average)
                energies[j-1][idxMutation] = averageMUT
            else:
               energies[j-1][idxMutation] = 0
    
    resultsALL = pd.DataFrame(energies, columns=aa)  
    rosetta_positions = list(range(1, len(seq)+1))           
    resultsALL = resultsALL.set_index([pd.Index(rosetta_positions)])
    scaled_results= resultsALL.multiply(scaling_factor)
    scaled_results['Seq'] = seq
    scaled_results['rosetta_#'] = rosetta_positions  
    scaled_results['PDB_#'] = pdb_positions
    return scaled_results 

def alanine_comparison_pre_post(alanines_pre,alanines_post): #datasets to be compared neeed to have the same indexes 
    ala_Pre_Post = pd.concat([alanines_pre.loc[:,"ALA"], alanines_post.loc[:,"ALA"]], axis=1)
    ala_Pre_Post.columns = ['Ala_Pre','Ala_Post']
    ala_Pre_Post["difference"] = ala_Pre_Post['Ala_Post'].subtract(ala_Pre_Post['Ala_Pre'])            
    hot_spots = pd.DataFrame()
    for index, row in ala_Pre_Post.iterrows():
        if row['Ala_Pre']<= 1.3 and row['Ala_Post']>=-1.3:
            if -1<= row['Ala_Pre']<= 1 and -1<= row['Ala_Post']<= 1: #if both pre and post values are in the neutral zone, do not consider them
                continue
            else:
                if row["difference"] >= 0.7: #there are some values that do not fall into the neutral zone but are not good enough because they are close between pre and post. For example pre = 0.8 and post 1.3. Also this removes the native alanines
                    hot_spots = hot_spots.append(row)
    
    hotspots_pre = []
    hotspots_post = []
    
    for idx,row in hot_spots.iterrows():        
        hotspots_pre.append(alanines_pre.loc[idx,"rosetta_#"])
        hotspots_post.append(alanines_post.loc[idx,"rosetta_#"])
        
    return hotspots_pre, hotspots_post, ala_Pre_Post, hot_spots                                                                                       
                                                
################################################
root = os.getcwd()   
### Command line arguments 
parser = argparse.ArgumentParser(description='Comparison ddg alanine scanning prefusion vs postfusion')
parser.add_argument('arg_file', help = "File containing all input arguments")
args = parser.parse_args()
# Input arguments from arg_file
arg_dict = {}
with open(args.arg_file, "r") as fopen:
    for line in fopen:
        if line and line[0].isalpha():
            line = line.split("=")
            arg_dict[line[0].strip()] = line[1].strip()
   
### Postfusion alanine ddg results
post_directory = arg_dict["post_dir"]
pdb_post = arg_dict["pdb_post"]
seq_post_monomer, pdb_positions_post,_ = sequence_from_pdb(pdb_post, chain=arg_dict["post_monomer_ch"]) 
folder = arg_dict["post_ala_folder"].split(",")
# Calculate ddg per experiment. Then filter out outlier ddgs and take average between all ddgs of the same mutation
data_post_dict = collect_data(post_directory, folder, seq_post_monomer)
results_post = results(seq_post_monomer, data_post_dict, pdb_positions_post)
os.chdir(root)

### Prefusion alanine ddg results
pre_directory = arg_dict["pre_dir"]
pdb_pre = arg_dict["pdb_pre"]
seq_pre_monomer, pdb_positions_pre,_ = sequence_from_pdb(pdb_pre, chain=arg_dict["pre_monomer_ch"]) 
folder = arg_dict["pre_ala_folder"].split(",")
# Calculate ddg per experiment. Then filter out outlier ddgs and take average between all ddgs of the same mutation
data_pre_dict = collect_data(pre_directory, folder, seq_pre_monomer)
results_pre = results(seq_pre_monomer, data_pre_dict, pdb_positions_pre)
os.chdir(root)

### only consider residue positions that are shared between pre- and postfusion structures
alignment = read_alignment(arg_dict["alignment"])
alignment = alignment.rename(columns={arg_dict["pre_name_alignment"]: "prefusion", arg_dict["post_name_alignment"]: "postfusion"})
alignment = add_data_according_to_alignment("postfusion", alignment, list(range(1, len(seq_post_monomer)+1)), column_name="rosetta_#" )
alignment = add_data_according_to_alignment("prefusion", alignment, list(range(1, len(seq_pre_monomer)+1)), column_name="rosetta_#" )
  
idx_shared = []
for index, row in alignment.iterrows():
    if alignment.loc[index,"prefusion"] != "-" and alignment.loc[index,"postfusion"] != "-":
        idx_shared.append(index)
                                                                                         
alanines_post = results_post.loc[alignment.loc[idx_shared,"postfusion_rosetta_#"].tolist(),["ALA","rosetta_#","Seq"]]
alanines_pre = results_pre.loc[alignment.loc[idx_shared,"prefusion_rosetta_#"].tolist(),["ALA","rosetta_#","Seq"]]

### Find positions where ddg stabilizing_pre/ destabilizing_post or neutral_pre/ destabilizing_post      
# comparison ala scan pre-post
alanines_post_newIdx = alanines_post.set_index(pd.Index(alanines_pre.index.tolist()),inplace=False) # we need to have the two dataframes with the same indexes to be able to substract their values. The indexes of prefusion were chosen here                                          
hotspots_pre, hotspots_post, ala_pre_post, hotspots_data = alanine_comparison_pre_post(alanines_pre,alanines_post_newIdx) #hotspots are in rosetta numbering (first monomer only)

# Include positions naturally bearing alanines (because native alanines are never tested)
for index,row in alignment.loc[idx_shared, :].iterrows(): 
    if alignment.loc[index, "prefusion"] == "A":
        hotspots_pre.append(alignment.loc[index, "prefusion_rosetta_#"])
        hotspots_post.append(alignment.loc[index, "postfusion_rosetta_#"])                                           

# Include the fusion peptide residues where alanine scanning shows a stabilizing effect
fusion_peptide = arg_dict["fusion_peptide"].split("-")
fp_idx_start = results_pre[results_pre["PDB_#"] == int(fusion_peptide[0])].index.tolist()[0]
fp_idx_end = results_pre[results_pre["PDB_#"] == int(fusion_peptide[1])].index.tolist()[0]
fusion_peptide = results_pre.loc[fp_idx_start:fp_idx_end,["ALA","rosetta_#","Seq"]]
hot_spots_FP = fusion_peptide[(fusion_peptide['ALA']<= -0.7) | (fusion_peptide['Seq']== "A")] #-0.7 because the selection process for the other parts of the protein stablished a minimum energetic difference (post-pre) of 0.7. Since this region is not present in the postfusion state, it is asummed that the region has a ddg value of 0 in the postusion state

hotspots_pre+= hot_spots_FP.loc[:, "rosetta_#"].tolist()
hotspots_pre.sort()
hotspots_post.sort()

### write mut_file to perform all-amino acids scanning at hotspots       
output_dir_pre = os.path.join(root, "cartesian_all-aa_pre")
output_dir_post = os.path.join(root, "cartesian_all-aa_post")
subprocess.call(["mkdir",output_dir_pre]); subprocess.call(["mkdir",output_dir_post])
mut_file_trimer(hotspots_pre, seq_pre_monomer, output_dir_pre, alanine=True, all_aa=True) 
mut_file_trimer(hotspots_post, seq_post_monomer, output_dir_post, alanine=True, all_aa=True) 


            
            
 




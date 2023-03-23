#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Karen J. Gonzalez
"""
import os
import pandas as pd
import numpy as np
from statistics import mean 
import math
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

def insert_rows(alignment, state, results):
    hole = False
    start = '' 
    alignment_idx = []
    for index, row in alignment.iterrows():
        if alignment.loc[index,state] == "-":
            if hole == False:           
                start = alignment.loc[index-1,state+"_rosetta_#"]
                alignment_idx.append(index)   
                hole = True  
            else:
                alignment_idx.append(index)  
        else:
            if start != "":
                end =  alignment.loc[index,state+"_rosetta_#"]
                #add extra residues to results
                extra = len(alignment_idx) #this gives the total extra residues that need to be added. 
                extra = [[0 for i in range(24)] for j in range(extra)]   # 24 is the number of columns in the results dataframe           
                idxs = ["00_"+str(i) for i in alignment_idx]
                extra = pd.concat([results.loc[:start, :], pd.DataFrame(extra,index=idxs, columns=results.columns)], ignore_index=False) 
                results = pd.concat([extra,results.loc[end:,:]], ignore_index=False)
                #reset variables
                hole = False   
                alignment_idx = []
                start =''
                
    if start != "": #if the start variable ends non-empty is because the extra residues are at the end of the sequence.
        #add extra residues to results
        extra = len(alignment_idx) #this gives the total extra residues that need to be added. 
        extra = [[0 for i in range(24)] for j in range(extra)]   # 24 is the number of columns in the results dataframe           
        idxs = ["00_"+str(i) for i in alignment_idx]
        extra = pd.DataFrame(extra,index=idxs, columns=results.columns)
        results = pd.concat([results, extra], ignore_index=False) 
    return results                                                                                 
                                                                                           
def pssm_matrix_V2(results_pre_all, results_post_all, surface_resi, threshold_difference=0.7):
    aa = ['ARG','HIS','LYS','ASP','GLU','SER','THR','ASN','GLN','CYS','GLY','PRO','ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP']
    pssm_matrix = pd.DataFrame(columns=aa, index=results_pre_all.index.tolist())
    results_post_all_newIdx = results_post_all.set_index(pd.Index(results_pre_all.index.tolist()),inplace=False)#we need to have the two dataframes with the same indexes
    differences = results_post_all_newIdx.loc[:, :"TRP"].subtract(results_pre_all.loc[:, :"TRP"])   
    for index, row in differences.iterrows():
        if index in surface_resi: #discard positions that are surface-exposed in both pre and postfusion structures
                pssm_matrix.loc[index,:] = 0.05 # pseudocount to avoid errors because of log(0)
                continue
        for key, value in differences.iteritems():                   
            if results_pre_all.loc[index,key] >= 1: #discard positions that destabilize the prefusion state
                pssm_matrix.loc[index,key] = 0.05 # pseudocount to avoid errors because of log(0)
            elif differences.loc[index,key] < threshold_difference: #discard positions that have a difference in energy pre vs post less than x units
                pssm_matrix.loc[index,key] = 0.05 # add pseudocounts.  
            else:
                pssm_matrix.loc[index,key] = differences.loc[index,key] + 0.05 #add pseudocounts not to alter the preference
    
    pssm_log = pssm_matrix.copy()
    for index, row in pssm_matrix.iterrows():
        sum_row = sum(row[:])                
        for key, value in pssm_matrix.iteritems():
            pssm_matrix.loc[index,key] = pssm_matrix.loc[index,key]/sum_row  
            pssm_log.loc[index,key] = math.log2(pssm_matrix.loc[index,key]/0.05)                      
    return pssm_matrix, pssm_log, differences 

def pssm_output(root, pssm_log, pssm_matrix, results_all, pssm_name="pre.pssm"):
    new_dir = os.path.join(root,"combinatorial_design")
    subprocess.call(["mkdir", new_dir])
    os.chdir(new_dir)
    aa_psiblast = "ARNDCQEGHILKMFPSTWYV" #order of aa in the psi blast
    aa_psiblast_ext = expandSequence(aa_psiblast)
    pssm_log = pssm_log[aa_psiblast_ext] #change the order of the columns   
    pssm_log.index = results_all.index.tolist() #to syncronize with the indexes in the results dataframe
    pssm_matrix = pssm_matrix[aa_psiblast_ext]  
    pssm_matrix.index = results_all.index.tolist()#to syncronize with the indexes in the results dataframe
    #write PSSM
    pssm = open(pssm_name, 'w')
    pssm.write("\n")
    pssm.write("log_transformed probabilities, probabilities\n")
    pssm.write("            ")
    for i in aa_psiblast:
        pssm.write(i+'   ')
    pssm.write(" ")    
    for i in aa_psiblast:
        pssm.write(i+'   ')    
    pssm.write("\n")  
    
    for index, row in results_all.iterrows():
        if type(index) == str: #discard the extra residues that were added to be able to process the results
            continue
        else:
            pssm.write('    ')
            pssm.write(str(results_all.loc[index,"rosetta_#"])+" "+results_all.loc[index,"Seq"]+"    ")
            for key, value in pssm_matrix.iteritems():
                pssm.write(str(round(pssm_log.loc[index, key], 4))+'   ')
                pssm.write(" ")  
            for key, value in pssm_matrix.iteritems():   
                pssm.write(str(round(100*pssm_matrix.loc[index, key], 4))+'   ')
            pssm.write("\n") 
    pssm.close()   

def resfile_V2_2(root, differences, results_all_pre, results_all_post, surface_resi, threshold_difference=0.7, threshold_stability = 0.7, resfile_name='resfile_pre'):
    new_dir = os.path.join(root,"combinatorial_design")
    os.chdir(new_dir)
    results_all_post_newIdx = results_all_post.set_index(pd.Index(results_all_pre.index.tolist()),inplace=False)    
    resfile = open(resfile_name,'w')
    resfile.write("nataa\nstart\n")  
    positions_pre = []
    positions_post = []
    for index, row in differences.iterrows():
        if index not in surface_resi:
            amino_acid = ''
            for key, value in differences.iteritems():
                if results_all_pre.loc[index,key] <=1: # discard all substitutions that destabilize the prefusion state
                    if differences.loc[index,key] > threshold_difference: #this corresponds to the values with a ddg difference greater than threshold_difference
                        if results_all_pre.loc[index,key] < 0 and results_all_post_newIdx.loc[index,key] < 0: # if both states are stabilized, only accept the substitution if the stabilizing effect in prefusion is higher than in postfusion by more than the threshold_stability
                            if differences.loc[index,key] > threshold_stability:
                                amino_acid += contractSeq([key])[0]
                        else:
                            amino_acid += contractSeq([key])[0]
                    elif contractSeq([key])[0] == results_all_pre.loc[index,"Seq"]:#the native amino acids. Since contract gives a list, we need to put [0]
                        amino_acid += contractSeq([key])[0]
            if len(amino_acid) > 1:   #This is to only include residues that have more substitutions than just their native aa    
                resfile.write(str(results_all_pre.loc[index,"PDB_#"])+'\t'+results_all_pre.loc[index,"chain"]+'\tPIKAA\t'+amino_acid+"\n") 
                positions_pre.append(index)
                position_post = results_all_post_newIdx.loc[index,"rosetta_#"]
                if position_post != 0:
                    positions_post.append(position_post)
    resfile.close()   
    return positions_pre, positions_post  

def resfile_control(root, positions, results_all, name_resfile="resfile_control"):
    new_dir = os.path.join(root,"combinatorial_design")
    os.chdir(new_dir)
    resfile = open(name_resfile,'w')
    resfile.write("nataa\nstart\n")     
    for index in positions:
        resfile.write(str(results_all.loc[index,"PDB_#"])+'\t'+results_all.loc[index,"chain"]+'\tPIKAA\t')
        amino_acid = results_all.loc[index,"Seq"] #this corresponds to the native amino acids.
        resfile.write(amino_acid)
        resfile.write("\n")    
    resfile.close()  

###############################################################################
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
    
### Load all- amino acids scanning results         
## Postfusion ddg results
post_directory = arg_dict["post_dir"]
pdb_post = arg_dict["pdb_post"]
seq_post_monomer, pdb_positions_post, chain_monomer_post = sequence_from_pdb(pdb_post, chain=arg_dict["post_monomer_ch"]) 
folder = arg_dict["post_ddg_folder"].split(",")
# Calculate ddg per experiment. Then filter out outlier ddgs and take average between all ddgs of the same mutation
data_post_dict = collect_data(post_directory, folder, seq_post_monomer)
results_post = results(seq_post_monomer, data_post_dict, pdb_positions_post)
results_post["chain"] = chain_monomer_post
os.chdir(root)

## Prefusion ddg results
pre_directory = arg_dict["pre_dir"]
pdb_pre = arg_dict["pdb_pre"]
seq_pre_monomer, pdb_positions_pre, chain_monomer_pre = sequence_from_pdb(pdb_pre, chain=arg_dict["pre_monomer_ch"]) 
folder = arg_dict["pre_ddg_folder"].split(",")
# Calculate ddg per experiment. Then filter out outlier ddgs and take average between all ddgs of the same mutation
data_pre_dict = collect_data(pre_directory, folder, seq_pre_monomer)
results_pre = results(seq_pre_monomer, data_pre_dict, pdb_positions_pre)
results_pre["chain"] = chain_monomer_pre
os.chdir(root)

## Since the pre and postfusion structures contain different number of residues, we need to add these extra residues to each results_dataframe
#Identify extra residues based on sequence alignment 
alignment = read_alignment(arg_dict["alignment"])
alignment = alignment.rename(columns={arg_dict["pre_name_alignment"]: "prefusion", arg_dict["post_name_alignment"]: "postfusion"})
alignment = add_data_according_to_alignment("postfusion", alignment, list(range(1, len(seq_post_monomer)+1)), column_name="rosetta_#" )
alignment = add_data_according_to_alignment("prefusion", alignment, list(range(1, len(seq_pre_monomer)+1)), column_name="rosetta_#" )
alignment = add_data_according_to_alignment("postfusion", alignment, pdb_positions_post, column_name="PDB_#" )
alignment = add_data_according_to_alignment("prefusion", alignment, pdb_positions_pre, column_name="PDB_#" )

#Add extra residues to results dataframe                                          
results_post_all = insert_rows(alignment, "postfusion",results_post)
results_pre_all = insert_rows(alignment, "prefusion",results_pre)

## Discard amino acid positions that are surface-exposed in both pre- and postfusion structures
surface_resi =[]
if arg_dict["discard_surf_resi"] == "True":
    from pyrosetta import *
    from rosetta import *
    from rosetta.core.scoring.sasa import SasaCalc
    init()   
    pose_pre = pose_from_pdb(pdb_pre)
    pose_post = pose_from_pdb(pdb_post) 
    #set up a residue selector operation to select only surface residues
    resi_selector_surface = core.select.residue_selector.LayerSelector()
    resi_selector_surface.set_layers(False,False,True) #pick_core: False, pick_boundary: False, pick_surface: True
    #Apply the residue selector operation to pre and post poses
    surface_resi_pre = resi_selector_surface.apply(pose_pre)
    surface_resi_pre = [i for i in surface_resi_pre] #this gives a list of Trues and Falses
    surface_resi_post = resi_selector_surface.apply(pose_post)
    surface_resi_post = [i for i in surface_resi_post]
    #add surface exposure data to the alignment dataframe
    alignment = add_data_according_to_alignment("postfusion", alignment, surface_resi_post, column_name="surface" )
    alignment = add_data_according_to_alignment("prefusion", alignment, surface_resi_pre, column_name="surface" )
    # Identify surface-exposed residues (prefusion- rosetta numbering)
    for index, row in alignment.iterrows():
        if alignment.loc[index, "prefusion_surface"] == True and alignment.loc[index, "postfusion_surface"] == True:
            surface_resi.append(alignment.loc[index, "prefusion_rosetta_#"])
                                            
## Write PSSM
pssm_matrix, pssm_log, differences = pssm_matrix_V2(results_pre_all, results_post_all, surface_resi, threshold_difference=2 )
pssm_output(root, pssm_log, pssm_matrix, results_pre_all, pssm_name="pre.pssm")
pssm_output(root, pssm_log, pssm_matrix, results_post_all, pssm_name="post.pssm")

#Write resfile for prefusion combinatorial design
designable_idx_pre,designable_idx_post  = resfile_V2_2(root, differences, results_pre_all, results_post_all, surface_resi, threshold_difference=2,threshold_stability = 4, resfile_name='resfile_pre')
#Write control resfile for prefusion
resfile_control(root, designable_idx_pre, results_pre_all, name_resfile="resfile_control_pre")






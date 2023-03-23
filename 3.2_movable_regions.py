#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Karen J. Gonzalez
"""
import os
import pandas as pd
import math
import argparse
import subprocess
from pyrosetta import *
from rosetta import *

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

def RMSD(vector1,vector2):
    euclidean_distance = ((vector1[0]-vector2[0])**2) + ((vector1[1]-vector2[1])**2) + ((vector1[2]-vector2[2])**2)
    rmsd = math.sqrt((euclidean_distance/2))
    return rmsd

def ca_coordinates(pdb, chain="A"):  
    data = []
    rosetta_position = 0
    with open(pdb,'r') as fopen:
        for line in fopen:
            if line.startswith("ATOM") and line[13:15] == "CA":
                if line[21] in chain:
                    position = int(line[22:26])
                    rosetta_position+=1
                    sequence = contractSeq([line[17:20]])  
                    data.append([float(line[30:38]), float(line[38:46]), float(line[46:54]), sequence[0], position, rosetta_position])     
        fopen.close()
    coordinates = pd.DataFrame(data, columns=["x","y","z","seq","PDB_#","rosetta_#"])
    coordinates.index = coordinates["rosetta_#"].tolist()
    return coordinates

def round_half_up(n, decimals=0):
    multiplier = 10 ** decimals
    return math.floor(n*multiplier + 0.5 )/multiplier

        
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
                                           
################################################
root = os.getcwd()    
### Command line arguments 
parser = argparse.ArgumentParser(description='Finding highly movable regions')
parser.add_argument('arg_file', help = "File containing all input arguments")
args = parser.parse_args()
# input arguments from arg_file
arg_dict = {}
with open(args.arg_file, "r") as fopen:
    for line in fopen:
        if line and line[0].isalpha():
            line = line.split("=")
            arg_dict[line[0].strip()] = line[1].strip()

### Obtain ca coordinates
pdb_post = arg_dict["pdb_post"]
monomer_chain_post = arg_dict["post_monomer_ch"]
coordinates_post = ca_coordinates(pdb_post, chain= monomer_chain_post)

pdb_pre = arg_dict["pdb_pre"]
monomer_chain_pre = arg_dict["pre_monomer_ch"]
coordinates_pre = ca_coordinates(pdb_pre, chain= monomer_chain_pre)
# Only consider residue positions that are shared between pre- and postfusion structures
alignment = read_alignment(arg_dict["alignment"])
alignment = alignment.rename(columns={arg_dict["pre_name_alignment"]: "prefusion", arg_dict["post_name_alignment"]: "postfusion"})
alignment = add_data_according_to_alignment("prefusion", alignment, coordinates_pre.loc[:,"rosetta_#"].tolist(), column_name="rosetta_#")
alignment = add_data_according_to_alignment("postfusion", alignment, coordinates_post.loc[:,"rosetta_#"].tolist(), column_name="rosetta_#")
                                                   
idx_shared = []
for index, row in alignment.iterrows():
    if alignment.loc[index,"prefusion"] != "-" and alignment.loc[index,"postfusion"] != "-":
        idx_shared.append(index)
                                                    
# Trim coordinates to shared positions only
coordinates_pre_shared = coordinates_pre.loc[alignment.loc[idx_shared,"prefusion_rosetta_#"].tolist(), :]
coordinates_post_shared = coordinates_post.loc[alignment.loc[idx_shared,"postfusion_rosetta_#"].tolist(),:].set_index(pd.Index(coordinates_pre_shared.index.tolist()), inplace=False) # we need to have the two dataframes with the same indexes. The indexes of prefusion were chosen here  

### Calculate per residue root mean square deviation 
distance = []
for index, row in coordinates_pre_shared.iterrows():
    data1 = coordinates_pre_shared.loc[index, ["x","y","z"]].tolist()
    data2 = coordinates_post_shared.loc[index, ["x","y","z"]].tolist()
    rmsd = RMSD(data1,data2)
    distance.append(rmsd)
    
### Find positions that move more than ~10 Angs
mobile_pre = []
mobile_post = []

for i in range(len(distance)):
    if round_half_up(distance[i]) >= 10 :
        mobile_pre.append(coordinates_pre_shared.iloc[i,-1])
        mobile_post.append(coordinates_post_shared.iloc[i,-1])

### Add flanking regions
'''We have also included regions flanking the mobile areas where the secondary structure differs between the pre- and postfusion structures.'''
''' we search until finding a block of X (custom) amino acids matching their secondary strucures. '''

# Identify secondary structure of input PDBs
init() #initialize pyrosetta
pose_pre = pose_from_pdb(pdb_pre)
pose_post = pose_from_pdb(pdb_post)
dssp = pyrosetta.rosetta.protocols.moves.DsspMover()
dssp.apply(pose_pre)    
dssp.apply(pose_post)
ss_pre = pose_pre.secstruct() #secondary structure prediction for the prefusion state
ss_post = pose_post.secstruct() #secondary structure prediction for the postfusion state

# Add secondary structure information to aligned sequences
alignment = add_data_according_to_alignment("prefusion", alignment, ss_pre, column_name="ss")
alignment = add_data_according_to_alignment("postfusion", alignment, ss_post, column_name="ss")

# Restrict the analysis to shared regions
alignment_shared = alignment.loc[idx_shared,:]

# Find flanking regions 
flanking_pre = []
flanking_post  = []
resi_flanking = int(arg_dict["flanking"])
for i in range(len(mobile_pre)):
    if i+1 < len(mobile_pre):
        if mobile_pre[i] +1 != mobile_pre[i+1]  : #find where the mobile area is discontinuos
            #analyze upstream neighbors
            neighbor_up = mobile_pre[i] +1     
            if neighbor_up in alignment_shared['prefusion_rosetta_#'].tolist():
                idx_neighbor_up = alignment_shared[alignment_shared['prefusion_rosetta_#'] == neighbor_up].index.tolist()[0]
                for idx in range(idx_neighbor_up, alignment_shared.index[-1]+1):
                    if alignment_shared.loc[idx:idx+resi_flanking, "prefusion_ss"].tolist() == alignment_shared.loc[idx:idx+resi_flanking, "postfusion_ss"].tolist():
                        flanking_pre += alignment_shared.loc[idx_neighbor_up:idx+resi_flanking, 'prefusion_rosetta_#'].tolist()
                        flanking_post += alignment_shared.loc[idx_neighbor_up:idx+resi_flanking, 'postfusion_rosetta_#'].tolist()                                                  
                        break
                              
# Update mobile regions with flanking residues                                 
mobile_pre += flanking_pre ; mobile_post += flanking_post        

# Include the fusion peptide residues 
fusion_peptide = arg_dict["fusion_peptide"].split("-")
fp_idx_start = coordinates_pre[coordinates_pre["PDB_#"] == int(fusion_peptide[0])].index.tolist()[0]
fp_idx_end = coordinates_pre[coordinates_pre["PDB_#"] == int(fusion_peptide[1])].index.tolist()[0]
fusion_peptide = coordinates_pre.loc[fp_idx_start:fp_idx_end,"rosetta_#"].tolist()
mobile_pre+= fusion_peptide

mobile_pre = list(set(mobile_pre)); mobile_pre.sort()
mobile_post = list(set(mobile_post)); mobile_post.sort()
        
#### Write mut_file to perform all-amino acids scanning at hotspots       
output_dir_pre = os.path.join(root, "cartesian_all-aa_pre")
output_dir_post = os.path.join(root, "cartesian_all-aa_post")
subprocess.call(["mkdir",output_dir_pre]); subprocess.call(["mkdir",output_dir_post])
mut_file_trimer(mobile_pre, coordinates_pre.loc[:,"seq"].tolist(), output_dir_pre, alanine=True, all_aa=True)
mut_file_trimer(mobile_post, coordinates_post.loc[:,"seq"].tolist(), output_dir_post, alanine=True, all_aa=True) 
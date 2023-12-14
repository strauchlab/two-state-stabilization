#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Karen J. Gonzalez
"""
import os
import argparse
import subprocess

######################## Functions ########################

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

def mut_file_trimer(positions, seq, path_output, alanine=True, all_aa=False, include_original=False):
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
parser = argparse.ArgumentParser(description='Writting mut_files for alanine scanning')
parser.add_argument('pdb', help = "File name or directory of the PDB to be analyzed")
parser.add_argument('ch', help = "Chain(s) ID of one protomer")
args = parser.parse_args()

## Extract sequence information
seq_monomer, _,_ = sequence_from_pdb(args.pdb, chain=args.ch) 
rosetta_positions = list(range(1,len(seq_monomer)+1))

### write mut_file to perform alanine scanning     
output_dir = os.path.join(root, "mut_files")
subprocess.call(["mkdir",output_dir])
mut_file_trimer(rosetta_positions, seq_monomer, output_dir) 



            
            
 




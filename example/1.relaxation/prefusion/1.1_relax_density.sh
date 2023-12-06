#!/bin/bash

$Rosetta/main/source/bin/relax.linuxgccrelease \
-database $Rosetta/main/database \
-in:file:s 5wb0_BIOMT_expanded_ignorechain_ABC_INPUT.pdb \ #the monomeric structure as symmetry is used to define the trimer
-relax:fast \
-relax:jump_move true \
-edensity:mapfile 5wb0_phases_2mFo-DFc.ccp4 \
-edensity:mapreso 2.82 \ # x= INPUT map resolution 
-edensity:fastdens_wt 50.0 \ # y= density weight according to map resolution. 
-symmetry_definition 5wb0_M2.symm \ #symmetry definition file
-ignore_unrecognized_res \
-out::nstruct 5 \ #number of desired output structures
-ex1 -ex2 \
 

 

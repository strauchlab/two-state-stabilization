#!/bin/bash

$Rosetta/main/source/bin/relax.linuxgccrelease \
-database $Rosetta/main/database \
-in:file:s 5l1x_single_ABC_INPUT.pdb \ #the monomeric structure as symmetry is used to define the trimer
-relax:fast \
-relax:jump_move true \
-edensity:mapfile 5l1x_phases_2mFo-DFc.ccp4 \
-edensity:mapreso 3.35 \ # x= INPUT map resolution 
-edensity:fastdens_wt 50.0 \ # y= density weight according to map resolution. 
-symmetry_definition 5l1x_M2.symm \ #symmetry definition file
-ignore_unrecognized_res \
-out::nstruct 5 \ #number of desired output structures
-ex1 -ex2 \
 

 

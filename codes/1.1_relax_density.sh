#!/bin/bash

$Rosetta/main/source/bin/relax.linuxgccrelease \
-database $Rosetta/main/database \
-in:file:s INPUT_Monomer.pdb \ #the monomeric structure as symmetry is used to define the trimer
-relax:fast \
-relax:jump_move true \
-edensity:mapfile INPUT_map.ccp4 \
-edensity:mapreso x \ # x= INPUT map resolution 
-edensity:fastdens_wt y \ # y= density weight according to map resolution. 
-symmetry_definition INPUT_symmetry_file.symm \ #symmetry definition file
-ignore_unrecognized_res \
-out::nstruct 100 \ #number of desired output structures
-ex1 -ex2 \
 

 

#!/bin/bash

$Rosetta/main/source/bin/relax.linuxgccrelease \
-database $Rosetta/main/database \
-in:file:s 5l1x_INPUT.pdb \
-relax:fast \
-relax:jump_move true \
-edensity:mapfile 5l1x_phases_2mFo-DFc.ccp4 \
-edensity:mapreso 3.35 \
-edensity:fastdens_wt 50.0 \
-symmetry_definition 5l1x_M2.symm \
-ignore_unrecognized_res \
-out::nstruct 100 \
-ex1 -ex2 \
 

 

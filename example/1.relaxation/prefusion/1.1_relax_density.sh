#!/bin/bash

$Rosetta/main/source/bin/relax.linuxgccrelease \
-database $Rosetta/main/database \
-in:file:s 5wb0_INPUT.pdb \
-relax:fast \
-relax:jump_move true \
-edensity:mapfile 5wb0_phases_2mFo-DFc.ccp4 \
-edensity:mapreso 2.82 \
-edensity:fastdens_wt 50.0 \
-symmetry_definition 5wb0_M2.symm \
-ignore_unrecognized_res \
-out::nstruct 100 \
-ex1 -ex2 \
 

 

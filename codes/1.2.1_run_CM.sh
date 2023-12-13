#!/bin/bash

$Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease \
-database $Rosetta/main/database \
-in:file:fasta INPUT.fasta \
-parser:protocol 1.2.0_rosetta_comparative_modeling.xml \ 
-relax:dualspace \
-relax:jump_move true \
-edensity:mapfile INPUT_map.map \ #cryoEM map
-edensity:mapreso x \ # x= INPUT map resolution 
-edensity::cryoem_scatterers \
-beta \
-hybridize:stage1_probability 1.0 \
-ignore_unrecognized_res \
-out:nstruct 100 \ #number of desired output structures

 
 

#!/bin/bash

$Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease \
-database $Rosetta/main/database \
-in:file: INPUT_Monomer.pdb \ #the monomeric structure as symmetry is used to define the trimer
-parser: protocol 1.2.2_second_relax_density.xml \
-parser:script_vars map=INPUT.map symm=INPUT_symmetry_file.symm \
-edensity:mapreso x \ # x= INPUT map resolution 
-edensity:cryoem_scatterers \
-beta \
-ignore_unrecognized_res \
-crystal_refine \
-out:nstruct 100 \ #number of desired output structures


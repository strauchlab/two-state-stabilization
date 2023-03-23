#!/bin/bash

$Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease \
-database $Rosetta/main/database \
-in:file:s INPUT_Monomer.pdb \ #the monomeric structure as symmetry is used to define the trimer
-parser:protocol 6.0_combinatorial_design.xml \
-parser:script_vars symm_file=INPUT_symmetry.symm resfile=resfile_pre pssm=pre.pssm \
@6.1_flags -overwrite

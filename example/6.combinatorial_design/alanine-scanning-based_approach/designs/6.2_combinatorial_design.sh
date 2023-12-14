#!/bin/bash

$Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease \
-database $Rosetta/main/database \
-in:file:s 5wb0_INPUT_0001_AB_0007_ABC_A.pdb \
-parser:protocol 6.0_combinatorial_design.xml \
-parser:script_vars symm_file=5wb0_M2.symm resfile=resfile_pre pssm=pre.pssm \
@6.1_flags -overwrite

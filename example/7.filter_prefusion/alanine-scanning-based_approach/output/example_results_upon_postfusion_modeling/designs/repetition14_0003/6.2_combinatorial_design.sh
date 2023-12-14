#!/bin/bash

$Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease \
-database $Rosetta/main/database \
-in:file:s 5l1x_INPUT_0005_AB_0007_ABC_A.pdb \
-parser:protocol 6.0_combinatorial_design.xml \
-parser:script_vars symm_file=5l1x_M2.symm resfile=resfile_post pssm=post.pssm \
@6.1_flags -overwrite





#!/bin/bash

$Rosetta/main/source/bin/cartesian_ddg.static.linuxgccrelease \
-database $Rosetta/main/database \
-s 5wb0_INPUT_0001_AB_0007.pdb \
-ddg:iterations 3 \
-ddg::cartesian \
-ddg::dump_pdbs false \
-ddg:bbnbrs 1 \
-beta_cart \
-fa_max_dis 9.0 \
-ddg:mut_file ./mut_files/mut_file1 


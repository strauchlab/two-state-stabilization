#!/bin/bash

$Rosetta/main/source/bin/cartesian_ddg.static.linuxgccrelease \
-database $Rosetta/main/database \
-s INPUT.pdb \ #trimeric structure as cartesian_ddg is not aware of symmetry
-ddg:iterations 3 \
-ddg::cartesian \
-ddg::dump_pdbs false \
-ddg:bbnbrs 1 \
-beta_cart \
-fa_max_dis 9.0 \
-ddg:mut_file mut_file \ #mut_file contains the position to be mutated by alanine. You need to iterate over all mut_files/positions



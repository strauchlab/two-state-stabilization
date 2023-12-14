#!/bin/bash

$Rosetta/main/source/bin/relax.linuxgccrelease \
-database $Rosetta/main/database \
-s 5l1x_INPUT_0005_AB.pdb \
-use_input_sc \
-constrain_relax_to_start_coords \
-ignore_unrecognized_res \
-nstruct 15 \
-relax:cartesian -score:weights ref2015_cart \
-relax:min_type lbfgs_armijo_nonmonotone \
-relax:script 2.0_cart2.script \
-symmetry_definition 5l1x_M2.symm 

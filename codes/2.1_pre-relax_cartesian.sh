#!/bin/bash

$Rosetta/main/source/bin/relax.linuxgccrelease \
-database $Rosetta/main/database \
-s INPUT_Monomer.pdb \ #the monomeric structure as symmetry is used to define the trimer. This should be the best structure obtained after the density relaxation
-use_input_sc \
-constrain_relax_to_start_coords \
-ignore_unrecognized_res \
-nstruct 15 \ #number of desired output structures
-relax:cartesian -score:weights ref2015_cart \
-relax:min_type lbfgs_armijo_nonmonotone \
-relax:script 2.0_cart2.script \
-symmetry_definition INPUT_symmetry_file.symm \ #symmetry definition file

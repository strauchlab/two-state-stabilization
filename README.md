# Two-state-stabilization
## Requirements:

1. Rosetta version: 2020.10.post.dev+12.master.c7b9c3e c7b9c3e4aeb1febab211d63da2914b119622e69b  
   Instructions on how to install Rosetta can be found here: https://new.rosettacommons.org/demos/latest/tutorials/install_build/install_build
3. PyRosetta-4 2019 [Rosetta PyRosetta4.conda.linux.CentOS.python37.Release 2019.47+release.3d995b15922374726493453159734beedd7e28be 2019-11-20T17:52:20]  
   Instructions on how to install PyRosetta can be found here: https://www.pyrosetta.org/downloads
5. Python3.7
6. Python packages: 
	- os
 	- pandas 
 	- numpy
 	- subprocess
 	- argparse
 	- pyrosetta 
 	- rosetta
 	- statistics
 	- math
 	- re

> **NOTE:** All scripts were tested in Ubuntu 18.04.5 LTS

## Step-by-step scripts
--------------------- 
To avoid unexpected results, the input pdbs must be cleaned before starting the analysis. For this, use the script $Rosetta/tools/protein_tools/scripts/clean_pdb.py

### 1. Structural relaxation of both pre- and postfusion structures guided by density data: 

Script used for RSV F and hMPV F:

	-1.1_relax_density.sh

Scripts used for SARS-CoV-2 S:

	Comparative modeling:
	-1.2.0_rosetta_comparative_modeling.xml
	-1.2.1_run_CM.sh

	Second relaxation:
	-1.2.2_second_relax_density.xml
	-1.2.3_run_second_RD.sh

Symmetry definition files are generated using the command line: "perl $Rosetta/main/source/src/apps/public/symmetry/make_symmdef_file.pl -m NCS -a A -i B -p INPUT.pdb > symmetry.symm"  
Before using the symmetry file, it needs to be modified as described at https://faculty.washington.edu/dimaio/files/rosetta_density_tutorial_aug18.pdf

### 2. Alanine scanning with standard “Cartesian ddg”:

Protocol as described at https://www.rosettacommons.org/docs/latest/cartesian-ddG:

> **IMPORTANT:** If the input PDB has segments where the sequence is discontinuos within the same protomer, each segment must be labeled with a different chain ID. This is because the cartesian ddg application can generate artificial bonds for disconnected regions. 
In addition, make sure the residue numbering for each protomer starts from 1 (i.e., 1,2,3..,lenght of protomer), independently of the chain ID. This is because all our analysis are done with the rosetta numbering corresponding to the first monomer. 
Example:    

					      Protomer 1
			 --------------------- chain A ---------------------
	original PDB#    25 26 27 28 29 30 31 32 --- 50 51 52 53 54 55 56 57

					     Protomer 1
			 ------- chain A -------    ------- chain B -------
	modified PDB#     1  2  3  4  5  6  7  8 --- 9 10 11 12 13 14 15 16


Scripts to perform alanine scanning:

	Pre-relaxation in cartesian space:
	-2.0_cart2.script
	-2.1_pre-relax_cartesian.sh

	Write mut_files:
	-2.2_mut_file_alanine.py [-h] pdb ch 
This script takes as input the pdb to be analyzed and the chain(s) ID of one protomer. As output, it creates a folder with the mut_files needed to carry out alanine scanning. Note that mut_files are written in rosetta numbering. 

	Alanine scanning execution:
	-2.3_alanine_scanning.sh

Alanine scanning needs to be carried out on every residue of both pre- and postfusion structures. To ensure a robust analysis, we recommend performing the alanine scanning in at least two independent runs. On step 3 the results are combined and outliers are filtered out. 


### 3. Comparison prefusion vs postfusion alanine ddg results and selection of target positions for all-amino acids scanning: 

	-3.1_comparison_alanines.py [-h] arg_file 
This script identifies target positions to redesign based on stabilization of the prefusion state over the postfusion state. It takes as input an external file containing all arguments needed for running the script. As output, it creates two folders, one for each state (prefusion and postfusion), containing the mut_files required for all-amino acids scanning at the identified positions. Note that mut_files are written in rosetta numbering. 

	-3.2_movable_regions.py [-h] arg_file
If alanine scanning is not enough to identify significant designable spots, all- amino acids scanning can be carried out on all regions undergoing drastic conformational changes. This script identifies those movable regions. It takes as input an external file containing all arguments needed for running the script. As output, it creates two folders, one for each state (prefusion and postfusion), containing the mut_files to perform all-amino acids scanning at highly movable regions. NOTE: Input PDBs should be aligned prior to running the script.


### 4. Perform all-amino acids substitutions with the mut_files generated in step 3. 
The script to perform the substitutions is the same as the alanine scanning (-2.2_alanine_scanning.sh). To ensure a robust analysis, we recommend performing the all-amino acids substitutions in at least two independent runs. On step 5 the results are combined and outliers are filtered out. 


### 5. Comparison of prefusion vs postfusion all-amino acids scanning ddg results and selection of positions for combinatorial design:
	
	-5_analysis_all_aa_substitutions.py [-h] arg_file
This script selects the positions and the substitutions to be combined based on favorable ddg results for the prefusion state and neutral or destabilizing results for the postfusion state. The scripts takes as input an external file containing all arguments needed for running it, and it outputs a folder called "combinatorial_design". This folder contains the PSSM-like file for prefusion and postfusion, a resfile file for redesigning the prefusion structure, and a control resfile (prefusion).


### 6. Combinatorial design on prefusion state. 

	-6.0_combinatorial_design.xml
	-6.1_flags
	-6.2_combinatorial_design.sh
	
These scripts need to be executed twice, once for generating the designed sequences and another time for running a control experiment. Since the selection process is done based on energetic differences, a control experiment is required to compare against a structure that has undergone the same relaxation as the designed sequence (but without sequence changes). For the control experiment the resfile used is called "resfile_control_pre". The PSSM file is the same as for the designed sequences.

Independent jobs can be run in parallel to increase the number of designs while reducing the waiting time. Remember to save the results in different folders to avoid overwriting previous results.


### 7. Filter out designed sequences that are repeated and identify sequences improving the prefusion energy compared to the control (wild-type sequence)

	-7_filter_prefusion.py

After finding sequences improving the prefusion energy, these sequences need to be modeled in the postfusion state to determine which conformation is more favored by the designed sequence. Therefore, this script first identifies the candidate sequences and then generates several folders with individual resfiles to model the corresponding sequence in the postfusion state. Additionally, a control postfusion resfile is generated to perform the control experiment as explained on step 6. These folders are found in folder called "postfusion_designs"
All pdbs with improved prefusion energy are transferred to a folder called "selection_pre_E" created inside the prefusion results folder.

### 8. Identify sequences where prefusion improved and postfusion got worse. 

	-8_selection_by_energy.py

This script identifies which of the sequences improving the prefusion energy do not favor the postfusion structure. These sequences would be the leading canditates. To summarize the results, this script generates two types of files, one containing all the mutatations present on each leading design, and another set of csv files with per-residue energy differences between the pre- and postfusion structures at each mutated position ( negative scores indicate that the substitution is more stable than the native sequence in the specified structure. Therefore, we are looking for mutations with negative scores in the prefusion state and positive scores in the postfusion state). All numbering refers to rosetta numbering for the first protomer of the prefusion structure (if it is not specified otherwise).
Since not all mutations in a sequence contribute equally to stabilize the prefusion structure or destabilize postfusion, it is recommended to select only mutations with the most significant effects. This selection can be guided by the output energetic differences, analysis  of formation or disruption of hydrogen bonds and salt bridges (which can be done with the function "energy_terms" available in this script), or by manual inspection of each mutation.  

After few mutations are selected, we recommend to repeat the design process (steps 6 to 8) to verify whether the selected mutations stabilize the prefusion structure over postfusion.



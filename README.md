# Two-state-stabilization
## Requirements:

1. Rosetta version: 2020.10.post.dev+12.master.c7b9c3e c7b9c3e4aeb1febab211d63da2914b119622e69b  
   Instructions on how to install Rosetta can be found [here](https://new.rosettacommons.org/demos/latest/tutorials/install_build/install_build).
3. PyRosetta-4 2019 [Rosetta PyRosetta4.conda.linux.CentOS.python37.Release 2019.47+release.3d995b15922374726493453159734beedd7e28be 2019-11-20T17:52:20]  
   Instructions on how to install PyRosetta can be found [here](https://www.pyrosetta.org/downloads).
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

> [!NOTE]
> All scripts were tested in Ubuntu 18.04.5 LTS

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

Refer to [this tutorial](https://faculty.washington.edu/dimaio/files/rosetta_density_tutorial_aug18.pdf) (Example 2B: Symmetric refinement into cryoEM density) for more information on symmetry files.

Once the relaxation process is completed, select the best structure based on Rosetta energy, Molprobity scores, and agreement with density data. This chosen structure will serve as the input for step (2).

### 2. Alanine scanning with standard “Cartesian ddg”:

Follow [this protocol](https://www.rosettacommons.org/docs/latest/cartesian-ddG) to execute the cartesian ddg application. 

> [!IMPORTANT]
> If the input PDB has segments where the sequence is discontinuous within the same protomer, each segment must be labeled with a different chain ID. The cartesian ddg application can generate artificial bonds for disconnected regions.  
In addition, make sure the residue numbering for each protomer starts from 1 (i.e., 1,2,3.., length of protomer), independently of the chain ID. All our analyses are done with the rosetta numbering corresponding to the first monomer. 
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
The -2.2_mut_file_alanine.py script takes as input the pdb to be analyzed and the chain(s) ID of one protomer. An example of the command line to run the script is found at */example/2.alanine_scanning/prefusion/scanning/2.2_arg_hmpv*

As output, the script creates a folder with the mut_files needed for alanine scanning. Note that mut_files are written in rosetta numbering.  

	Alanine scanning execution:
	-2.3_alanine_scanning.sh

Alanine scanning must be carried out on every residue of pre- and postfusion structures. We recommend performing the alanine scanning in at least two independent runs to ensure a robust analysis. In step 3, the results are combined, and outliers are filtered out.  


### 3. Comparison prefusion-vs-postfusion alanine ddg results and selection of target positions for all-amino acid scanning: 

	-3.1_comparison_alanines.py [-h] arg_file 
This script identifies target positions to redesign based on stabilizing the prefusion over the postfusion state. It takes as input an external file containing all arguments needed for running the script. An example of this external file is found at */example/3.comparison_ala_scan/3.1/3.1_arg_hmpv*. 

As output, the script creates two folders, one for each state (prefusion and postfusion), containing the mut_files required for all-amino acid scanning at the identified positions. Note that mut_files are written in rosetta numbering.

	-3.2_mobile_regions.py [-h] arg_file
If alanine scanning cannot identify significant designable spots, all-amino acid scanning can be carried out on all regions undergoing drastic conformational changes. This script identifies those mobile regions. It takes as input an external file containing all arguments needed for running the script. An example of this external file is found at */example/3.comparison_ala_scan/3.2/3.2_arg_hmpv*. 

As output, the script creates two folders, one for each state (prefusion and postfusion), containing the mut_files to perform all-amino acid scanning at highly mobile regions. 

> [!NOTE]
> Input PDBs should be aligned prior to running the script.


### 4. Perform all-amino acid substitutions with the mut_files generated in step 3. 
Use the same script as alanine scanning (2.2_alanine_scanning.sh). We recommend performing the all-amino acid substitutions in at least two independent runs to ensure a robust analysis. In step 5, the results are combined, and outliers are filtered out.


### 5. Comparison of prefusion-vs-postfusion all-amino acid scanning ddg results and selection of positions for combinatorial design:
	
	-5_analysis_all_aa_substitutions.py [-h] arg_file
This script selects positions and substitutions for combinatorial design based on favorable ddG results for the prefusion state and neutral or destabilizing results for the postfusion state. The script takes as input an external file containing all arguments needed for running it. An example of this external file is found at */example/5.comparison_all_amino_acids_scan/alanine-scanning-based_approach/5_arg_hmpv*.

The script outputs a "combinatorial_design" folder. This folder contains the PSSM-like file for prefusion and postfusion, a resfile file for redesigning the prefusion structure, and a control resfile (prefusion).  

Run this script independently for alanine-scanning-based and mobile-regions-based approaches.

### 6. Combinatorial design on prefusion state. 

	-6.0_combinatorial_design.xml
	-6.1_flags
	-6.2_combinatorial_design.sh
	
These scripts must be executed twice - once to generate the designed sequences and once to run a control experiment without sequence design. The control experiment is used to define the reference energy as the selection process is based on energetic differences. The control resfile is called "resfile_control_pre", and the PSSM file would be the same as for the designed sequences.

Independent jobs can be run in parallel to increase the number of designs while reducing the waiting time. Remember to save the results in different folders to avoid overwriting issues.


### 7. Filter out designed sequences that are repeated and identify sequences improving the prefusion energy compared to the control (wild-type sequence)

	-7_filter_prefusion.py [-h] arg_file

After finding sequences improving the prefusion energy, these sequences must be modeled in the postfusion state to determine which conformation is more favored by the designed sequence. Therefore, this script identifies the prefusion candidate sequences and generates individual resfiles to model the sequence in the postfusion state. A control postfusion resfile is also output to perform the control experiment, as explained in step (6). The resfiles are split into several folders within the "postfusion_designs" folder. The PSSM files to model the postfusion structures were obtained in step (5).

The script takes as input an external file containing all arguments needed for running it. An example of this external file is found at */example/7.filter_prefusion/alanine-scanning-based_approach/7_arg_hmpv*

> [!NOTE]
> All pdbs with improved prefusion energy are transferred to a "selection_pre_E" folder, created within the prefusion results folder.

### 8. Identify sequences where prefusion improved and postfusion got worse. 

	-8_selection_by_energy.py [-h] arg_file

This script identifies sequences improving prefusion energy but not favoring the postfusion structure. The script takes as input an external file containing all arguments needed for running it. An example of this external file is found at */example/8.selection/8_arg_hmpv*

To summarize the results, the script generates two types of files, one containing all mutations present on each leading design and CSV files with per-residue energy differences between the pre- and postfusion structures at each mutated position. Negative scores indicate that the substitution is more stable than the native sequence in the specified structure. Therefore, look for mutations with negative scores in the prefusion state and positive scores in the postfusion state. All numbering refers to rosetta numbering for the first protomer of the prefusion structure (if it is not specified otherwise).

Since not all mutations in a sequence contribute equally to stabilizing the prefusion structure or destabilizing postfusion, selecting only mutations with the most significant effects is recommended. This selection can be guided by the output energetic differences, analysis of formation or disruption of hydrogen bonds and salt bridges (which can be done with the function "energy_terms" available in this script), or by manual inspection of each mutation.  

After selecting a few mutations, we recommend repeating the design process using only the selected positions (steps 6 to 8) to verify whether the selected mutations stabilize the prefusion structure over postfusion.


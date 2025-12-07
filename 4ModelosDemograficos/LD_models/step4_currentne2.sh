#!/bin/bash

#⣠⣀⡠⠤⠂⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
#⢇⠀⠈⠒⠢⡀⠀⠁⠂⠄⡀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣤⡶⣿⡁
#⠀⠣⡀⠀⠀⠘⢟⡀⠀⠀⠀⠑⠒⠒⠐⠄⠄⠔⣊⠝⠳⡤⠄⠀
#⠀⠀⠈⠢⡀⠀⠀⠈⠁⠂⠤⣀⠀⠀⠀⢄⡢⠚⠁⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠑⠢⠤⠤⠤⠤⢮⠣⠐⠂⠁⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀

# Snigenda workshop
# Author: Kaden Winspear - Dec. 2025
# Usage: sbatch step3_goneRunner.sh [Population] Example: sbatch step3_goneRunner.sh ENP 
# Purpose: Runs CurrentNe2 on populations.

set -euo pipefail

# Load variables
pop=$1
inputdir="/Users/snigenda/Documents/workshop_2025/input_files/filtered/LD_ROH/LD_stuff/${pop}"
outputdir="/Users/snigenda/Documents/workshop_2025/results/gone_ne"
mkdir -p "${outputdir}"

# Go to output directory
cd "${outputdir}"

# Run CurrentNe2 on the input.
./Users/snigenda/Documents/workshop_2025/programs/currentNe2/currentne2 ${inputdir}/${pop}.ped 3 -o ${outputdir}/${pop}.currentNe.txt

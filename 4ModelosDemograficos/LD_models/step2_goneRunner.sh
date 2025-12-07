#!/bin/bash -l

#⣠⣀⡠⠤⠂⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
#⢇⠀⠈⠒⠢⡀⠀⠁⠂⠄⡀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣤⡶⣿⡁
#⠀⠣⡀⠀⠀⠘⢟⡀⠀⠀⠀⠑⠒⠒⠐⠄⠄⠔⣊⠝⠳⡤⠄⠀
#⠀⠀⠈⠢⡀⠀⠀⠈⠁⠂⠤⣀⠀⠀⠀⢄⡢⠚⠁⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠑⠢⠤⠤⠤⠤⢮⠣⠐⠂⠁⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀

# Snigenda workshop
# Author: Kaden Winspear - Dec. 2025
# Usage: sbatch step3_goneRunner.sh [Population] Example: sbatch step3_goneRunner.sh ENP 
# Purpose: Runs gone2 on populations.

set -euo pipefail

# Load variables
pop=$1
inputdir="/Users/snigenda/Documents/workshop_2025/input_files/filtered/LD_ROH/LD_stuff/${pop}"
outputdir="/Users/snigenda/Documents/workshop_2025/results/gone_ne"
mkdir -p "${outputdir}"

# Go to output directory
cd "${outputdir}"

# Run GONE2 on the input.
./Users/snigenda/Documents/workshop_2025/programs/GONE2/gone2 -t 6  -r 0.8 "${inputdir}/${pop}_Chr1-3_1M.ped" -o "${outputdir}"


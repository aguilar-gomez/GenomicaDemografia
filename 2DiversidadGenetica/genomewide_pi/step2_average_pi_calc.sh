#!/bin/bash


#⣠⣀⡠⠤⠂⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
#⢇⠀⠈⠒⠢⡀⠀⠁⠂⠄⡀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣤⡶⣿⡁
#⠀⠣⡀⠀⠀⠘⢟⡀⠀⠀⠀⠑⠒⠒⠐⠄⠄⠔⣊⠝⠳⡤⠄⠀
#⠀⠀⠈⠢⡀⠀⠀⠈⠁⠂⠤⣀⠀⠀⠀⢄⡢⠚⠁⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠑⠢⠤⠤⠤⠤⢮⠣⠐⠂⠁⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀

# Snigenda workshop
# Author:Kaden Winspear - Dec. 2025
# Purpose: This script takes in Pixy output & calculates average pi per pop.

mamba activate pixy

HOMEDIR=/Users/snigenda/Documents/workshop_2025
workdir="${HOMEDIR}/results/genomic_diversity/pixy"
input="${workdir}/Chrom19thru21_pi.txt"

cd ${workdir}

# extract header and add to separated txt files
head -n 1 ${input} > ENP_pi.txt
head -n 1 ${input} > GOC_pi.txt

# split by population (column 1)
awk 'NR>1 && $1=="ENP"' ${input} >> ENP_pi.txt
awk 'NR>1 && $1=="GOC"' ${input} >> GOC_pi.txt

# compute means of avg_pi (column 5)
awk 'NR>1 {sum+=$5;n++} END {print "ENP mean pi:",sum/n}' ENP_pi.txt > mean_pi.txt
awk 'NR>1 {sum+=$5;n++} END {print "GOC mean pi:",sum/n}' GOC_pi.txt >> mean_pi.txt

mamba deactivate


#!/bin/bash

# Author: Kaden Winspear
# Used to filter VCFs prior to running Pixy for simple filtering.
# sbatch [script]

#========WinBin===========
#        .
#       ":"
#     ___:____     |"\/"|
#   ,'        `.    \  /
#   |  O        \___/  |
# ~^~^~^~^~^~^~^~^~^~^~^~^~
#=======Genomics===========

mamba activate pixy

######## Directories ########

HOMEDIR=/Users/snigenda/Documents/workshop_2025
VCFDIR=$HOMEDIR/input_files/all_sites
OUTDIR=$HOMEDIR/input_files/filtered/pixy

mkdir -p ${OUTDIR}

# In vcf file
VCF="${VCFDIR}/ENP_GOC_chr19-21.vcf.gz"


# out file
outprefix="${OUTDIR}/Pixy_filtered"

# Run filtering
vcftools --gzvcf ${VCF} --remove-filtered-all --remove-indels --max-missing 0.8 --min-meanDP 20 --max-meanDP 200 --recode --out "${outprefix}" # takes like 14 minutes

# Compress
bgzip -c "${outprefix}.recode.vcf" > "${outprefix}.vcf.gz" # takes like 5 minutes

# Index
tabix -p vcf "${outprefix}.vcf.gz"

mamba deactivate


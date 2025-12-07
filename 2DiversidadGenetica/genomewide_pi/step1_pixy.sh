
#!/bin/bash


#⣠⣀⡠⠤⠂⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
#⢇⠀⠈⠒⠢⡀⠀⠁⠂⠄⡀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣤⡶⣿⡁
#⠀⠣⡀⠀⠀⠘⢟⡀⠀⠀⠀⠑⠒⠒⠐⠄⠄⠔⣊⠝⠳⡤⠄⠀
#⠀⠀⠈⠢⡀⠀⠀⠈⠁⠂⠤⣀⠀⠀⠀⢄⡢⠚⠁⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠑⠢⠤⠤⠤⠤⢮⠣⠐⠂⠁⠀⠀⠀⠀⠀⠀⠀⠀
#⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀

# Snigenda workshop
# Author:Kaden Winspear - Dec. 2025

mamba activate pixy

HOMEDIR=/Users/snigenda/Documents/workshop_2025
VCFDIR=${HOMEDIR}/input_files/filtered/pixy
VCF=${VCFDIR}/Passing_biallelic_all_19-21.vcf.gz
POP_FILE=${HOMEDIR}/input_files/config/popmap.txt
OUTPUT_DIR=${HOMEDIR}/results/genomic_diversity/pixy


# Run Pixy [10000 window]
pixy --stats pi --vcf $VCF --populations $POP_FILE --window_size 10000 --n_cores 6 --output_folder $OUTPUT_DIR --output_prefix "Chrom19thru21" # takes like 25 min

mamba deactivate


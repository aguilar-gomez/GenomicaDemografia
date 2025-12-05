#load necessary modules
. /u/local/Modules/default/init/modules.sh
#load plink #(v1.90b6.24)
module load plink/1.90b624
module load vcftools

###############################
#Apply thinning and compress the result
###############################
vcftools --vcf tmp_maf_filtered.recode.vcf --recode --thin 10000 --recode-INFO-all \
  --stdout | bgzip -c > briceiOutgroupAll_maf2_Thinned10Kb.vcf.gz
###############################
# Prepare data for ADMIXTURE
###############################
# Convert PLINK files to the correct format
plink --bfile finwhale \
      --make-bed \
      --allow-extra-chr \
      --out finwhale_admix


ADMIXTURE=/u/project/klohmuel/mkenfiel/OutgroupSpecies/RicesWhale/admixture/dist/admixture_linux-1.3.0/admixture
#run admixture for a K of 1-10, using cross-validation, with 10 threads
for K in 1 2 3 4 5 6 7 8 9 10; 
do $ADMIXTURE --cv -j10  briceiOutgroupAll_maf2_Thinned10Kb.bed $K | tee log${K}.out;
done

############ AFTER JOB ############
#Which K iteration is optimal according to ADMIXTURE ?
grep -h CV log*.out > log.errors.txt
#grep -h CV subsetLog*.out > subsetLog.errors.txt

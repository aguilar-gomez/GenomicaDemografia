#load necessary modules
. /u/local/Modules/default/init/modules.sh
#load plink #(v1.90b6.24)
module load plink/1.90b624
module load vcftools

###############################
# Apply thinning and compress the result
###############################

vcftools \
  --gzvcf subsetfinWhale_PASS_VAR_biallelicSNPs.vcf.gz \
  --thin 10000 \
  --recode \
  --recode-INFO-all \
  --stdout | bgzip -c > subsetfinWhale_PASS_VAR_biallelicSNPs_Thinned10Kb.vcf.gz

# Index
bcftools index subsetfinWhale_PASS_VAR_biallelicSNPs_Thinned10Kb.vcf.gz

###############################
# Prepare data for ADMIXTURE
###############################

# Convert to PLINK bed files
plink --vcf subsetfinWhale_PASS_VAR_biallelicSNPs_Thinned10Kb.vcf.gz \
      --make-bed \
      --allow-extra-chr \
      --out finwhale_admix \
      --double-id \
      --keep-allele-order


###############################
# Run ADMIXTURE K=1–10
###############################

ADMIXTURE=/u/project/klohmuel/mkenfiel/OutgroupSpecies/RicesWhale/admixture/dist/admixture_linux-1.3.0/admixture

for K in 1 2 3 4 5 6 7 8 9 10
do
    $ADMIXTURE --cv -j10 finwhale_admix.bed $K | tee log${K}.out
done

############ AFTER JOB ############
# Which K iteration is optimal according to ADMIXTURE?
grep -h CV log*.out > log.errors.txt

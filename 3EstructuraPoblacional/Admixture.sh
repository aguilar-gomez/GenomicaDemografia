#load necessary modules
. /u/local/Modules/default/init/modules.sh
#load plink #(v1.90b6.24)
module load plink/1.90b624

#use plink to convert vcf directly to bed format:
plink --vcf briceiOutgroupAll_maf2_Thinned10Kb.vcf.gz --double-id --make-bed --out briceiOutgroupAll_maf2_Thinned10Kb

ADMIXTURE=/u/project/klohmuel/mkenfiel/OutgroupSpecies/RicesWhale/admixture/dist/admixture_linux-1.3.0/admixture
#run admixture for a K of 1-10, using cross-validation, with 10 threads
for K in 1 2 3 4 5 6 7 8 9 10; 
do $ADMIXTURE --cv -j10  briceiOutgroupAll_maf2_Thinned10Kb.bed $K | tee log${K}.out;
done

############ AFTER JOB ############
#Which K iteration is optimal according to ADMIXTURE ?
grep -h CV log*.out > log.errors.txt
#grep -h CV subsetLog*.out > subsetLog.errors.txt

#Keep only Variants that pass filters
VCF=subset_ENP_GOC_merged.vcf.gz
bcftools view -f PASS -e 'INFO/VariantType="NO_VARIATION"' $VCF -Oz -o subsetfinWhale_PASS_VAR.vcf.gz

#Count variants
bcftools view $VCF | grep -v '^#' | wc -l

bcftools view subsetfinWhale_PASS_VAR.vcf.gz | grep -v '^#' | wc -l

#Keep only biallelic SNPs

#RUN PCA on plink

#Prepare data for admixture


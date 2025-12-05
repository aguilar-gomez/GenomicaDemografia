###############################
# Keep only variants that pass filters
###############################
VCF=subset_ENP_GOC_merged.vcf.gz

bcftools view -f PASS -e 'INFO/VariantType="NO_VARIATION"' $VCF \
    -Oz -o subsetfinWhale_PASS_VAR.vcf.gz
###############################
# Count variants (before and after filtering)
###############################
echo "Original:"
bcftools view $VCF | grep -v '^#' | wc -l

echo "After filtering:"
bcftools view subsetfinWhale_PASS_VAR.vcf.gz | grep -v '^#' | wc -l

###############################
# Keep only biallelic SNPs
###############################
bcftools view -m2 -M2 -v snps subsetfinWhale_PASS_VAR.vcf.gz \
    -Oz -o subsetfinWhale_PASS_VAR_biallelicSNPs.vcf.gz

bcftools index subsetfinWhale_PASS_VAR_biallelicSNPs.vcf.gz




###############################
# Convert to PLINK and run PCA
###############################
plink --vcf subsetfinWhale_PASS_VAR_biallelicSNPs.vcf.gz \
      --make-bed \
      --out finwhale \
      --allow-extra-chr \
      --double-id \
      --keep-allele-order

plink --bfile finwhale \
      --pca \
      --maf 0.05 \
      --allow-extra-chr \
      --out finwhale_pca


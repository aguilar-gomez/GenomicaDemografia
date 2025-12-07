###############################
# Convert to PLINK and run PCA
###############################
plink --vcf Passing_biallelic_var_19-21.vcf.gz \
      --make-bed \
      --out finwhale \
      --allow-extra-chr \
      --double-id \
      --keep-allele-order

#Total genotyping rate is 0.992468.
#1,325,140 variants and 10 people pass filters and QC.

plink --bfile finwhale \
      --pca \
      --maf 0.05 \
      --allow-extra-chr \
      --out finwhale_pca

#1,141,686 variants removed due to minor allele threshold(s)
#183,454 variants and 10 people pass filters and QC.

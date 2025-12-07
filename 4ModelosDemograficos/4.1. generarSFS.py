pip3 install vcfpy --user
#########################################
import vcfpy
import numpy as np
import dadi

# ---- Define populations ----
pop1 = ["ENPWA14", "ENPCA08", "ENPBC18", "ENPAK29", "ENPBC16"]
pop2 = ["GOC091", "GOC010", "GOC082", "GOC050", "GOC038"]

vcf_path = "Passing_biallelic_var_19-21.vcf.gz"

# ---- Load VCF ----
reader = vcfpy.Reader.from_path(vcf_path)

# Check that samples exist
vcf_samples = reader.header.samples.names
for s in pop1 + pop2:
    if s not in vcf_samples:
        raise ValueError("Sample {} not found in VCF".format(s))

# ---- Prepare SFS arrays ----
n1 = len(pop1)
n2 = len(pop2)
sfs = np.zeros((n1+1, n2+1), dtype=int)

# ---- Iterate over SNPs ----
for record in reader:
    # Skip non-biallelic
    if len(record.ALT) != 1:
        continue

    # Count alt alleles in each population
    ac1 = 0
    ac2 = 0

    for s in pop1:
        gt = record.call_for_sample[s].gt_bases
        if gt:
            alleles = gt.split("/")
            ac1 += alleles.count(record.ALT[0].value)

    for s in pop2:
        gt = record.call_for_sample[s].gt_bases
        if gt:
            alleles = gt.split("/")
            ac2 += alleles.count(record.ALT[0].value)

    # Fill SFS cell
    if 0 <= ac1 <= n1*2 and 0 <= ac2 <= n2*2:
        # Convert from chromosome counts (0..2n) to frequency bins (0..n)
        # dadi expects haploid counts → bin by haploid copies
        sfs[ac1, ac2] += 1

# ---- Create dadi Spectrum ----
fs = dadi.Spectrum(sfs, mask_corners=True)

# ---- Save ----
fs.to_file("SFS_2D.fs")
print("Saved 2D SFS → SFS_2D.fs")


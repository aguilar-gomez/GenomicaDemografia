#!/usr/bin/env python3

import allel
import numpy as np
import dadi

# -----------------------------
# INPUT
# -----------------------------
vcf_file = "Passing_biallelic_var_19-21.vcf.gz"
pop1 = ["ENPWA14", "ENPCA08", "ENPBC18" ,"ENPAK29", "ENPBC16" ]     # Cambia a tus IDs
pop2 = ["GOC091",  "GOC010",  "GOC082" , "GOC050" , "GOC038" ]     # Para 2D SFS (opcional)

# -----------------------------
# LOAD VCF
# -----------------------------
callset = allel.read_vcf(vcf_file, fields=["samples", "calldata/GT"])
samples = callset["samples"]
gt = allel.GenotypeArray(callset["calldata/GT"])

# -----------------------------
# POPULATION INDEXING
# -----------------------------
# Función para encontrar índices de individuos
def get_idx(pop_list):
    return [np.where(samples == i)[0][0] for i in pop_list]

idx1 = get_idx(pop1)
idx2 = get_idx(pop2)

# -----------------------------
# 1D SFS
# -----------------------------
ac1 = gt[:, idx1].count_alleles()
# Solo referencia/alternativa (bialélico)
ac1 = ac1[:, :2]

# Conteo de alelo alternativo
allele_counts1 = ac1[:, 1]
n_chrom = len(idx1) * 2

# Bins del SFS
sfs1d = np.histogram(allele_counts1, bins=np.arange(n_chrom+2))[0]

# Convertir a dadi Spectrum
spectrum1d = dadi.Spectrum(sfs1d)

spectrum1d.to_file("SFS_1D.fs")
print("Archivo generado: SFS_1D.fs")


# -----------------------------
# 2D SFS  (si quieres 2 poblaciones)
# -----------------------------
ac1 = gt[:, idx1].count_alleles()[:, 1]
ac2 = gt[:, idx2].count_alleles()[:, 1]

n1 = len(idx1)*2
n2 = len(idx2)*2

sfs2d = np.histogram2d(
    ac1,
    ac2,
    bins=[np.arange(n1+2), np.arange(n2+2)]
)[0]

spectrum2d = dadi.Spectrum(sfs2d)

spectrum2d.to_file("SFS_2D.fs")
print("Archivo generado: SFS_2D.fs")

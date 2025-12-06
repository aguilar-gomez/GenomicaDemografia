#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Dadi demo script
----------------
Genera SFS falsos (1D y 2D), define modelos demográficos 
y ejecuta optimización de parámetros.

Requiere:
    pip install dadi
"""

import dadi
import numpy as np
import matplotlib.pyplot as plt


# -------------------------------------------------------------------
# 1. GENERAR SFS FALSOS (SIMULADOS)
# -------------------------------------------------------------------
# 1D SFS con 20 cromosomas muestreados (10 individuos diploides)
ns_1d = [20]
sfs1d = dadi.Spectrum(np.random.randint(1, 50, size=(21,)))
sfs1d = sfs1d.fold()  # simular unifariante
sfs1d.to_file("fake_1D.fs")

print("Se generó fake_1D.fs")

# 2D SFS con dos poblaciones: 20 y 18 cromosomas
ns_2d = [20, 18]
sfs2d = dadi.Spectrum(np.random.randint(1, 20, size=(21, 19)))
sfs2d = sfs2d.fold()
sfs2d.to_file("fake_2D.fs")

print("Se generó fake_2D.fs")


# -------------------------------------------------------------------
# 2. DEFINIR MODELOS DEMOGRÁFICOS 1D
# -------------------------------------------------------------------

def model_1epoch(params, ns, pts):
    """
    Modelo 1 época: Ne constante
    params = [Nu]
    """
    Nu = params[0]
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, T=1, nu=Nu)
    return dadi.Spectrum.from_phi(phi, ns, (xx,))

def model_2epoch(params, ns, pts):
    """
    Modelo 2 épocas: cambio de tamaño hace T
    params = [Nu, T]
    """
    Nu, T = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, T, nu=Nu)
    return dadi.Spectrum.from_phi(phi, ns, (xx,))

def model_3epoch(params, ns, pts):
    """
    Modelo 3 épocas: cuello de botella + recuperación
    params = [NuB, TB, NuF, TF]
    """
    NuB, TB, NuF, TF = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, TB, nu=NuB)
    phi = dadi.Integration.one_pop(phi, xx, TF, nu=NuF)
    return dadi.Spectrum.from_phi(phi, ns, (xx,))


# -------------------------------------------------------------------
# 3. DEFINIR MODELOS DEMOGRÁFICOS 2D
# -------------------------------------------------------------------

def split_no_mig(params, ns, pts):
    """
    Split sin migración
    params=[T]
    """
    T = params[0]
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(phi)
    phi = dadi.Integration.two_pops(phi, xx, T, 1, 1, m12=0, m21=0)
    return dadi.Spectrum.from_phi(phi, ns, (xx, xx))

def split_mig_sym(params, ns, pts):
    """
    Split con migración simétrica
    params=[T, m]
    """
    T, m = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(phi)
    phi = dadi.Integration.two_pops(phi, xx, T, 1, 1, m12=m, m21=m)
    return dadi.Spectrum.from_phi(phi, ns, (xx, xx))

def IM_model(params, ns, pts):
    """
    Isolation-with-Migration
    params=[T, nu1, nu2, m12, m21]
    """
    T, nu1, nu2, m12, m21 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(phi)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1, nu2, m12, m21)
    return dadi.Spectrum.from_phi(phi, ns, (xx, xx))

def secondary_contact(params, ns, pts):
    """
    Modelo de contacto secundario
    params=[Tiso, Tmig, m12, m21]
    """
    Tiso, Tmig, m12, m21 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(phi)
    phi = dadi.Integration.two_pops(phi, xx, Tiso, 1, 1, 0, 0)
    phi = dadi.Integration.two_pops(phi, xx, Tmig, 1, 1, m12, m21)
    return dadi.Spectrum.from_phi(phi, ns, (xx, xx))


# -------------------------------------------------------------------
# 4. OPTIMIZAR UN EJEMPLO 1D
# -------------------------------------------------------------------

# Cargar SFS falso
data_1d = dadi.Spectrum.from_file("fake_1D.fs")

# Grid points
pts_l = [40, 50, 60]

# Escoger modelo (ejemplo: 2 épocas)
func = model_2epoch
func_ex = dadi.Numerics.make_extrap_log_func(func)

# Parámetros iniciales, límites:
p0 = [0.5, 0.1]
lower = [1e-3, 1e-4]
upper = [10, 10]

print("\nOptimizando 1D modelo 2 épocas...")
popt = dadi.Inference.optimize_log(p0, data_1d, func_ex, pts_l,
                                   lower_bound=lower, upper_bound=upper,
                                   verbose=len(p0))

model = func_ex(popt, data_1d.sample_sizes, pts_l)

ll = dadi.Inference.ll_multinom(model, data_1d)
theta = dadi.Inference.optimal_sfs_scaling(model, data_1d)

print("\nRESULTADOS 1D:")
print("Parametros óptimos:", popt)
print("Log-likelihood:", ll)
print("Theta:", theta)

print("\nGraficando resultados 1D...\n")

# Plot SFS observado
plt.figure(figsize=(6,4))
dadi.Plotting.plot_1d_fs(data_1d)
plt.title("SFS observado")
plt.tight_layout()
plt.show()

# Comparación Modelo vs Datos (Poisson)
plt.figure(figsize=(6,4))
dadi.Plotting.plot_1d_comp_Poisson(model, data_1d)
plt.title("Modelo vs Datos (Poisson)")
plt.tight_layout()
plt.show()

# Comparación Modelo vs Datos (Multinomial)
plt.figure(figsize=(6,4))
dadi.Plotting.plot_1d_comp_multinom(model, data_1d)
plt.title("Modelo vs Datos (Multinomial)")
plt.tight_layout()
plt.show()


# -------------------------------------------------------------------
# 5. OPTIMIZAR UN EJEMPLO 2D
# -------------------------------------------------------------------

# Cargar SFS 2D
data_2d = dadi.Spectrum.from_file("fake_2D.fs")

# Definir modelo (split simétrico)
func2 = split_mig_sym
func2_ex = dadi.Numerics.make_extrap_log_func(func2)

# Grid sizes para extrapolación
pts_l = [40, 50, 60]

# Parámetros iniciales y límites
p0 = [0.5, 0.1]   # [T, m]
lower = [1e-3, 0]
upper = [10, 10]

print("\nOptimizando 2D modelo split + migración...")

popt2 = dadi.Inference.optimize_log(
    p0, data_2d, func2_ex, pts_l,
    lower_bound=lower, upper_bound=upper,
    verbose=len(p0)
)

# Modelo optimizado
model2 = func2_ex(popt2, data_2d.sample_sizes, pts_l)

# Likelihood y theta
ll2 = dadi.Inference.ll_multinom(model2, data_2d)
theta2 = dadi.Inference.optimal_sfs_scaling(model2, data_2d)

print("\nRESULTADOS 2D:")
print("Parámetros óptimos:", popt2)
print("Log-likelihood:", ll2)
print("Theta:", theta2)

# ------------------------------------------------------------
# Gráficas con funciones disponibles en tu versión de dadi
# ------------------------------------------------------------

# 1) SFS observado vs modelo
plt.figure(figsize=(6,5))
dadi.Plotting.plot_2d_comp_Poisson(model2, data_2d)
plt.title("Comparación 2D modelo vs datos (Poisson)")
plt.tight_layout()
plt.savefig("2D_model_fit.png", dpi=300)
plt.close()

# 2) SFS observado solamente
plt.figure(figsize=(6,5))
dadi.Plotting.plot_single_2d_sfs(data_2d)
plt.title("SFS observado (2D)")
plt.tight_layout()
plt.savefig("2D_SFS_observed.png", dpi=300)
plt.close()

print("\nListo: análisis 2D optimizado y gráficas guardadas.")

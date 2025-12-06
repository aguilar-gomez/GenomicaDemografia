#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_dadi_full_py36.py
Versión limpia y compatible con Python 3.6 y dadi antiguo.
- Genera SFS falsos (1D y 2D)
- Define modelos 1D (1,2,3 epochs)
- Define modelos 2D (split, split+mig sim, IM, secondary contact)
- Ajusta parámetros con dadi.Inference.optimize_log
- Guarda plots en PNG en la carpeta dadi_demo_out
Requiere: dadi, numpy, matplotlib
"""

from __future__ import print_function
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")   # backend no interactivo para clusters
import matplotlib.pyplot as plt

import dadi

# ---------------------------
# Config
# ---------------------------
OUTDIR = "dadi_demo_out"
os.makedirs(OUTDIR, exist_ok=True)

# Grid points para extrapolación
pts_l = [40, 50, 60]

# reproducibilidad
np.random.seed(42)

# ---------------------------
# 1) GENERAR SFS FALSOS (compatible con dadi antiguo)
# ---------------------------
# 1D: tamaño de muestra = 20 cromosomas (10 individuos diploides)
ns_1d = [20]
arr1d = np.abs(np.random.poisson(5, size=(ns_1d[0] + 1,))).astype(float)
sfs1d = dadi.Spectrum(arr1d)   # constructor sin argumento 'folded'
sfs1d = sfs1d.fold()           # doblar aquí
sfs1d.to_file(os.path.join(OUTDIR, "fake_1D.fs"))

# 2D: tamaños de muestra 20 y 18 cromosomas
ns_2d = [20, 18]
arr2d = np.abs(np.random.poisson(3, size=(ns_2d[0] + 1, ns_2d[1] + 1))).astype(float)
sfs2d = dadi.Spectrum(arr2d)
sfs2d = sfs2d.fold()
sfs2d.to_file(os.path.join(OUTDIR, "fake_2D.fs"))

print("SFS falsos generados en:", OUTDIR)

# ---------------------------
# 2) MODELOS 1D
# ---------------------------

def model_1epoch(params, ns, pts):
    # params = [Nu]
    Nu = params[0]
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, T=1.0, nu=Nu)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

def model_2epoch(params, ns, pts):
    # params = [Nu, T]
    Nu, T = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, T=T, nu=Nu)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

def model_3epoch(params, ns, pts):
    # params = [NuB, TB, NuF, TF]
    NuB, TB, NuF, TF = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, T=TB, nu=NuB)
    phi = dadi.Integration.one_pop(phi, xx, T=TF, nu=NuF)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs

# ---------------------------
# 3) MODELOS 2D (compatibles con versiones antiguas)
# ---------------------------

def split_no_mig(params, ns, pts):
    # params = [T]
    T = params[0]
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    # phi_1D_to_2D requiere (phi, xx) en versiones antiguas
    phi = dadi.PhiManip.phi_1D_to_2D(phi, xx)
    # Algunas versiones tienen split_2D; si existe, úsala
    if hasattr(dadi.PhiManip, "split_2D"):
        phi = dadi.PhiManip.split_2D(phi, xx)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1=1.0, nu2=1.0, m12=0.0, m21=0.0)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

def split_mig_sym(params, ns, pts):
    # params = [T, m]
    T, m = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(phi, xx)
    if hasattr(dadi.PhiManip, "split_2D"):
        phi = dadi.PhiManip.split_2D(phi, xx)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1=1.0, nu2=1.0, m12=m, m21=m)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

def IM_model(params, ns, pts):
    # params = [T, nu1, nu2, m12, m21]
    T, nu1, nu2, m12, m21 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(phi, xx)
    if hasattr(dadi.PhiManip, "split_2D"):
        phi = dadi.PhiManip.split_2D(phi, xx)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1=nu1, nu2=nu2, m12=m12, m21=m21)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

def secondary_contact(params, ns, pts):
    # params = [Tiso, Tmig, m12, m21]
    Tiso, Tmig, m12, m21 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(phi, xx)
    if hasattr(dadi.PhiManip, "split_2D"):
        phi = dadi.PhiManip.split_2D(phi, xx)
    # aislamiento
    phi = dadi.Integration.two_pops(phi, xx, Tiso, nu1=1.0, nu2=1.0, m12=0.0, m21=0.0)
    # contacto secundario
    phi = dadi.Integration.two_pops(phi, xx, Tmig, nu1=1.0, nu2=1.0, m12=m12, m21=m21)
    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

# ---------------------------
# 4) AJUSTE 1D (ejemplo)
# ---------------------------
print("\n=== AJUSTE 1D: modelo de 2 épocas (ejemplo) ===\n")
data_1d = dadi.Spectrum.from_file(os.path.join(OUTDIR, "fake_1D.fs"))
func_1d = model_2epoch
func_1d_ex = dadi.Numerics.make_extrap_log_func(func_1d)

# params iniciales [Nu, T]
p0 = [0.5, 0.2]
lower = [1e-3, 1e-4]
upper = [20, 20]

popt1 = dadi.Inference.optimize_log(p0, data_1d, func_1d_ex, pts_l,
                                    lower_bound=lower, upper_bound=upper,
                                    verbose=len(p0))

model1 = func_1d_ex(popt1, data_1d.sample_sizes, pts_l)
ll1 = dadi.Inference.ll_multinom(model1, data_1d)
theta1 = dadi.Inference.optimal_sfs_scaling(model1, data_1d)

print("Parámetros óptimos (1D):", popt1)
print("Log-likelihood (1D):", ll1)
print("Theta (1D):", theta1)

# Guardar plots 1D
plt.figure(figsize=(6, 4))
if hasattr(dadi.Plotting, "plot_1d_fs"):
    dadi.Plotting.plot_1d_fs(data_1d)
else:
    print("plot_1d_fs no disponible")
plt.title("SFS observado (1D)")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "1D_SFS_observed.png"), dpi=300)
plt.close()

plt.figure(figsize=(6, 4))
if hasattr(dadi.Plotting, "plot_1d_comp_Poisson"):
    dadi.Plotting.plot_1d_comp_Poisson(model1, data_1d)
elif hasattr(dadi.Plotting, "plot_1d_comp_multinom"):
    dadi.Plotting.plot_1d_comp_multinom(model1, data_1d)
else:
    print("No hay función de comparación 1D disponible")
plt.title("Comparación: modelo vs datos (1D)")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "1D_model_comp.png"), dpi=300)
plt.close()

# ---------------------------
# 5) AJUSTE 2D (ejemplo)
# ---------------------------
print("\n=== AJUSTE 2D: modelo split + migración simétrica (ejemplo) ===\n")
data_2d = dadi.Spectrum.from_file(os.path.join(OUTDIR, "fake_2D.fs"))

func2 = split_mig_sym
func2_ex = dadi.Numerics.make_extrap_log_func(func2)

# params iniciales [T, m]
p0_2 = [0.5, 0.2]
lower2 = [1e-4, 0.0]
upper2 = [20, 20]

# Intento de optimización con manejo de errores:
try:
    popt2 = dadi.Inference.optimize_log(p0_2, data_2d, func2_ex, pts_l,
                                        lower_bound=lower2, upper_bound=upper2,
                                        verbose=len(p0_2))
except Exception as e:
    print("Primera optimización 2D falló:", e)
    # Reintentar con otros puntos iniciales (prueba 2)
    p0_try = [0.1, 0.01]
    print("Reintentando con p0 =", p0_try)
    try:
        popt2 = dadi.Inference.optimize_log(p0_try, data_2d, func2_ex, pts_l,
                                            lower_bound=lower2, upper_bound=upper2,
                                            verbose=len(p0_try))
    except Exception as e2:
        print("Segunda optimización 2D también falló:", e2)
        print("Abandonando optimización 2D y saliendo con error.")
        raise e2

# Si llegamos aquí, tenemos popt2
model2 = func2_ex(popt2, data_2d.sample_sizes, pts_l)

# Cálculo de likelihood y theta utilizando funciones disponibles
try:
    ll2 = dadi.Inference.ll_multinom(model2, data_2d)
except Exception:
    # fallback: si ll_multinom falla, intentar ll (menos común)
    try:
        ll2 = dadi.Inference.ll(model2, data_2d)
    except Exception:
        ll2 = float("nan")

try:
    theta2 = dadi.Inference.optimal_sfs_scaling(model2, data_2d)
except Exception:
    theta2 = float("nan")

print("Parámetros óptimos (2D):", popt2)
print("Log-likelihood (2D):", ll2)
print("Theta (2D):", theta2)

# ---------------------------
# Guardar plots 2D (compatibles con tu versión)
# ---------------------------
# 1) Comparación modelo vs datos (Poisson o multinom)
plt.figure(figsize=(6, 6))
if hasattr(dadi.Plotting, "plot_2d_comp_Poisson"):
    dadi.Plotting.plot_2d_comp_Poisson(model2, data_2d)
elif hasattr(dadi.Plotting, "plot_2d_comp_multinom"):
    dadi.Plotting.plot_2d_comp_multinom(model2, data_2d)
else:
    print("No hay función de comparación 2D disponible en dadi.Plotting")
plt.title("2D: Modelo vs Datos")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "2D_model_comp.png"), dpi=300)
plt.close()

# 2) SFS observado (2D)
plt.figure(figsize=(6, 6))
if hasattr(dadi.Plotting, "plot_single_2d_sfs"):
    dadi.Plotting.plot_single_2d_sfs(data_2d)
else:
    print("plot_single_2d_sfs no disponible")
plt.title("2D: SFS observado")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "2D_SFS_observed.png"), dpi=300)
plt.close()

print("\nResultados y figuras guardadas en:", OUTDIR)
print("Fin del script.")


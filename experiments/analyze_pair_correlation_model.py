#!/usr/bin/env python3
"""
First pair-correlation experiment built on the q=997 wheel model.

This is a count-level Gaussian pair-covariance model, not yet a full
Hardy-Littlewood point-process model.

Input:
    first_1000_prime_square_interval_index.csv

Outputs:
    pair_correlation_model_comparison.csv
    pair_correlation_rho_summary.csv
    pair_correlation_range_sweep.csv
"""
import math
import numpy as np
import pandas as pd
from scipy.special import expi
from pathlib import Path

ROOT = Path(__file__).resolve().parent
df = pd.read_csv(ROOT/"first_1000_prime_square_interval_index.csv")
s = df[df["p_n"] >= 4001].copy().reset_index(drop=True)

n = s["الترتيب_n"].to_numpy(dtype=int)
p = s["p_n"].to_numpy(dtype=int)
p_next = s["p_التالي"].to_numpy(dtype=int)
A = s["p_n^2"].to_numpy(dtype=np.int64)
B = s["p_التالي^2"].to_numpy(dtype=np.int64)
H = (B-A).astype(float)
obs = s["عدد_الأوليات_بين_p2_ومربع_التالي"].to_numpy(dtype=float)

# Exact q=997 wheel survivors.
max_x = int(B.max())
small = np.ones(998,dtype=bool)
small[:2] = False
for k in range(2,math.isqrt(997)+1):
    if small[k]:
        small[k*k:998:k] = False
wheel = np.flatnonzero(small)

surv = np.ones(max_x+1,dtype=bool)
surv[:2] = False
for q in wheel:
    q=int(q)
    surv[q*q:max_x+1:q] = False
C = np.fromiter(
    (np.count_nonzero(surv[int(a)+1:int(b)]) for a,b in zip(A,B)),
    dtype=np.int64,count=len(A)
)

base_mu = expi(np.log(B.astype(float))) - expi(np.log(A.astype(float)))
alpha = obs.sum()/base_mu.sum()
mu = alpha*base_mu
var = mu*(1-mu/C)
z_obs=(obs-mu)/np.sqrt(var)

def V(h,X):
    if h <= 0:
        return 0.0
    return h*math.log(X/h)

def pair_rho(i,j):
    a,b=float(A[i]),float(B[i])
    c,d=float(A[j]),float(B[j])
    Xloc=math.sqrt(((a+b)/2)*((c+d)/2))
    cov=0.5*(V(d-a,Xloc)+V(c-b,Xloc)-V(c-a,Xloc)-V(d-b,Xloc))
    return cov/math.sqrt(V(b-a,Xloc)*V(d-c,Xloc))

def build_corr(max_lag):
    N=len(A)
    R=np.eye(N)
    for lag in range(1,max_lag+1):
        for i in range(N-lag):
            r=pair_rho(i,i+lag)
            R[i,i+lag]=r
            R[i+lag,i]=r
    return R

# The repository CSVs contain the full Monte Carlo sweep.
for L in [1,2,3,5,10,20,50]:
    R=build_corr(L)
    print(L, np.linalg.eigvalsh(R).min())

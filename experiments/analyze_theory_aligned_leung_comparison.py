#!/usr/bin/env python3
"""
Theory-aligned comparison with Leung (2026) and the count-level
pair-Gaussian short-interval model.

Inputs:
    primes_up_to_1001st_prime_square.csv
    negative_correlation_details.csv

Outputs:
    leung2026_weighted_interval_summary.csv
    leung2026_lag_correlation_comparison.csv
    leung2026_triple_ordering_comparison.csv
    prime_square_pair_gaussian_higher_order_diagnostics.csv
    prime_square_pair_gaussian_sign_triples.csv
"""

import math
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
PRIME_FILE = ROOT / "primes_up_to_1001st_prime_square.csv"
DETAIL_FILE = ROOT / "negative_correlation_details.csv"

primes = pd.read_csv(PRIME_FILE, usecols=["prime"])["prime"].to_numpy(np.int64)
max_y = 62_800_000
powers, weights = [], []
for p0 in primes:
    p = int(p0)
    if p > max_y:
        break
    v = p
    lp = math.log(p)
    while v <= max_y:
        powers.append(v)
        weights.append(lp)
        if v > max_y // p:
            break
        v *= p

order = np.argsort(powers)
pp = np.asarray(powers, dtype=np.int64)[order]
ww = np.asarray(weights, dtype=float)[order]
cw = np.concatenate([[0.0], np.cumsum(ww)])

def psi(y):
    idx = np.searchsorted(pp, np.asarray(y, dtype=float), side="right")
    return cw[idx]

gamma = 0.5772156649015328606
c0 = 1 - gamma - math.log(2*math.pi)

def Delta(t):
    c = 0.0 if t == 1 else (t-1)*math.log(t-1)
    return 0.5*((t+1)*math.log(t+1)-2*t*math.log(t)+c)

deltas = [0.001,0.002,0.004,0.008]
r = 4
samples = 8000
u = np.linspace(math.log(20_000_000.0), math.log(55_000_000.0), samples)
x = np.exp(u)
rng = np.random.default_rng(20260901)

summary_rows, lag_rows, ordering_rows = [], [], []

for delta in deltas:
    Vdelta = delta*math.log(1/delta) + c0*delta
    Z = np.empty((samples,r))
    for j in range(r):
        lo=(1+(j-0.5)*delta)*x
        hi=(1+(j+0.5)*delta)*x
        E=(psi(hi)-psi(lo)-delta*x)/np.sqrt(x)
        Z[:,j]=E/math.sqrt(Vdelta)

    corr=np.corrcoef(Z,rowvar=False)
    Dlead=math.log(1/delta)
    Dsec=Dlead+c0

    Csec=np.eye(r)
    for j in range(r):
        for k in range(j+1,r):
            d=abs(j-k)
            Csec[j,k]=Csec[k,j]=-Delta(d)/Dsec

    G=rng.multivariate_normal(np.zeros(r),Csec,size=400_000)
    emp_runs=1+np.sum((Z[:,1:]>=0)!=(Z[:,:-1]>=0),axis=1)
    th_runs=1+np.sum((G[:,1:]>=0)!=(G[:,:-1]>=0),axis=1)

    summary_rows.append({
        "delta":delta,
        "empirical_variance_z_avg":float(Z.var(axis=0,ddof=1).mean()),
        "empirical_runs4_mean":float(emp_runs.mean()),
        "gaussian_theory_runs4_mean":float(th_runs.mean()),
    })

    for d in range(1,r):
        vals=[corr[j,j+d] for j in range(r-d)]
        lag_rows.append({
            "delta":delta,
            "lag_in_interval_units":d,
            "empirical_corr_mean":float(np.mean(vals)),
            "leung_leading_prediction":float(-Delta(d)/Dlead),
            "leung_secondary_prediction":float(-Delta(d)/Dsec),
        })

    oe=np.argsort(Z[:,:3],axis=1)
    og=np.argsort(G[:,:3],axis=1)
    for perm in [(0,1,2),(0,2,1),(1,0,2),(1,2,0),(2,0,1),(2,1,0)]:
        pa=np.asarray(perm)
        ordering_rows.append({
            "delta":delta,
            "ordering":"".join(map(str,perm)),
            "empirical_frequency":float(np.mean(np.all(oe==pa,axis=1))),
            "gaussian_theory_frequency":float(np.mean(np.all(og==pa,axis=1))),
        })

pd.DataFrame(summary_rows).to_csv(
    ROOT/"leung2026_weighted_interval_summary.csv",index=False,encoding="utf-8-sig")
pd.DataFrame(lag_rows).to_csv(
    ROOT/"leung2026_lag_correlation_comparison.csv",index=False,encoding="utf-8-sig")
pd.DataFrame(ordering_rows).to_csv(
    ROOT/"leung2026_triple_ordering_comparison.csv",index=False,encoding="utf-8-sig")

d = pd.read_csv(DETAIL_FILE)
z_obs=d["z_A"].to_numpy(float)
p=d["p_n"].to_numpy(np.int64)
pn=d["p_next"].to_numpy(np.int64)
A=p*p
B=pn*pn
centers=.5*(A+B)
N=len(A)

def V(h,X):
    return 0.0 if h<=0 else h*math.log(X/h)

def rho_pair(i,j):
    a,b=float(A[i]),float(B[i])
    c,e=float(A[j]),float(B[j])
    Xloc=math.sqrt(float(centers[i])*float(centers[j]))
    cov=.5*(V(e-a,Xloc)+V(c-b,Xloc)-V(c-a,Xloc)-V(e-b,Xloc))
    return cov/math.sqrt(V(b-a,Xloc)*V(e-c,Xloc))

R=np.eye(N)
for lag in range(1,51):
    for i in range(N-lag):
        rr=rho_pair(i,i+lag)
        R[i,i+lag]=R[i+lag,i]=rr

chol=np.linalg.cholesky(R)
trip=np.column_stack([z_obs[:-2],z_obs[1:-1],z_obs[2:]])
obs_mono=int(np.sum(((trip[:,0]<trip[:,1])&(trip[:,1]<trip[:,2])) |
                    ((trip[:,0]>trip[:,1])&(trip[:,1]>trip[:,2]))))
obs_sign=((trip[:,0]>=0).astype(int)*4+
          (trip[:,1]>=0).astype(int)*2+
          (trip[:,2]>=0).astype(int))
obs_counts=np.bincount(obs_sign,minlength=8)

REPS=50_000
BATCH=250
rng=np.random.default_rng(20260901)
mono=np.empty(REPS,dtype=np.int16)
runs=np.empty(REPS,dtype=np.int16)
patterns=np.empty((REPS,8),dtype=np.int16)

for st in range(0,REPS,BATCH):
    m=min(BATCH,REPS-st)
    G=rng.standard_normal((m,N))@chol.T
    g0,g1,g2=G[:,:-2],G[:,1:-1],G[:,2:]
    mono[st:st+m]=np.sum(((g0<g1)&(g1<g2))|((g0>g1)&(g1>g2)),axis=1)
    sg=G>=0
    runs[st:st+m]=1+np.sum(sg[:,1:]!=sg[:,:-1],axis=1)
    codes=sg[:,:-2].astype(np.int8)*4+sg[:,1:-1].astype(np.int8)*2+sg[:,2:].astype(np.int8)
    for k in range(8):
        patterns[st:st+m,k]=np.sum(codes==k,axis=1)

def two_p(obs,sims):
    c=sims.mean()
    return (1+np.sum(np.abs(sims-c)>=abs(obs-c)))/(len(sims)+1)

diag=pd.DataFrame([
    ["runs",256,float(runs.mean()),float(runs.std(ddof=1)),two_p(256,runs)],
    ["monotone_triples",obs_mono,float(mono.mean()),float(mono.std(ddof=1)),two_p(obs_mono,mono)],
],columns=["metric","observed","pair_gaussian_mean","pair_gaussian_sd","two_sided_p"])
diag.to_csv(ROOT/"prime_square_pair_gaussian_higher_order_diagnostics.csv",
            index=False,encoding="utf-8-sig")

rows=[]
for k in range(8):
    rows.append([f"{k:03b}",int(obs_counts[k]),float(patterns[:,k].mean()),
                 float(patterns[:,k].std(ddof=1)),two_p(obs_counts[k],patterns[:,k])])
pd.DataFrame(rows,columns=["sign_pattern","observed_count","pair_gaussian_mean",
                           "pair_gaussian_sd","two_sided_p"]).to_csv(
    ROOT/"prime_square_pair_gaussian_sign_triples.csv",
    index=False,encoding="utf-8-sig")

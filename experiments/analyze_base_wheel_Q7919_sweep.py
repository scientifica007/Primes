#!/usr/bin/env python3
"""
Fixed-Q base-wheel sweep.

For each fixed wheel cutoff q0, remove multiples of all primes <= q0.
For primes q0 < r <= 7919, choose one forbidden residue class modulo r
uniformly at random. Then thin within each prime-square interval so that
the smooth expected prime count mu_i is preserved.

This is a controlled point-level random-sieve experiment designed to
measure the tradeoff between marginal dispersion and pair correlation.

Usage:
    python analyze_base_wheel_Q7919_sweep.py

Edit OUTER/INNER for more or fewer Monte Carlo draws.
"""
import math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import expi

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "first_1000_prime_square_interval_index.csv"

BASES = [2,3,5,7,11,17,29,47,71,97,107,127,139,151,173,181,199,
         223,229,239,251,277,397,499,631,797,997]
Q = 7919
OUTER = 5
INNER = 20
SEED = 20260831

df = pd.read_csv(DATA)
s = df[df["p_n"] >= 4001].copy().reset_index(drop=True)
A = s["p_n^2"].to_numpy(dtype=np.int64)
B = s["p_التالي^2"].to_numpy(dtype=np.int64)
obs = s["عدد_الأوليات_بين_p2_ومربع_التالي"].to_numpy(dtype=np.int64)

start,end=int(A[0]),int(B[-1])
N=end-start
m=len(A)
starts=A-start
ends=B-start

def primes_upto(limit):
    z=np.ones(limit+1,dtype=bool)
    z[:2]=False
    for k in range(2,math.isqrt(limit)+1):
        if z[k]:
            z[k*k:limit+1:k]=False
    return np.flatnonzero(z)

primes=primes_upto(Q)
shape=expi(np.log(B.astype(float)))-expi(np.log(A.astype(float)))
alpha=obs.sum()/shape.sum()
mu=alpha*shape

interval_id=np.empty(N,dtype=np.uint16)
for i,(a,b) in enumerate(zip(starts,ends)):
    interval_id[a:b]=i

boundary=np.unique(np.concatenate([A,B]))-start
boundary=boundary[(boundary>=0)&(boundary<N)]

def lag1_rows(Z):
    x=Z[:,:-1]
    y=Z[:,1:]
    xc=x-x.mean(axis=1,keepdims=True)
    yc=y-y.mean(axis=1,keepdims=True)
    return np.sum(xc*yc,axis=1)/np.sqrt(
        np.sum(xc*xc,axis=1)*np.sum(yc*yc,axis=1)
    )

def simulate(base_q):
    rng=np.random.default_rng(SEED+1009*base_q)

    base=np.ones(N,dtype=bool)
    for q0 in primes[primes<=base_q]:
        q=int(q0)
        base[(-start)%q::q]=False
    base[boundary]=False

    baseC=np.bincount(interval_id[base],minlength=m).astype(np.int64)
    var0=mu*(1-mu/baseC)
    zobs=(obs-mu)/np.sqrt(var0)
    Robs=float(np.sum(zobs*zobs)/(m-1))
    lagobs=float(np.corrcoef(zobs[:-1],zobs[1:])[0,1])
    sg=zobs>=0
    runsobs=int(1+np.sum(sg[1:]!=sg[:-1]))

    tail=primes[(primes>base_q)&(primes<=Q)]
    survival=float(np.prod(1-1/tail.astype(float))) if len(tail) else 1.0
    theta=mu/(baseC*survival)

    Rs,Ls,Us=[],[],[]
    for _ in range(OUTER):
        mask=base.copy()
        counts=baseC.copy()

        for rv in tail:
            r=int(rv)
            forbidden=int(rng.integers(0,r))
            off=(forbidden-start)%r
            sl=slice(off,N,r)
            alive=mask[sl]
            if alive.any():
                counts-=np.bincount(interval_id[sl][alive],minlength=m)
                mask[sl]=False

        Y=rng.binomial(counts,theta,size=(INNER,m)).astype(float)
        Z=(Y-mu)/np.sqrt(var0)
        Rs.extend((np.sum(Z*Z,axis=1)/(m-1)).tolist())
        Ls.extend(lag1_rows(Z).tolist())
        S=Z>=0
        Us.extend((1+np.sum(S[:,1:]!=S[:,:-1],axis=1)).tolist())

    Rs=np.asarray(Rs); Ls=np.asarray(Ls); Us=np.asarray(Us)
    Rm,Rsd=float(Rs.mean()),float(Rs.std(ddof=1))
    Lm,Lsd=float(Ls.mean()),float(Ls.std(ddof=1))
    Um,Usd=float(Us.mean()),float(Us.std(ddof=1))
    zR=(Rm-Robs)/Rsd
    zL=(Lm-lagobs)/Lsd
    zU=(Um-runsobs)/Usd

    return {
        "base_wheel_q0":base_q,
        "Q":Q,
        "simulations":len(Rs),
        "R_observed":Robs,
        "R_mc_mean":Rm,
        "R_mc_sd":Rsd,
        "R_error_z":zR,
        "lag1_observed":lagobs,
        "lag1_mc_mean":Lm,
        "lag1_mc_sd":Lsd,
        "lag1_error_z":zL,
        "runs_observed":runsobs,
        "runs_mc_mean":Um,
        "runs_mc_sd":Usd,
        "runs_error_z":zU,
        "balance_R_lag1":math.hypot(zR,zL),
        "balance_R_lag1_runs":math.sqrt(zR*zR+zL*zL+zU*zU),
    }

out=pd.DataFrame([simulate(q) for q in BASES])
out.to_csv(ROOT/"base_wheel_Q7919_dense_sweep.csv",
           index=False,encoding="utf-8-sig")
print(out.to_string(index=False))

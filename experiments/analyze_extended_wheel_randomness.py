#!/usr/bin/env python3
"""
Extended wheel-conditioning study for prime counts between consecutive prime squares.

Requires:
    first_1000_prime_square_interval_index.csv

Outputs fixed-sample panels that let the wheel cutoff grow much farther than 53.
"""
import math
import numpy as np
import pandas as pd
from scipy.special import expi
from pathlib import Path

ROOT = Path(__file__).resolve().parent
df = pd.read_csv(ROOT/"first_1000_prime_square_interval_index.csv")
base = df[df["p_n"] >= 59].copy().reset_index(drop=True)

p_all = base["p_n"].to_numpy(dtype=int)
a_all = base["p_n^2"].to_numpy(dtype=np.int64)
b_all = base["p_التالي^2"].to_numpy(dtype=np.int64)
obs_all = base["عدد_الأوليات_بين_p2_ومربع_التالي"].to_numpy(dtype=np.int64)
max_x = int(b_all.max())

def primes_upto(N):
    s=np.ones(N+1,dtype=bool)
    s[:2]=False
    for k in range(2,math.isqrt(N)+1):
        if s[k]:
            s[k*k:N+1:k]=False
    return np.flatnonzero(s)

wheel_ps=primes_upto(3989)
targets=[2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,
         71,97,127,173,251,313,397,499,631,797,997,
         1259,1583,1999,2511,3163,3989]
checkpoints=sorted({int(wheel_ps[wheel_ps<=t][-1]) for t in targets})

surv=np.ones(max_x+1,dtype=bool)
surv[:2]=False

def counts(mask):
    return np.fromiter(
        (np.count_nonzero(mask[int(a)+1:int(b)]) for a,b in zip(a_all,b_all)),
        dtype=np.int64,count=len(a_all)
    )

counts_by_q={}
for q in wheel_ps:
    q=int(q)
    surv[q*q:max_x+1:q]=False
    if q in set(checkpoints):
        counts_by_q[q]=counts(surv)

C0=(b_all-a_all-1).astype(np.int64)
panels=[("p>=59",59,53),("p>=257",257,251),
        ("p>=1009",1009,997),("p>=4001",4001,3989)]

rows=[]
for name,pmin,qmax in panels:
    mask=p_all>=pmin
    a,b,obs=a_all[mask],b_all[mask],obs_all[mask]
    base_mu=expi(np.log(b.astype(float)))-expi(np.log(a.astype(float)))
    alpha=obs.sum()/base_mu.sum()
    mu=alpha*base_mu
    stages=[0]+[q for q in checkpoints if q<=qmax]
    r0=None

    for q in stages:
        C=C0[mask] if q==0 else counts_by_q[q][mask]
        prob=mu/C
        valid=bool(np.all(prob<1))
        if valid:
            V=mu*(1-prob)
            z=(obs-mu)/np.sqrt(V)
            R=float(np.sum(z*z)/(len(obs)-1))
            if r0 is None: r0=R
            explained=(R-r0)/(1-r0)
        else:
            R=np.nan
            explained=np.nan

        rows.append({
            "panel":name,
            "sample_min_p":pmin,
            "sample_intervals":len(obs),
            "wheel_through_q":q,
            "model_valid_all_intervals":valid,
            "invalid_interval_count":int(np.sum(prob>=1)),
            "dispersion_ratio_R":R,
            "fraction_initial_dispersion_deficit_explained":explained,
            "surviving_candidates_total":int(C.sum()),
            "actual_primes_total":int(obs.sum()),
            "prime_fraction_among_survivors":float(obs.sum()/C.sum()),
            "max_conditional_probability_mu_over_C":float(prob.max()),
            "median_conditional_probability_mu_over_C":float(np.median(prob)),
        })

pd.DataFrame(rows).to_csv(
    ROOT/"extended_wheel_dispersion_panels.csv",index=False,encoding="utf-8-sig"
)

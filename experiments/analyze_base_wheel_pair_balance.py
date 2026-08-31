#!/usr/bin/env python3
# Reproducible sweep for one fixed base wheel.
# Usage example:
#   python analyze_base_wheel_pair_balance.py 97 40 10
#
# Run separately for q0 = 29, 97, 251, 997 and combine the CSV outputs.

import sys, math, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import expi

base_q = int(sys.argv[1])
OUTER = int(sys.argv[2]) if len(sys.argv) > 2 else 40
INNER = int(sys.argv[3]) if len(sys.argv) > 3 else 10
SEED = 20260831 + base_q

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "first_1000_prime_square_interval_index.csv"
df = pd.read_csv(DATA)
s = df[df["p_n"] >= 4001].copy().reset_index(drop=True)
A = s["p_n^2"].to_numpy(dtype=np.int64)
B = s["p_التالي^2"].to_numpy(dtype=np.int64)
obs = s["عدد_الأوليات_بين_p2_ومربع_التالي"].to_numpy(dtype=np.int64)

start, end = int(A[0]), int(B[-1])
N = end - start
starts = A - start
ends = B - start
m = len(A)

def primes_upto(limit):
    z = np.ones(limit+1, dtype=bool)
    z[:2] = False
    for k in range(2, math.isqrt(limit)+1):
        if z[k]:
            z[k*k:limit+1:k] = False
    return np.flatnonzero(z)

primes = primes_upto(7919)

interval_id = np.empty(N, dtype=np.uint16)
for i,(a,b) in enumerate(zip(starts,ends)):
    interval_id[a:b] = i

shape = expi(np.log(B.astype(float))) - expi(np.log(A.astype(float)))
alpha = obs.sum()/shape.sum()
mu = alpha*shape

stage_grid = [29,47,71,97,127,173,251,397,499,631,797,997,
              1259,1583,1999,2503,3163,3989,5003,6311,7919]
stages = sorted(set(int(primes[primes <= Q][-1]) for Q in stage_grid if Q >= base_q))
stage_set = set(stages)

base = np.ones(N, dtype=bool)
for q0 in primes[primes <= base_q]:
    q = int(q0)
    base[(-start) % q::q] = False
for x in np.unique(np.concatenate([A,B])):
    j = int(x)-start
    if 0 <= j < N:
        base[j] = False

baseC = np.bincount(interval_id[base], minlength=m).astype(np.int64)

base_prob = mu/baseC
base_var = mu*(1-base_prob)
zobs = (obs-mu)/np.sqrt(base_var)
Robs = float(np.sum(zobs*zobs)/(m-1))
lagobs = float(np.corrcoef(zobs[:-1],zobs[1:])[0,1])
sgn=zobs>=0
runsobs=int(1+np.sum(sgn[1:]!=sgn[:-1]))

prod_by_stage={base_q:1.0}
prod=1.0
for r0 in primes[(primes>base_q)&(primes<=max(stages))]:
    r=int(r0); prod*=1-1/r
    if r in stage_set:
        prod_by_stage[r]=prod

stores={Q:{"R":[],"lag":[],"runs":[]} for Q in stages}

def lag1_rows(Z):
    A0=Z[:,:-1]; B0=Z[:,1:]
    Ac=A0-A0.mean(axis=1,keepdims=True)
    Bc=B0-B0.mean(axis=1,keepdims=True)
    return np.sum(Ac*Bc,axis=1)/np.sqrt(np.sum(Ac*Ac,axis=1)*np.sum(Bc*Bc,axis=1))

rng=np.random.default_rng(SEED)
tail=primes[(primes>base_q)&(primes<=max(stages))]
t0=time.time()

for outer in range(OUTER):
    mask=base.copy()
    counts=baseC.copy()

    if base_q in stage_set:
        theta=mu/baseC
        Y=rng.binomial(counts,theta,size=(INNER,m)).astype(float)
        Z=(Y-mu)/np.sqrt(base_var)
        stores[base_q]["R"].extend((np.sum(Z*Z,axis=1)/(m-1)).tolist())
        stores[base_q]["lag"].extend(lag1_rows(Z).tolist())
        S=Z>=0
        stores[base_q]["runs"].extend((1+np.sum(S[:,1:]!=S[:,:-1],axis=1)).tolist())

    for r0 in tail:
        r=int(r0)
        forbidden=int(rng.integers(0,r))
        off=(forbidden-start)%r
        sl=slice(off,N,r)
        alive=mask[sl]
        if alive.any():
            counts -= np.bincount(interval_id[sl][alive], minlength=m)
            mask[sl]=False

        if r in stage_set:
            prod=prod_by_stage[r]
            theta=mu/(baseC*prod)
            if np.any(theta>=1):
                continue
            Y=rng.binomial(counts,theta,size=(INNER,m)).astype(float)
            Z=(Y-mu)/np.sqrt(base_var)
            stores[r]["R"].extend((np.sum(Z*Z,axis=1)/(m-1)).tolist())
            stores[r]["lag"].extend(lag1_rows(Z).tolist())
            S=Z>=0
            stores[r]["runs"].extend((1+np.sum(S[:,1:]!=S[:,:-1],axis=1)).tolist())

rows=[]
for Q in stages:
    ar=np.asarray(stores[Q]["R"],float)
    al=np.asarray(stores[Q]["lag"],float)
    au=np.asarray(stores[Q]["runs"],float)
    if len(ar)==0:
        continue
    Rm,Rsd=float(ar.mean()),float(ar.std(ddof=1))
    Lm,Lsd=float(al.mean()),float(al.std(ddof=1))
    Um,Usd=float(au.mean()),float(au.std(ddof=1))
    zR=(Rm-Robs)/Rsd
    zL=(Lm-lagobs)/Lsd
    zU=(Um-runsobs)/Usd
    prod=1.0 if Q==base_q else prod_by_stage[Q]
    theta=mu/(baseC*prod)
    rows.append({
        "base_wheel_q0":base_q,
        "tail_cutoff_Q":Q,
        "simulations":len(ar),
        "base_survivors_total":int(baseC.sum()),
        "tail_expected_survival_product":float(prod),
        "theta_mean":float(theta.mean()),
        "theta_max":float(theta.max()),
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
        "balance_score_R_lag1":float(math.hypot(zR,zL)),
        "balance_score_R_lag1_runs":float(math.sqrt(zR*zR+zL*zL+zU*zU)),
        "runtime_seconds":time.time()-t0,
    })

out=pd.DataFrame(rows)
path=ROOT / f"base_wheel_{base_q}_pair_balance.csv"
out.to_csv(path,index=False,encoding="utf-8-sig")
best=out.sort_values("balance_score_R_lag1").iloc[0]
print(f"base={base_q} runtime={time.time()-t0:.2f}s reps/stage={OUTER*INNER}")
print(best[["base_wheel_q0","tail_cutoff_Q","R_observed","R_mc_mean","R_error_z",
            "lag1_observed","lag1_mc_mean","lag1_error_z","balance_score_R_lag1"]].to_string())
print(path)

#!/usr/bin/env python3
"""
Wheel-constrained randomness experiment.

Compares prime counts between consecutive prime squares with independent
Bernoulli models defined only on positions surviving successively larger
small-prime wheels.
"""

import numpy as np
import pandas as pd
from scipy.special import expi
from pathlib import Path

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "first_1000_prime_square_interval_index.csv"
REPS = 5000
BATCH = 250
SEED = 20260831

df = pd.read_csv(INPUT)
sample = df[df["p_n"] >= 59].copy()

n = sample["الترتيب_n"].to_numpy(dtype=int)
p = sample["p_n"].to_numpy(dtype=int)
a = sample["p_n^2"].to_numpy(dtype=np.int64)
b = sample["p_التالي^2"].to_numpy(dtype=np.int64)
obs = sample["عدد_الأوليات_بين_p2_ومربع_التالي"].to_numpy(dtype=np.int64)

base_mu = expi(np.log(b.astype(float))) - expi(np.log(a.astype(float)))
alpha = obs.sum() / base_mu.sum()
mu = alpha * base_mu

wheel_primes = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53]
stages = [None] + wheel_primes

def products_with_sign(primes, max_value):
    products = [(1,1)]
    for q in primes:
        current = list(products)
        for d, sign in current:
            nd = d*q
            if nd <= max_value:
                products.append((nd,-sign))
    return products

def survivor_counts(primes):
    hi = b-1
    out = np.zeros(len(a), dtype=np.int64)
    prods = products_with_sign(primes, int(hi.max()))
    for d, sign in prods:
        out += sign*((hi//d)-(a//d))
    return out, len(prods)

def lag1(z):
    return np.corrcoef(z[:-1],z[1:])[0,1]

def runs(z):
    s = z >= 0
    return int(1 + np.sum(s[1:] != s[:-1]))

rng = np.random.default_rng(SEED)
rows = []

for stage in stages:
    primes_used = [] if stage is None else wheel_primes[:wheel_primes.index(stage)+1]
    C, n_products = survivor_counts(primes_used)
    prob = mu/C
    if np.any(prob >= 1):
        raise RuntimeError("invalid probability")

    var = C*prob*(1-prob)
    z = (obs-mu)/np.sqrt(var)
    D_obs = np.sum(z*z)
    lag_obs = lag1(z)
    runs_obs = runs(z)

    D_sim = np.empty(REPS)
    lag_sim = np.empty(REPS)
    runs_sim = np.empty(REPS, dtype=np.int32)

    for start in range(0, REPS, BATCH):
        m = min(BATCH, REPS-start)
        Y = rng.binomial(C, prob, size=(m,len(C))).astype(float)

        # Re-fit the one global density parameter in each bootstrap sample.
        alpha_sim = Y.sum(axis=1)/base_mu.sum()
        MU = alpha_sim[:,None]*base_mu[None,:]
        P = MU/C[None,:]
        V = MU*(1-P)
        Z = (Y-MU)/np.sqrt(V)

        D_sim[start:start+m] = np.sum(Z*Z, axis=1)

        A, Bz = Z[:,:-1], Z[:,1:]
        Ac = A-A.mean(axis=1, keepdims=True)
        Bc = Bz-Bz.mean(axis=1, keepdims=True)
        lag_sim[start:start+m] = np.sum(Ac*Bc,axis=1)/np.sqrt(
            np.sum(Ac*Ac,axis=1)*np.sum(Bc*Bc,axis=1)
        )

        S = Z >= 0
        runs_sim[start:start+m] = 1 + np.sum(S[:,1:] != S[:,:-1],axis=1)

    p_under = (1 + np.sum(D_sim <= D_obs))/(REPS+1)
    p_lag = (1 + np.sum(np.abs(lag_sim) >= abs(lag_obs)))/(REPS+1)
    p_runs = (
        1 + np.sum(np.abs(runs_sim-runs_sim.mean()) >= abs(runs_obs-runs_sim.mean()))
    )/(REPS+1)

    rows.append({
        "wheel_through_prime": 0 if stage is None else stage,
        "wheel_primes": "" if stage is None else " ".join(map(str,primes_used)),
        "inclusion_exclusion_products": n_products,
        "total_surviving_candidates": int(C.sum()),
        "candidate_fraction_of_all_positions": float(C.sum()/np.sum(b-a-1)),
        "min_conditional_prime_probability": float(prob.min()),
        "max_conditional_prime_probability": float(prob.max()),
        "dispersion_D_over_df": float(D_obs/(len(obs)-1)),
        "lag1_autocorrelation": float(lag_obs),
        "runs_count": runs_obs,
        "mc_p_underdispersion": float(p_under),
        "mc_p_lag1_two_sided": float(p_lag),
        "mc_p_runs_two_sided": float(p_runs),
        "mc_mean_D_over_df": float(D_sim.mean()/(len(obs)-1)),
        "mc_sd_D_over_df": float(D_sim.std(ddof=1)/(len(obs)-1)),
    })

out = pd.DataFrame(rows)
r0 = out.loc[0,"dispersion_D_over_df"]
out["fraction_initial_dispersion_deficit_explained"] = (
    (out["dispersion_D_over_df"]-r0)/(1-r0)
)
out.to_csv(ROOT/"wheel_constrained_randomness_progression.csv",
           index=False, encoding="utf-8-sig")
print(out.to_string(index=False))

"""دراسة تجريبية لتمثيل الأعداد الأولية في قواعد عد مختلفة.

البرنامج يفصل بين:
1) بنية suffix/modular الدقيقة على جميع الأوليات حتى N.
2) خصائص رقمية أعمق (انتقالات الأرقام، مجموع الأرقام، المجموع المتناوب)
   على عينة متطابقة الحجم من الأوليات والمركبات لتجنب كلفة O(N * عدد القواعد).

يتطلب NumPy وMatplotlib.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from math import gcd, isqrt
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np

SPECIAL_BASES = (6, 30, 210, 2310)
DEFAULT_SAMPLE_SIZE = 200_000
DENSE_SUFFIX_LIMIT = 2_000_000
DENSE_TRANSITION_BASE_LIMIT = 128


def prime_factors(n: int) -> list[int]:
    out: list[int] = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        out.append(n)
    return out


def euler_phi(n: int) -> int:
    result = n
    for p in prime_factors(n):
        result -= result // p
    return result


def prime_sieve(limit: int) -> tuple[np.ndarray, np.ndarray]:
    if limit < 2:
        sieve = np.zeros(limit + 1, dtype=np.bool_)
        return np.empty(0, dtype=np.int64), sieve
    sieve = np.ones(limit + 1, dtype=np.bool_)
    sieve[:2] = False
    sieve[4::2] = False
    for p in range(3, isqrt(limit) + 1, 2):
        if sieve[p]:
            sieve[p * p :: 2 * p] = False
    primes = np.flatnonzero(sieve).astype(np.int64, copy=False)
    return primes, sieve


def exact_candidate_count(limit: int, base: int) -> int:
    """عدد n في [2,N] التي تحقق gcd(n,base)=1."""
    if limit < 2:
        return 0
    q, r = divmod(limit, base)
    count = q * euler_phi(base)
    count += sum(1 for x in range(1, r + 1) if gcd(x, base) == 1)
    return count - 1  # استبعاد n=1


def regular_primes(primes: np.ndarray, base: int) -> np.ndarray:
    factors = prime_factors(base)
    if not factors or primes.size == 0:
        return primes
    mask = np.ones(primes.size, dtype=np.bool_)
    for p in factors:
        mask &= primes != p
    return primes[mask]


def entropy_from_counts(counts: np.ndarray) -> float:
    counts = counts[counts > 0]
    total = int(counts.sum())
    if total <= 0:
        return 0.0
    probs = counts.astype(np.float64) / total
    return float(-(probs * np.log2(probs)).sum())


def residue_total_counts(residues: np.ndarray, modulus: int, limit: int) -> np.ndarray:
    residues = residues.astype(np.int64, copy=False)
    first = residues.copy()
    need = first < 2
    if np.any(need):
        delta = 2 - first[need]
        first[need] += ((delta + modulus - 1) // modulus) * modulus
    totals = np.zeros(first.size, dtype=np.int64)
    valid = first <= limit
    totals[valid] = ((limit - first[valid]) // modulus) + 1
    return totals


def residue_to_digits(residue: int, base: int, width: int) -> str:
    digits = [0] * width
    x = int(residue)
    for i in range(width - 1, -1, -1):
        x, d = divmod(x, base)
        digits[i] = d
    return ",".join(str(d) for d in digits)


def suffix_distribution(primes: np.ndarray, limit: int, base: int, width: int) -> dict:
    ps = regular_primes(primes, base)
    modulus = base**width
    residues = ps % modulus

    if modulus <= DENSE_SUFFIX_LIMIT:
        dense = np.bincount(residues, minlength=modulus)
        nz = np.flatnonzero(dense)
        counts = dense[nz].astype(np.int64, copy=False)
        values = nz.astype(np.int64, copy=False)
    else:
        values, counts = np.unique(residues, return_counts=True)
        values = values.astype(np.int64, copy=False)
        counts = counts.astype(np.int64, copy=False)

    n = int(ps.size)
    theoretical_support = euler_phi(base) * (base ** (width - 1))
    observed_support = int(values.size)
    expected = n / theoretical_support if theoretical_support else 0.0
    h = entropy_from_counts(counts)
    hmax = math.log2(theoretical_support) if theoretical_support > 1 else 0.0
    kl_uniform_bits = hmax - h if hmax else 0.0

    if expected > 0 and theoretical_support > 1:
        sum_sq = float(np.square(counts.astype(np.float64)).sum())
        chi2 = sum_sq / expected - n
        reduced_chi2 = chi2 / (theoretical_support - 1)
    else:
        reduced_chi2 = 0.0

    order_desc = np.argsort(counts)[::-1]
    order_asc = np.argsort(counts)
    top = [
        {
            "pattern": residue_to_digits(int(values[i]), base, width),
            "residue": int(values[i]),
            "count": int(counts[i]),
        }
        for i in order_desc[:5]
    ]
    rare = [
        {
            "pattern": residue_to_digits(int(values[i]), base, width),
            "residue": int(values[i]),
            "count": int(counts[i]),
        }
        for i in order_asc[:5]
    ]

    candidate_count = exact_candidate_count(limit, base)
    base_prime_rate = n / candidate_count if candidate_count else 0.0
    totals = residue_total_counts(values, modulus, limit)
    usable = totals >= 100
    if base_prime_rate > 0 and np.any(usable):
        rates = counts[usable] / totals[usable]
        lifts = rates / base_prime_rate
        max_lift = float(lifts.max())
        min_lift = float(lifts.min())
        lift_cv = float(lifts.std() / lifts.mean()) if lifts.mean() else 0.0
    else:
        max_lift = min_lift = lift_cv = float("nan")

    return {
        "base": base,
        "width": width,
        "modulus": modulus,
        "prime_count": n,
        "theoretical_support": theoretical_support,
        "observed_support": observed_support,
        "zero_allowed_patterns": theoretical_support - observed_support,
        "expected_per_pattern": expected,
        "entropy_bits": h,
        "normalized_entropy": h / hmax if hmax else 1.0,
        "kl_from_uniform_bits": kl_uniform_bits,
        "reduced_chi2": reduced_chi2,
        "base_prime_rate": base_prime_rate,
        "max_prime_rate_lift": max_lift,
        "min_prime_rate_lift": min_lift,
        "lift_cv": lift_cv,
        "top_patterns": top,
        "rare_patterns": rare,
    }


def sample_composites(sieve: np.ndarray, limit: int, size: int, rng: np.random.Generator) -> np.ndarray:
    chunks: list[np.ndarray] = []
    have = 0
    while have < size:
        draw = rng.integers(2, limit + 1, size=max(10_000, (size - have) * 2), dtype=np.int64)
        comp = draw[~sieve[draw]]
        if comp.size:
            chunks.append(comp)
            have += int(comp.size)
    return np.concatenate(chunks)[:size]


def _add_sparse(target: dict[int, int], keys: np.ndarray) -> None:
    if keys.size == 0:
        return
    unique, counts = np.unique(keys, return_counts=True)
    for k, c in zip(unique.tolist(), counts.tolist()):
        target[int(k)] = target.get(int(k), 0) + int(c)


def _conditional_entropy_dense(matrix: np.ndarray) -> float:
    total = int(matrix.sum())
    if total == 0:
        return 0.0
    rows = matrix.sum(axis=1)
    h = 0.0
    for i in np.flatnonzero(rows):
        row = matrix[i]
        nz = row[row > 0]
        probs = nz.astype(np.float64) / rows[i]
        h_row = float(-(probs * np.log2(probs)).sum())
        h += (rows[i] / total) * h_row
    return h


def _conditional_entropy_sparse(mapping: dict[int, int], base: int) -> float:
    total = sum(mapping.values())
    if total == 0:
        return 0.0
    row_totals: dict[int, int] = defaultdict(int)
    for key, count in mapping.items():
        row_totals[key // base] += count
    h = 0.0
    for key, count in mapping.items():
        row_total = row_totals[key // base]
        p_joint = count / total
        h -= p_joint * math.log2(count / row_total)
    return h


def _js_dense(a: np.ndarray, b: np.ndarray) -> float:
    af = a.astype(np.float64).ravel()
    bf = b.astype(np.float64).ravel()
    sa, sb = af.sum(), bf.sum()
    if sa == 0 or sb == 0:
        return 0.0
    p, q = af / sa, bf / sb
    m = 0.5 * (p + q)
    maskp = p > 0
    maskq = q > 0
    return float(0.5 * np.sum(p[maskp] * np.log2(p[maskp] / m[maskp])) + 0.5 * np.sum(q[maskq] * np.log2(q[maskq] / m[maskq])))


def _js_sparse(a: dict[int, int], b: dict[int, int]) -> float:
    sa, sb = sum(a.values()), sum(b.values())
    if sa == 0 or sb == 0:
        return 0.0
    result = 0.0
    for key in set(a) | set(b):
        p = a.get(key, 0) / sa
        q = b.get(key, 0) / sb
        m = 0.5 * (p + q)
        if p:
            result += 0.5 * p * math.log2(p / m)
        if q:
            result += 0.5 * q * math.log2(q / m)
    return result


def _js_1d(a: np.ndarray, b: np.ndarray) -> float:
    length = max(a.size, b.size)
    aa = np.zeros(length, dtype=np.float64)
    bb = np.zeros(length, dtype=np.float64)
    aa[: a.size] = a
    bb[: b.size] = b
    sa, sb = aa.sum(), bb.sum()
    if sa == 0 or sb == 0:
        return 0.0
    p, q = aa / sa, bb / sb
    m = 0.5 * (p + q)
    maskp, maskq = p > 0, q > 0
    return float(0.5 * np.sum(p[maskp] * np.log2(p[maskp] / m[maskp])) + 0.5 * np.sum(q[maskq] * np.log2(q[maskq] / m[maskq])))


def digit_features(values: np.ndarray, base: int) -> dict:
    x = values.astype(np.int64, copy=True)
    digit_sum = np.zeros(x.size, dtype=np.int32)
    alternating = np.zeros(x.size, dtype=np.int32)
    sign = 1
    position = 0
    dense_mode = base <= DENSE_TRANSITION_BASE_LIMIT
    if dense_mode:
        transitions = np.zeros((base, base), dtype=np.int64)
        internal = np.zeros((base, base), dtype=np.int64)
    else:
        transitions = {}
        internal = {}

    while np.any(x > 0):
        digit = x % base
        digit_sum += digit.astype(np.int32)
        alternating += sign * digit.astype(np.int32)
        q = x // base
        mask = q > 0
        if np.any(mask):
            next_digit = q[mask] % base
            keys = digit[mask] * base + next_digit
            if dense_mode:
                bc = np.bincount(keys, minlength=base * base).reshape(base, base)
                transitions += bc
                if position >= 1:
                    internal += bc
            else:
                _add_sparse(transitions, keys)
                if position >= 1:
                    _add_sparse(internal, keys)
        x = q
        sign = -sign
        position += 1

    sum_counts = np.bincount(digit_sum)
    return {
        "digit_sum": digit_sum,
        "alternating": alternating,
        "sum_counts": sum_counts,
        "transitions": transitions,
        "internal": internal,
        "positions": position,
    }


def compare_digit_structure(prime_sample: np.ndarray, composite_sample: np.ndarray, base: int) -> tuple[dict, list[dict]]:
    pf = digit_features(prime_sample, base)
    cf = digit_features(composite_sample, base)
    if isinstance(pf["transitions"], np.ndarray):
        trans_js = _js_dense(pf["transitions"], cf["transitions"])
        internal_js = _js_dense(pf["internal"], cf["internal"])
        trans_h = _conditional_entropy_dense(pf["transitions"])
        internal_h = _conditional_entropy_dense(pf["internal"])
    else:
        trans_js = _js_sparse(pf["transitions"], cf["transitions"])
        internal_js = _js_sparse(pf["internal"], cf["internal"])
        trans_h = _conditional_entropy_sparse(pf["transitions"], base)
        internal_h = _conditional_entropy_sparse(pf["internal"], base)

    digit_sum_js = _js_1d(pf["sum_counts"], cf["sum_counts"])
    logb = math.log2(base) if base > 1 else 1.0
    row = {
        "base": base,
        "sample_size": int(prime_sample.size),
        "prime_digit_sum_mean": float(pf["digit_sum"].mean()),
        "prime_digit_sum_std": float(pf["digit_sum"].std()),
        "composite_digit_sum_mean": float(cf["digit_sum"].mean()),
        "composite_digit_sum_std": float(cf["digit_sum"].std()),
        "digit_sum_js_bits": digit_sum_js,
        "transition_js_bits": trans_js,
        "internal_transition_js_bits": internal_js,
        "prime_transition_cond_entropy": trans_h,
        "prime_transition_cond_entropy_norm": trans_h / logb,
        "prime_internal_cond_entropy": internal_h,
        "prime_internal_cond_entropy_norm": internal_h / logb,
        "digit_sum_mod_b_minus_1_mismatches": int(np.count_nonzero((pf["digit_sum"] - prime_sample) % max(1, base - 1))),
        "alternating_mod_b_plus_1_mismatches": int(np.count_nonzero((pf["alternating"] - prime_sample) % (base + 1))),
    }

    tests: list[dict] = []
    for relation, modulus, pvals, cvals in (
        ("digit_sum", base - 1, pf["digit_sum"], cf["digit_sum"]),
        ("alternating_sum", base + 1, pf["alternating"], cf["alternating"]),
    ):
        if modulus <= 1:
            continue
        for q in prime_factors(modulus):
            tests.append(
                {
                    "base": base,
                    "relation": relation,
                    "modulus": modulus,
                    "prime_factor": q,
                    "prime_zero_fraction": float(np.mean((pvals % q) == 0)),
                    "composite_zero_fraction": float(np.mean((cvals % q) == 0)),
                }
            )
    return row, tests


def primorial_rows(limit: int) -> list[dict]:
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23]
    product = 1
    rows: list[dict] = []
    for p in primes:
        product *= p
        if product < 6:
            continue
        phi = euler_phi(product)
        density = phi / product
        rows.append(
            {
                "largest_prime": p,
                "primorial": product,
                "phi": phi,
                "density": density,
                "removed_fraction": 1 - density,
                "candidate_sites_at_N": limit * density,
                "bit_bytes_at_N": limit * density / 8,
                "relative_to_odd_bitset": 2 * density,
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    if not rows:
        return
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def make_plots(outdir: Path, base_rows: list[dict], suffix_rows: list[dict], digit_rows: list[dict], primorial: list[dict]) -> None:
    small = [r for r in base_rows if r["base"] <= 100]
    x = [r["base"] for r in small]
    y = [r["phi_over_b"] for r in small]
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(x, y, marker=".", linestyle="none")
    ax.set_xlabel("Base b")
    ax.set_ylabel("phi(b) / b")
    ax.set_title("Candidate density by base")
    fig.tight_layout()
    fig.savefig(outdir / "candidate_density_by_base.png", dpi=160)
    plt.close(fig)

    suffix1 = [r for r in suffix_rows if r["width"] == 1 and r["base"] <= 100]
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot([r["base"] for r in suffix1], [r["kl_from_uniform_bits"] for r in suffix1], marker=".", linestyle="none")
    ax.set_xlabel("Base b")
    ax.set_ylabel("KL from uniform over allowed last digits (bits)")
    ax.set_title("Finite-range last-digit deviation after modular conditioning")
    fig.tight_layout()
    fig.savefig(outdir / "last_digit_kl_by_base.png", dpi=160)
    plt.close(fig)

    dsmall = [r for r in digit_rows if r["base"] <= 100]
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot([r["base"] for r in dsmall], [r["internal_transition_js_bits"] for r in dsmall], marker=".", linestyle="none")
    ax.set_xlabel("Base b")
    ax.set_ylabel("JS divergence: prime vs composite internal digit transitions (bits)")
    ax.set_title("Non-terminal digit-transition discrimination")
    fig.tight_layout()
    fig.savefig(outdir / "internal_transition_js_by_base.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot([r["primorial"] for r in primorial], [r["density"] for r in primorial], marker="o")
    ax.set_xscale("log")
    ax.set_xlabel("Primorial M")
    ax.set_ylabel("phi(M) / M")
    ax.set_title("Candidate density along primorial wheels")
    fig.tight_layout()
    fig.savefig(outdir / "primorial_density.png", dpi=160)
    plt.close(fig)


def run(limit: int, bases_max: int, sample_size: int, output_dir: Path, seed: int) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    primes, sieve = prime_sieve(limit)
    bases = sorted(set(range(2, bases_max + 1)) | set(SPECIAL_BASES))

    rng = np.random.default_rng(seed)
    matched = min(sample_size, int(primes.size))
    if primes.size > matched:
        idx = rng.choice(primes.size, size=matched, replace=False)
        prime_sample = primes[idx]
    else:
        prime_sample = primes.copy()
    composite_sample = sample_composites(sieve, limit, matched, rng) if matched else np.empty(0, dtype=np.int64)

    base_rows: list[dict] = []
    suffix_rows: list[dict] = []
    extreme_rows: list[dict] = []
    digit_rows: list[dict] = []
    divisibility_rows: list[dict] = []

    for base in bases:
        phi = euler_phi(base)
        candidates = exact_candidate_count(limit, base)
        factors = prime_factors(base)
        possible_last = [r for r in range(base) if gcd(r, base) == 1]
        base_rows.append(
            {
                "base": base,
                "prime_factors": "*".join(map(str, factors)),
                "radical": math.prod(factors),
                "phi": phi,
                "phi_over_b": phi / base,
                "candidate_count_exact": candidates,
                "candidate_density_exact": candidates / max(1, limit - 1),
                "possible_last_digit_count": len(possible_last),
                "possible_last_digits": ",".join(map(str, possible_last)),
                "candidate_bit_bytes_at_N": candidates / 8,
                "relative_candidate_sites_vs_odd": candidates / ((limit - 1) / 2),
            }
        )

        for width in (1, 2, 3):
            m = suffix_distribution(primes, limit, base, width)
            suffix_rows.append({k: v for k, v in m.items() if k not in ("top_patterns", "rare_patterns")})
            for kind in ("top_patterns", "rare_patterns"):
                for rank, item in enumerate(m[kind], 1):
                    extreme_rows.append(
                        {
                            "base": base,
                            "width": width,
                            "kind": "top" if kind == "top_patterns" else "rare_nonzero",
                            "rank": rank,
                            **item,
                        }
                    )

        drow, tests = compare_digit_structure(prime_sample, composite_sample, base)
        digit_rows.append(drow)
        divisibility_rows.extend(tests)

    primorial = primorial_rows(limit)

    write_csv(output_dir / "base_summary.csv", base_rows)
    write_csv(output_dir / "suffix_metrics.csv", suffix_rows)
    write_csv(output_dir / "suffix_extremes.csv", extreme_rows)
    write_csv(output_dir / "digit_structure.csv", digit_rows)
    write_csv(output_dir / "digit_divisibility_tests.csv", divisibility_rows)
    write_csv(output_dir / "primorial_chain.csv", primorial)
    make_plots(output_dir, base_rows, suffix_rows, digit_rows, primorial)

    small_bases = [r for r in base_rows if r["base"] <= 100]
    best_prefilter = sorted(small_bases, key=lambda r: (r["phi_over_b"], r["base"]))[:15]
    reliable_suffix1 = [r for r in suffix_rows if r["width"] == 1 and r["base"] <= 100 and r["expected_per_pattern"] >= 100]
    strongest_last_digit_bias = sorted(reliable_suffix1, key=lambda r: r["kl_from_uniform_bits"], reverse=True)[:15]
    strongest_internal_js = sorted([r for r in digit_rows if r["base"] <= 100], key=lambda r: r["internal_transition_js_bits"], reverse=True)[:15]
    strongest_digit_sum_js = sorted([r for r in digit_rows if r["base"] <= 100], key=lambda r: r["digit_sum_js_bits"], reverse=True)[:15]

    summary = {
        "limit": limit,
        "prime_count": int(primes.size),
        "last_prime": int(primes[-1]) if primes.size else None,
        "sample_size_each_class": matched,
        "bases": bases,
        "best_prefilter_bases_le_100": best_prefilter,
        "strongest_last_digit_bias_le_100": strongest_last_digit_bias,
        "strongest_internal_transition_js_le_100": strongest_internal_js,
        "strongest_digit_sum_js_le_100": strongest_digit_sum_js,
        "special_bases": {str(b): next(r for r in base_rows if r["base"] == b) for b in SPECIAL_BASES},
        "primorial_chain": primorial,
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Prime representations across bases")
    parser.add_argument("limit", type=int, nargs="?", default=10_000_000)
    parser.add_argument("--bases-max", type=int, default=100)
    parser.add_argument("--sample-size", type=int, default=DEFAULT_SAMPLE_SIZE)
    parser.add_argument("--output-dir", type=Path, default=Path("base-study-output"))
    parser.add_argument("--seed", type=int, default=20260830)
    args = parser.parse_args()
    if args.limit < 10:
        parser.error("limit must be >= 10")
    if args.bases_max < 2:
        parser.error("bases-max must be >= 2")
    if args.sample_size < 1:
        parser.error("sample-size must be >= 1")

    summary = run(args.limit, args.bases_max, args.sample_size, args.output_dir, args.seed)
    print(f"N={summary['limit']:,}")
    print(f"Prime count={summary['prime_count']:,}")
    print(f"Last prime={summary['last_prime']:,}")
    print(f"Matched digit-analysis sample={summary['sample_size_each_class']:,} primes + composites")
    print("Best prefilter bases <=100:")
    for row in summary["best_prefilter_bases_le_100"][:10]:
        print(f"  b={row['base']:>3}: phi(b)/b={row['phi_over_b']:.9f}, factors={row['prime_factors']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

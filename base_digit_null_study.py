"""خط أساس null لتمييز الأنماط الرقمية الحقيقية عن ضجيج العينة.

يقارن:
- prime vs composite
- prime-half A vs prime-half B
- composite-half A vs composite-half B

ثم يطرح متوسط المقارنتين الداخليتين من JS بين الأوليات والمركبات.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from base_representation_study import (
    SPECIAL_BASES,
    _js_1d,
    _js_dense,
    _js_sparse,
    digit_features,
    prime_sieve,
    sample_composites,
)


def js_trans(a: dict, b: dict, key: str) -> float:
    x, y = a[key], b[key]
    if isinstance(x, np.ndarray):
        return _js_dense(x, y)
    return _js_sparse(x, y)


def analyze(limit: int, bases_max: int, sample_size: int, seed: int) -> list[dict]:
    primes, sieve = prime_sieve(limit)
    rng = np.random.default_rng(seed)
    size = min(sample_size, int(primes.size))
    pidx = rng.choice(primes.size, size=size, replace=False)
    ps = primes[pidx]
    cs = sample_composites(sieve, limit, size, rng)
    half = size // 2
    rows: list[dict] = []

    for base in sorted(set(range(2, bases_max + 1)) | set(SPECIAL_BASES)):
        pf = digit_features(ps, base)
        cf = digit_features(cs, base)
        p1 = digit_features(ps[:half], base)
        p2 = digit_features(ps[half : 2 * half], base)
        c1 = digit_features(cs[:half], base)
        c2 = digit_features(cs[half : 2 * half], base)

        cross_t = js_trans(pf, cf, "transitions")
        null_p_t = js_trans(p1, p2, "transitions")
        null_c_t = js_trans(c1, c2, "transitions")
        cross_i = js_trans(pf, cf, "internal")
        null_p_i = js_trans(p1, p2, "internal")
        null_c_i = js_trans(c1, c2, "internal")
        cross_s = _js_1d(pf["sum_counts"], cf["sum_counts"])
        null_p_s = _js_1d(p1["sum_counts"], p2["sum_counts"])
        null_c_s = _js_1d(c1["sum_counts"], c2["sum_counts"])

        rows.append(
            {
                "base": base,
                "sample_size": size,
                "cross_transition_js": cross_t,
                "null_prime_transition_js": null_p_t,
                "null_composite_transition_js": null_c_t,
                "excess_transition_js": cross_t - 0.5 * (null_p_t + null_c_t),
                "cross_internal_js": cross_i,
                "null_prime_internal_js": null_p_i,
                "null_composite_internal_js": null_c_i,
                "excess_internal_js": cross_i - 0.5 * (null_p_i + null_c_i),
                "cross_digit_sum_js": cross_s,
                "null_prime_digit_sum_js": null_p_s,
                "null_composite_digit_sum_js": null_c_s,
                "excess_digit_sum_js": cross_s - 0.5 * (null_p_s + null_c_s),
            }
        )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("limit", type=int, nargs="?", default=10_000_000)
    parser.add_argument("--bases-max", type=int, default=100)
    parser.add_argument("--sample-size", type=int, default=200_000)
    parser.add_argument("--seed", type=int, default=20260831)
    parser.add_argument("--output", type=Path, default=Path("base-study-output/digit_null_baseline.csv"))
    parser.add_argument("--summary", type=Path, default=Path("base-study-output/digit_null_summary.json"))
    args = parser.parse_args()

    rows = analyze(args.limit, args.bases_max, args.sample_size, args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    small = [r for r in rows if r["base"] <= 100]
    summary = {
        "strongest_excess_internal_js": sorted(small, key=lambda r: r["excess_internal_js"], reverse=True)[:20],
        "strongest_excess_digit_sum_js": sorted(small, key=lambda r: r["excess_digit_sum_js"], reverse=True)[:20],
        "special": {str(b): next(r for r in rows if r["base"] == b) for b in SPECIAL_BASES},
    }
    args.summary.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("Top excess internal-transition JS:")
    for r in summary["strongest_excess_internal_js"][:10]:
        print(f"  b={r['base']:>3} excess={r['excess_internal_js']:.9f} cross={r['cross_internal_js']:.9f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Focused timing benchmark for the larger practical primorial wheels.

This run intentionally omits tracemalloc so timing is not distorted by allocation tracing.
Memory was measured separately in benchmark_primorial_wheels.py.
"""
from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path

from primorial_wheel_segmented_method import consume_primes_primorial_wheel

WHEELS = [
    (2, 3, 5, 7),
    (2, 3, 5, 7, 11),
    (2, 3, 5, 7, 11, 13),
]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("limit", type=int, nargs="?", default=100_000_000)
    parser.add_argument("--segment-candidates", type=int, default=32_768)
    parser.add_argument("--csv", type=Path, default=Path("primorial-focus-results.csv"))
    args = parser.parse_args()

    rows = []
    signature = None
    for ps in WHEELS:
        t0 = time.perf_counter()
        count, checksum, last, stats = consume_primes_primorial_wheel(
            args.limit, ps, args.segment_candidates
        )
        elapsed = time.perf_counter() - t0
        sig = (count, checksum, last)
        if signature is None:
            signature = sig
        elif sig != signature:
            raise AssertionError(f"wheel {stats.wheel_modulus} disagrees")
        row = {
            "N": args.limit,
            "M": stats.wheel_modulus,
            "phi": stats.wheel_phi,
            "phi_over_M": stats.wheel_density,
            "prime_count": count,
            "last_prime": last,
            "checksum": checksum,
            "time_seconds": elapsed,
            "segments": stats.segments_processed,
            "represented_candidates": stats.represented_candidates,
            "strike_attempts": stats.strike_attempts,
            "newly_removed": stats.newly_removed,
            "max_segment_storage_bytes": stats.max_segment_storage_bytes,
            "base_storage_bytes": stats.base_storage_bytes,
        }
        rows.append(row)
        print(
            f"M={stats.wheel_modulus:>6} density={stats.wheel_density:.9f} "
            f"time={elapsed:.6f}s segments={stats.segments_processed:,} "
            f"candidates={stats.represented_candidates:,} strikes={stats.strike_attempts:,}"
        )

    with args.csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

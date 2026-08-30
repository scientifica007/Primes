"""Benchmark للسلسلة 6,30,210,2310,30030 بعجلة عامة واحدة."""
from __future__ import annotations

import argparse
import csv
import time
import tracemalloc
from pathlib import Path

from primorial_wheel_segmented_method import consume_primes_primorial_wheel

WHEELS = [
    (2, 3),
    (2, 3, 5),
    (2, 3, 5, 7),
    (2, 3, 5, 7, 11),
    (2, 3, 5, 7, 11, 13),
]


def measure(limit: int, wheel_primes: tuple[int, ...], segment_candidates: int) -> dict:
    tracemalloc.start()
    t0 = time.perf_counter()
    count, checksum, last, stats = consume_primes_primorial_wheel(
        limit, wheel_primes, segment_candidates
    )
    elapsed = time.perf_counter() - t0
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return {
        "N": limit,
        "wheel_primes": "*".join(map(str, wheel_primes)),
        "M": stats.wheel_modulus,
        "phi": stats.wheel_phi,
        "phi_over_M": stats.wheel_density,
        "prime_count": count,
        "last_prime": last,
        "checksum": checksum,
        "time_seconds": elapsed,
        "tracemalloc_peak_bytes": peak,
        "segments": stats.segments_processed,
        "represented_candidates": stats.represented_candidates,
        "strike_attempts": stats.strike_attempts,
        "newly_removed": stats.newly_removed,
        "max_segment_storage_bytes": stats.max_segment_storage_bytes,
        "base_storage_bytes": stats.base_storage_bytes,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("limit", type=int, nargs="?", default=10_000_000)
    parser.add_argument("--segment-candidates", type=int, default=32_768)
    parser.add_argument("--csv", type=Path, default=Path("primorial-wheel-results.csv"))
    args = parser.parse_args()

    rows = []
    signature = None
    for ps in WHEELS:
        row = measure(args.limit, ps, args.segment_candidates)
        current = (row["prime_count"], row["last_prime"], row["checksum"])
        if signature is None:
            signature = current
        elif current != signature:
            raise AssertionError(f"wheel {row['M']} disagrees with reference signature")
        rows.append(row)
        print(
            f"M={row['M']:>6} phi/M={row['phi_over_M']:.9f} "
            f"time={row['time_seconds']:.6f}s peak={row['tracemalloc_peak_bytes']:,} "
            f"candidates={row['represented_candidates']:,} strikes={row['strike_attempts']:,}"
        )

    with args.csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

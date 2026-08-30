"""مقارنة تجريبية بين النسخة الأصلية والنسخة المحسّنة.

يقيس البرنامج زمن التنفيذ وذروة الذاكرة التي يرصدها tracemalloc، ويتحقق كذلك
من أن النسختين تنتجان قائمة الأعداد الأولية نفسها لكل حد تجريبي.
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
import tracemalloc
from dataclasses import asdict
from typing import Callable, Iterable, Tuple

from optimized_method import primes_optimized
from original_method import primes_original


PrimeFunction = Callable[[int], Tuple[list[int], object]]


def measure(fn: PrimeFunction, limit: int) -> dict:
    tracemalloc.start()
    start = time.perf_counter()
    primes, stats = fn(limit)
    elapsed = time.perf_counter() - start
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return {
        "limit": limit,
        "prime_count": len(primes),
        "elapsed_seconds": elapsed,
        "peak_memory_bytes": peak,
        "primes": primes,
        "stats": asdict(stats),
    }


def run_benchmark(limits: Iterable[int], csv_path: str | None = None) -> None:
    rows = []

    for limit in limits:
        original = measure(primes_original, limit)
        optimized = measure(primes_optimized, limit)

        if original["primes"] != optimized["primes"]:
            raise AssertionError(f"Different prime lists for limit={limit}")

        speedup = (
            original["elapsed_seconds"] / optimized["elapsed_seconds"]
            if optimized["elapsed_seconds"] > 0
            else float("inf")
        )

        row = {
            "limit": limit,
            "prime_count": original["prime_count"],
            "original_seconds": original["elapsed_seconds"],
            "optimized_seconds": optimized["elapsed_seconds"],
            "speedup": speedup,
            "original_peak_bytes": original["peak_memory_bytes"],
            "optimized_peak_bytes": optimized["peak_memory_bytes"],
        }
        rows.append(row)

        print(f"N={limit:,}")
        print(f"  primes: {row['prime_count']:,}")
        print(f"  original : {row['original_seconds']:.6f} s | peak={row['original_peak_bytes']:,} B")
        print(f"  optimized: {row['optimized_seconds']:.6f} s | peak={row['optimized_peak_bytes']:,} B")
        print(f"  speedup  : {row['speedup']:.2f}x")
        print(f"  original stats : {original['stats']}")
        print(f"  optimized stats: {optimized['stats']}")
        print()

    if csv_path:
        fieldnames = list(rows[0].keys()) if rows else []
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        print(f"CSV written to: {csv_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark the two prime-search methods")
    parser.add_argument(
        "limits",
        type=int,
        nargs="*",
        default=[1_000, 10_000, 100_000],
        help="Upper limits to test",
    )
    parser.add_argument("--csv", dest="csv_path", help="Optional CSV output path")
    args = parser.parse_args()

    if any(n < 0 for n in args.limits):
        parser.error("limits must be >= 0")

    run_benchmark(args.limits, args.csv_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())

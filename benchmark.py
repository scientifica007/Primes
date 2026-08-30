"""مقارنة تجريبية بين أربع طرق للبحث عن الأعداد الأولية.

الطرق:
1. original_method.py: الصياغة الأصلية القائمة على إنشاء القوائم وطرح المضاعفات.
2. retain_prime_method.py: إبقاء العدد الأولي باستخدام set.
3. retain_prime_compact_method.py: الفكرة نفسها ببت واحد لكل مرشح فردي.
4. optimized_method.py: غربال محسن يعتمد bytearray ويبدأ من p^2.

يقيس البرنامج زمن التنفيذ وذروة الذاكرة بواسطة tracemalloc، ويتحقق من أن الطرق
كلها تنتج قائمة الأعداد الأولية نفسها. عند repeats > 1 يُستخدم وسيط القياسات.
"""

from __future__ import annotations

import argparse
import csv
import gc
import statistics
import sys
import time
import tracemalloc
from dataclasses import asdict
from typing import Callable, Iterable, Tuple

from optimized_method import primes_optimized
from original_method import primes_original
from retain_prime_compact_method import primes_retain_compact
from retain_prime_method import primes_retain


PrimeFunction = Callable[[int], Tuple[list[int], object]]


def measure(fn: PrimeFunction, limit: int, repeats: int = 5) -> dict:
    if repeats < 1:
        raise ValueError("repeats must be >= 1")

    elapsed_samples = []
    peak_samples = []
    last_primes = []
    last_stats = None

    for _ in range(repeats):
        gc.collect()
        tracemalloc.start()
        start = time.perf_counter()
        primes, stats = fn(limit)
        elapsed = time.perf_counter() - start
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        elapsed_samples.append(elapsed)
        peak_samples.append(peak)
        last_primes = primes
        last_stats = stats

    return {
        "limit": limit,
        "prime_count": len(last_primes),
        "elapsed_seconds": statistics.median(elapsed_samples),
        "peak_memory_bytes": int(statistics.median(peak_samples)),
        "primes": last_primes,
        "stats": asdict(last_stats),
    }


def run_benchmark(
    limits: Iterable[int], csv_path: str | None = None, repeats: int = 5
) -> None:
    rows = []

    for limit in limits:
        original = measure(primes_original, limit, repeats)
        retained = measure(primes_retain, limit, repeats)
        compact = measure(primes_retain_compact, limit, repeats)
        optimized = measure(primes_optimized, limit, repeats)

        if not (
            original["primes"]
            == retained["primes"]
            == compact["primes"]
            == optimized["primes"]
        ):
            raise AssertionError(f"Different prime lists for limit={limit}")

        def speedup(base: dict, candidate: dict) -> float:
            if candidate["elapsed_seconds"] <= 0:
                return float("inf")
            return base["elapsed_seconds"] / candidate["elapsed_seconds"]

        row = {
            "limit": limit,
            "prime_count": original["prime_count"],
            "original_seconds": original["elapsed_seconds"],
            "retained_seconds": retained["elapsed_seconds"],
            "compact_seconds": compact["elapsed_seconds"],
            "optimized_seconds": optimized["elapsed_seconds"],
            "original_to_retained_speedup": speedup(original, retained),
            "original_to_compact_speedup": speedup(original, compact),
            "original_to_optimized_speedup": speedup(original, optimized),
            "original_peak_bytes": original["peak_memory_bytes"],
            "retained_peak_bytes": retained["peak_memory_bytes"],
            "compact_peak_bytes": compact["peak_memory_bytes"],
            "optimized_peak_bytes": optimized["peak_memory_bytes"],
        }
        rows.append(row)

        print(f"N={limit:,} | repeats={repeats} | median")
        print(f"  primes   : {row['prime_count']:,}")
        print(f"  original : {row['original_seconds']:.6f} s | peak={row['original_peak_bytes']:,} B")
        print(f"  retained : {row['retained_seconds']:.6f} s | peak={row['retained_peak_bytes']:,} B")
        print(f"  compact  : {row['compact_seconds']:.6f} s | peak={row['compact_peak_bytes']:,} B")
        print(f"  optimized: {row['optimized_seconds']:.6f} s | peak={row['optimized_peak_bytes']:,} B")
        print(f"  original / retained speedup : {row['original_to_retained_speedup']:.2f}x")
        print(f"  original / compact speedup  : {row['original_to_compact_speedup']:.2f}x")
        print(f"  original / optimized speedup: {row['original_to_optimized_speedup']:.2f}x")
        print(f"  original stats : {original['stats']}")
        print(f"  retained stats : {retained['stats']}")
        print(f"  compact stats  : {compact['stats']}")
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
    parser = argparse.ArgumentParser(description="Benchmark the four prime-search methods")
    parser.add_argument(
        "limits",
        type=int,
        nargs="*",
        default=[1_000, 10_000, 100_000],
        help="Upper limits to test",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=5,
        help="Number of repetitions per method and limit; median is reported",
    )
    parser.add_argument("--csv", dest="csv_path", help="Optional CSV output path")
    args = parser.parse_args()

    if any(n < 0 for n in args.limits):
        parser.error("limits must be >= 0")
    if args.repeats < 1:
        parser.error("repeats must be >= 1")

    run_benchmark(args.limits, args.csv_path, args.repeats)
    return 0


if __name__ == "__main__":
    sys.exit(main())

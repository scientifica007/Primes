"""مقارنة تجريبية بين تمثيلات مختلفة لمنهج البحث عن الأعداد الأولية.

تقيس المقارنة الزمن وذروة الذاكرة بواسطة tracemalloc، وتشمل:
1. الصياغة الأصلية.
2. إبقاء العدد الأولي باستخدام set.
3. فضاء مرشحين مضغوط مع list Python كناتج.
4. فضاء مرشحين مضغوط مع array مضغوطة كناتج.
5. واجهة streaming لا تحتفظ بقائمة الناتج.
6. الغربال المحسن المعتمد على bytearray.
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
from typing import Callable, Iterable, Sequence, Tuple

from optimized_method import primes_optimized
from original_method import primes_original
from retain_prime_compact_method import primes_retain_compact
from retain_prime_method import primes_retain
from retain_prime_packed_output_method import (
    iter_primes_retain_packed,
    primes_retain_packed,
)


PrimeFunction = Callable[[int], Tuple[Sequence[int], object]]


def measure(fn: PrimeFunction, limit: int, repeats: int = 5) -> dict:
    if repeats < 1:
        raise ValueError("repeats must be >= 1")

    elapsed_samples = []
    peak_samples = []
    last_primes: Sequence[int] = []
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


def same_values(reference: Sequence[int], candidate: Sequence[int]) -> bool:
    """مقارنة تسلسلين دون تحويل array إلى list إضافية."""
    if len(reference) != len(candidate):
        return False
    return all(a == b for a, b in zip(reference, candidate))


def measure_stream(limit: int, expected: Sequence[int], repeats: int = 5) -> dict:
    """قياس التوليد المتدفق دون الاحتفاظ بقائمة ناتج داخل فترة القياس."""
    elapsed_samples = []
    peak_samples = []
    last_count = 0

    for _ in range(repeats):
        gc.collect()
        tracemalloc.start()
        start = time.perf_counter()

        count = 0
        for count, prime in enumerate(iter_primes_retain_packed(limit), start=1):
            index = count - 1
            if index >= len(expected) or prime != expected[index]:
                tracemalloc.stop()
                raise AssertionError(f"Streaming output differs for limit={limit}")

        if count != len(expected):
            tracemalloc.stop()
            raise AssertionError(f"Streaming output length differs for limit={limit}")

        elapsed = time.perf_counter() - start
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        elapsed_samples.append(elapsed)
        peak_samples.append(peak)
        last_count = count

    return {
        "prime_count": last_count,
        "elapsed_seconds": statistics.median(elapsed_samples),
        "peak_memory_bytes": int(statistics.median(peak_samples)),
    }


def run_benchmark(
    limits: Iterable[int], csv_path: str | None = None, repeats: int = 5
) -> None:
    rows = []

    for limit in limits:
        original = measure(primes_original, limit, repeats)
        retained = measure(primes_retain, limit, repeats)
        compact = measure(primes_retain_compact, limit, repeats)
        packed = measure(primes_retain_packed, limit, repeats)
        optimized = measure(primes_optimized, limit, repeats)

        reference = original["primes"]
        for name, result in (
            ("retained", retained),
            ("compact", compact),
            ("packed", packed),
            ("optimized", optimized),
        ):
            if not same_values(reference, result["primes"]):
                raise AssertionError(f"{name} differs for limit={limit}")

        stream = measure_stream(limit, reference, repeats)

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
            "packed_seconds": packed["elapsed_seconds"],
            "stream_seconds": stream["elapsed_seconds"],
            "optimized_seconds": optimized["elapsed_seconds"],
            "original_to_retained_speedup": speedup(original, retained),
            "original_to_compact_speedup": speedup(original, compact),
            "original_to_packed_speedup": speedup(original, packed),
            "original_to_optimized_speedup": speedup(original, optimized),
            "original_peak_bytes": original["peak_memory_bytes"],
            "retained_peak_bytes": retained["peak_memory_bytes"],
            "compact_peak_bytes": compact["peak_memory_bytes"],
            "packed_peak_bytes": packed["peak_memory_bytes"],
            "stream_peak_bytes": stream["peak_memory_bytes"],
            "optimized_peak_bytes": optimized["peak_memory_bytes"],
        }
        rows.append(row)

        print(f"N={limit:,} | repeats={repeats} | median")
        print(f"  primes   : {row['prime_count']:,}")
        print(f"  original : {row['original_seconds']:.6f} s | peak={row['original_peak_bytes']:,} B")
        print(f"  retained : {row['retained_seconds']:.6f} s | peak={row['retained_peak_bytes']:,} B")
        print(f"  compact  : {row['compact_seconds']:.6f} s | peak={row['compact_peak_bytes']:,} B")
        print(f"  packed   : {row['packed_seconds']:.6f} s | peak={row['packed_peak_bytes']:,} B")
        print(f"  stream   : {row['stream_seconds']:.6f} s | peak={row['stream_peak_bytes']:,} B")
        print(f"  optimized: {row['optimized_seconds']:.6f} s | peak={row['optimized_peak_bytes']:,} B")
        print(f"  original / retained speedup : {row['original_to_retained_speedup']:.2f}x")
        print(f"  original / compact speedup  : {row['original_to_compact_speedup']:.2f}x")
        print(f"  original / packed speedup   : {row['original_to_packed_speedup']:.2f}x")
        print(f"  original / optimized speedup: {row['original_to_optimized_speedup']:.2f}x")
        print(f"  original stats : {original['stats']}")
        print(f"  retained stats : {retained['stats']}")
        print(f"  compact stats  : {compact['stats']}")
        print(f"  packed stats   : {packed['stats']}")
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
    parser = argparse.ArgumentParser(description="Benchmark prime-search representations")
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

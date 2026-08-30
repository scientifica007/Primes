"""Benchmark خاص بتأثير التجزئة على ذاكرة البحث.

يقارن بين:
- streaming كامل المجال: يحتفظ ببت لكل عدد فردي حتى N.
- segmented streaming: يحتفظ بمقطع ثابت السعة فقط، إضافة إلى أوليات sqrt(N).

لا يحتفظ أي من المسارين بقائمة النتائج أثناء القياس. نقارن count/checksum/last
للتأكد من أن التسلسلين متكافئان دون إدخال list كبيرة في ذاكرة القياس.
"""

from __future__ import annotations

import argparse
import csv
import gc
import statistics
import time
import tracemalloc
from typing import Callable, Iterable, Iterator

from retain_prime_packed_output_method import iter_primes_retain_packed
from segmented_method import DEFAULT_SEGMENT_ODDS, iter_primes_segmented


IteratorFactory = Callable[[], Iterator[int]]


def consume(iterator: Iterator[int]) -> tuple[int, int, int]:
    count = 0
    checksum = 0
    last = 0
    for prime in iterator:
        count += 1
        checksum += prime
        last = prime
    return count, checksum, last


def measure(factory: IteratorFactory, repeats: int) -> dict:
    elapsed_samples = []
    peak_samples = []
    signature = (0, 0, 0)

    for _ in range(repeats):
        gc.collect()
        tracemalloc.start()
        start = time.perf_counter()
        current = consume(factory())
        elapsed = time.perf_counter() - start
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        if signature != (0, 0, 0) and current != signature:
            raise AssertionError("Different signatures across benchmark repeats")
        signature = current
        elapsed_samples.append(elapsed)
        peak_samples.append(peak)

    return {
        "count": signature[0],
        "checksum": signature[1],
        "last": signature[2],
        "seconds": statistics.median(elapsed_samples),
        "peak_bytes": int(statistics.median(peak_samples)),
    }


def run_benchmark(
    limits: Iterable[int],
    repeats: int = 5,
    segment_odd_count: int = DEFAULT_SEGMENT_ODDS,
    csv_path: str | None = None,
) -> None:
    rows = []

    for limit in limits:
        full = measure(lambda: iter_primes_retain_packed(limit), repeats)
        segmented = measure(
            lambda: iter_primes_segmented(limit, segment_odd_count), repeats
        )

        if (full["count"], full["checksum"], full["last"]) != (
            segmented["count"],
            segmented["checksum"],
            segmented["last"],
        ):
            raise AssertionError(f"Segmented result differs for N={limit}")

        memory_ratio = (
            full["peak_bytes"] / segmented["peak_bytes"]
            if segmented["peak_bytes"]
            else float("inf")
        )
        time_ratio = (
            full["seconds"] / segmented["seconds"]
            if segmented["seconds"]
            else float("inf")
        )

        row = {
            "limit": limit,
            "prime_count": full["count"],
            "last_prime": full["last"],
            "segment_odd_count": segment_odd_count,
            "full_stream_seconds": full["seconds"],
            "segmented_stream_seconds": segmented["seconds"],
            "full_over_segmented_time_ratio": time_ratio,
            "full_stream_peak_bytes": full["peak_bytes"],
            "segmented_stream_peak_bytes": segmented["peak_bytes"],
            "full_over_segmented_memory_ratio": memory_ratio,
        }
        rows.append(row)

        print(f"N={limit:,} | repeats={repeats} | segment_odds={segment_odd_count:,}")
        print(f"  primes     : {full['count']:,}")
        print(f"  last prime : {full['last']:,}")
        print(
            f"  full stream: {full['seconds']:.6f} s | "
            f"peak={full['peak_bytes']:,} B"
        )
        print(
            f"  segmented  : {segmented['seconds']:.6f} s | "
            f"peak={segmented['peak_bytes']:,} B"
        )
        print(f"  full / segmented memory: {memory_ratio:.2f}x")
        print(f"  full / segmented time  : {time_ratio:.2f}x")
        print()

    if csv_path:
        fieldnames = list(rows[0].keys()) if rows else []
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        print(f"CSV written to: {csv_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark segmented prime streaming")
    parser.add_argument(
        "limits",
        nargs="*",
        type=int,
        default=[100_000, 1_000_000],
        help="Upper limits to test",
    )
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument(
        "--segment-odds",
        type=int,
        default=DEFAULT_SEGMENT_ODDS,
        help="Odd candidates represented in one segment",
    )
    parser.add_argument("--csv", dest="csv_path")
    args = parser.parse_args()

    if any(limit < 0 for limit in args.limits):
        parser.error("limits must be >= 0")
    if args.repeats < 1:
        parser.error("repeats must be >= 1")
    if args.segment_odds < 1:
        parser.error("segment-odds must be >= 1")

    run_benchmark(args.limits, args.repeats, args.segment_odds, args.csv_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

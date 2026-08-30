"""مقارنة Wheel-30 مع Wheel-210 تحت نفس ميزانية البتات لكل مقطع."""

from __future__ import annotations

import argparse
import csv
import gc
import statistics
import time
import tracemalloc
from dataclasses import asdict
from typing import Callable, Iterable

from wheel30_segmented_method import consume_primes_wheel30_segmented
from wheel210_segmented_method import consume_primes_wheel210_segmented


Consumer = Callable[[int, int], tuple[int, int, int, object]]


def measure(consumer: Consumer, limit: int, segment_units: int, repeats: int) -> dict:
    elapsed = []
    peaks = []
    last_result = None
    for _ in range(repeats):
        gc.collect()
        tracemalloc.start()
        start = time.perf_counter()
        result = consumer(limit, segment_units)
        seconds = time.perf_counter() - start
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        elapsed.append(seconds)
        peaks.append(peak)
        last_result = result

    count, checksum, last, stats = last_result
    return {
        "count": count,
        "checksum": checksum,
        "last": last,
        "seconds": statistics.median(elapsed),
        "peak_bytes": int(statistics.median(peaks)),
        "stats": asdict(stats),
    }


def run(limits: Iterable[int], segment_units: int, repeats: int, csv_path: str | None) -> None:
    rows = []
    for limit in limits:
        wheel30 = measure(
            consume_primes_wheel30_segmented, limit, segment_units, repeats
        )
        wheel210 = measure(
            consume_primes_wheel210_segmented, limit, segment_units, repeats
        )

        signature30 = (wheel30["count"], wheel30["checksum"], wheel30["last"])
        signature210 = (wheel210["count"], wheel210["checksum"], wheel210["last"])
        if signature30 != signature210:
            raise AssertionError(f"Wheel-210 differs from Wheel-30 at N={limit}")

        c30 = wheel30["stats"]["represented_candidates"]
        c210 = wheel210["stats"]["represented_candidates"]
        reduction = 1.0 - c210 / c30 if c30 else 0.0

        row = {
            "limit": limit,
            "prime_count": wheel30["count"],
            "last_prime": wheel30["last"],
            "segment_units": segment_units,
            "wheel30_seconds": wheel30["seconds"],
            "wheel210_seconds": wheel210["seconds"],
            "wheel30_peak_bytes": wheel30["peak_bytes"],
            "wheel210_peak_bytes": wheel210["peak_bytes"],
            "wheel30_to_210_time_ratio": wheel30["seconds"] / wheel210["seconds"],
            "wheel30_to_210_memory_ratio": wheel30["peak_bytes"] / wheel210["peak_bytes"],
            "wheel30_candidates": c30,
            "wheel210_candidates": c210,
            "candidate_reduction_percent": reduction * 100.0,
            "wheel30_strike_attempts": wheel30["stats"]["strike_attempts"],
            "wheel210_strike_attempts": wheel210["stats"]["strike_attempts"],
            "wheel30_segments": wheel30["stats"]["segments_processed"],
            "wheel210_segments": wheel210["stats"]["segments_processed"],
        }
        rows.append(row)

        print(f"N={limit:,} | segment bits={segment_units:,} | repeats={repeats}")
        print(f"  primes: {row['prime_count']:,} | last={row['last_prime']:,}")
        print(
            f"  Wheel-30 : {row['wheel30_seconds']:.6f} s | peak={row['wheel30_peak_bytes']:,} B | "
            f"strikes={row['wheel30_strike_attempts']:,} | segments={row['wheel30_segments']:,}"
        )
        print(
            f"  Wheel-210: {row['wheel210_seconds']:.6f} s | peak={row['wheel210_peak_bytes']:,} B | "
            f"strikes={row['wheel210_strike_attempts']:,} | segments={row['wheel210_segments']:,}"
        )
        print(
            f"  represented candidates: {c30:,} -> {c210:,} "
            f"({row['candidate_reduction_percent']:.2f}% fewer)"
        )
        print(f"  Wheel-30 / Wheel-210 time ratio  : {row['wheel30_to_210_time_ratio']:.2f}x")
        print(f"  Wheel-30 / Wheel-210 memory ratio: {row['wheel30_to_210_memory_ratio']:.2f}x")
        print()

    if csv_path:
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else [])
            writer.writeheader()
            writer.writerows(rows)
        print(f"CSV written to: {csv_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark Wheel-30 vs Wheel-210")
    parser.add_argument("limits", type=int, nargs="*", default=[100_000, 1_000_000])
    parser.add_argument("--segment-units", type=int, default=32_768)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--csv", dest="csv_path")
    args = parser.parse_args()

    if any(n < 0 for n in args.limits):
        parser.error("limits must be >= 0")
    if args.segment_units < 1:
        parser.error("segment-units must be >= 1")
    if args.repeats < 1:
        parser.error("repeats must be >= 1")

    run(args.limits, args.segment_units, args.repeats, args.csv_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

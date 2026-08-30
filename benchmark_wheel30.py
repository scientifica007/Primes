"""مقارنة التجزئة الفردية مع التجزئة المبنية على Wheel-30.

كلا الطرفين يستخدم streaming ولا يحتفظ بقائمة ناتج كاملة. كما نستخدم العدد
نفسه من البتات القصوى في المقطع: segment_units بتًا، لكن Wheel-30 تجعل هذه
البتات تمثل فقط الأعداد coprime مع 30، ولذلك يغطي المقطع مجالًا عدديًا أكبر.
"""

from __future__ import annotations

import argparse
import csv
import gc
import statistics
import time
import tracemalloc
from dataclasses import asdict
from typing import Callable, Iterable

from segmented_method import consume_primes_segmented
from wheel30_segmented_method import consume_primes_wheel30_segmented


Consumer = Callable[[int, int], tuple[int, int, int, object]]


def measure(consumer: Consumer, limit: int, segment_units: int, repeats: int) -> dict:
    elapsed_samples = []
    peak_samples = []
    last_result = None

    for _ in range(repeats):
        gc.collect()
        tracemalloc.start()
        start = time.perf_counter()
        count, checksum, last, stats = consumer(limit, segment_units)
        elapsed = time.perf_counter() - start
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        elapsed_samples.append(elapsed)
        peak_samples.append(peak)
        last_result = (count, checksum, last, stats)

    count, checksum, last, stats = last_result
    return {
        "count": count,
        "checksum": checksum,
        "last": last,
        "seconds": statistics.median(elapsed_samples),
        "peak_bytes": int(statistics.median(peak_samples)),
        "stats": asdict(stats),
    }


def run(limits: Iterable[int], segment_units: int, repeats: int, csv_path: str | None) -> None:
    rows = []

    for limit in limits:
        odd = measure(consume_primes_segmented, limit, segment_units, repeats)
        wheel = measure(
            consume_primes_wheel30_segmented, limit, segment_units, repeats
        )

        if (odd["count"], odd["checksum"], odd["last"]) != (
            wheel["count"], wheel["checksum"], wheel["last"]
        ):
            raise AssertionError(f"Wheel-30 differs from odd segmented at N={limit}")

        odd_candidates = max(0, (limit - 1) // 2)
        wheel_candidates = wheel["stats"]["represented_candidates"]
        candidate_reduction = (
            1.0 - wheel_candidates / odd_candidates if odd_candidates else 0.0
        )

        row = {
            "limit": limit,
            "prime_count": odd["count"],
            "last_prime": odd["last"],
            "segment_units": segment_units,
            "odd_seconds": odd["seconds"],
            "wheel30_seconds": wheel["seconds"],
            "odd_peak_bytes": odd["peak_bytes"],
            "wheel30_peak_bytes": wheel["peak_bytes"],
            "odd_to_wheel_speed_ratio": odd["seconds"] / wheel["seconds"],
            "odd_to_wheel_memory_ratio": odd["peak_bytes"] / wheel["peak_bytes"],
            "odd_candidates": odd_candidates,
            "wheel30_candidates": wheel_candidates,
            "candidate_reduction_percent": candidate_reduction * 100.0,
            "odd_strike_attempts": odd["stats"]["strike_attempts"],
            "wheel30_strike_attempts": wheel["stats"]["strike_attempts"],
            "odd_segments": odd["stats"]["segments_processed"],
            "wheel30_segments": wheel["stats"]["segments_processed"],
        }
        rows.append(row)

        print(f"N={limit:,} | segment bits={segment_units:,} | repeats={repeats}")
        print(f"  primes: {row['prime_count']:,} | last={row['last_prime']:,}")
        print(
            f"  odd segmented    : {row['odd_seconds']:.6f} s | "
            f"peak={row['odd_peak_bytes']:,} B | "
            f"strikes={row['odd_strike_attempts']:,} | segments={row['odd_segments']:,}"
        )
        print(
            f"  Wheel-30 segmented: {row['wheel30_seconds']:.6f} s | "
            f"peak={row['wheel30_peak_bytes']:,} B | "
            f"strikes={row['wheel30_strike_attempts']:,} | segments={row['wheel30_segments']:,}"
        )
        print(
            f"  represented candidates: {row['odd_candidates']:,} -> "
            f"{row['wheel30_candidates']:,} "
            f"({row['candidate_reduction_percent']:.2f}% fewer)"
        )
        print(f"  odd / Wheel-30 time ratio  : {row['odd_to_wheel_speed_ratio']:.2f}x")
        print(f"  odd / Wheel-30 memory ratio: {row['odd_to_wheel_memory_ratio']:.2f}x")
        print()

    if csv_path:
        fieldnames = list(rows[0].keys()) if rows else []
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        print(f"CSV written to: {csv_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark odd segmented vs Wheel-30")
    parser.add_argument(
        "limits", type=int, nargs="*", default=[100_000, 1_000_000]
    )
    parser.add_argument(
        "--segment-units",
        type=int,
        default=32_768,
        help="Maximum candidate bits per segment for both methods",
    )
    parser.add_argument("--repeats", type=int, default=3)
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

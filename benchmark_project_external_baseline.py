"""Benchmark selected project configurations for comparison with external sieves.

This script deliberately uses the same project workload as our previous cache studies:
consume every generated prime and compute count + checksum + last prime.
"""
from __future__ import annotations

import argparse
import gc
import statistics
import time

from compact_primorial_wheel_method import consume_primes_compact_spec, make_compact_wheel


CONFIGS = (
    # Best memory/cache-balanced configuration found in the completed 100M sweep.
    (210, (2, 3, 5, 7), 64 * 1024),
    # Best speed configuration within the completed sweep through 1 MiB.
    (30030, (2, 3, 5, 7, 11, 13), 1024 * 1024),
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("limit", type=int, nargs="?", default=100_000_000)
    parser.add_argument("--repeats", type=int, default=2)
    args = parser.parse_args()

    if args.limit < 2:
        raise ValueError("limit must be >= 2")
    if args.repeats < 1:
        raise ValueError("repeats must be >= 1")

    reference = None
    for modulus, primes, segment_bytes in CONFIGS:
        spec = make_compact_wheel(primes)
        timings = []
        result = None
        stats = None
        for _ in range(args.repeats):
            gc.collect()
            start = time.perf_counter()
            count, checksum, last, current_stats = consume_primes_compact_spec(
                args.limit,
                spec,
                segment_bytes * 8,
            )
            elapsed = time.perf_counter() - start
            timings.append(elapsed)
            current = (count, checksum, last)
            if result is None:
                result = current
                stats = current_stats
            elif current != result:
                raise AssertionError("non-deterministic project result")

        if reference is None:
            reference = result
        elif result != reference:
            raise AssertionError("project configurations disagree")

        print(
            f"project wheel={modulus} segment_bytes={segment_bytes} "
            f"count={result[0]} checksum={result[1]} last={result[2]} "
            f"median_seconds={statistics.median(timings):.9f} "
            f"min_seconds={min(timings):.9f} max_seconds={max(timings):.9f} "
            f"working_set_bytes={stats.algorithmic_working_set_bytes} "
            f"strike_attempts={stats.strike_attempts}"
        )


if __name__ == "__main__":
    main()

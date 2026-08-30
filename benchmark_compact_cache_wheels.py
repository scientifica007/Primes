"""Benchmark compact primorial wheels across sieve-segment working-set sizes.

The benchmark deliberately reports two different memory quantities:

* ``segment_bytes``: the hot sieve bitset, chosen from cache-sized values;
* ``algorithmic_working_set_bytes``: compact wheel tables + base-prime storage
  + largest segment payload.  Python interpreter/object overhead is excluded.

A small legacy comparison is also available to quantify how much memory was
accidentally spent on Python tuples/ints before the compact representation.
"""
from __future__ import annotations

import argparse
import csv
import gc
import statistics
import sys
import time
from pathlib import Path

from compact_primorial_wheel_method import (
    consume_primes_compact_spec,
    make_compact_wheel,
)
from primorial_wheel_segmented_method import (
    consume_primes_primorial_wheel,
    make_wheel,
)


WHEELS = {
    30: (2, 3, 5),
    210: (2, 3, 5, 7),
    2310: (2, 3, 5, 7, 11),
    30030: (2, 3, 5, 7, 11, 13),
}
DEFAULT_SEGMENT_BYTES = (
    4 * 1024,
    8 * 1024,
    16 * 1024,
    32 * 1024,
    64 * 1024,
    128 * 1024,
    256 * 1024,
    512 * 1024,
    1024 * 1024,
)


def _deep_size(obj, seen=None) -> int:
    """Approximate unique Python object graph size for legacy wheel tables."""
    if seen is None:
        seen = set()
    oid = id(obj)
    if oid in seen:
        return 0
    seen.add(oid)
    size = sys.getsizeof(obj)
    if isinstance(obj, (tuple, list)):
        size += sum(_deep_size(x, seen) for x in obj)
    return size


def _time_compact(limit: int, spec, segment_candidates: int, repeats: int):
    timings = []
    result = None
    stats = None
    for _ in range(repeats):
        gc.collect()
        start = time.perf_counter()
        count, checksum, last, current_stats = consume_primes_compact_spec(
            limit, spec, segment_candidates
        )
        timings.append(time.perf_counter() - start)
        current = (count, checksum, last)
        if result is None:
            result = current
            stats = current_stats
        elif current != result:
            raise AssertionError("non-deterministic compact sieve result")
    return result, stats, timings


def run_benchmark(
    limit: int,
    segment_bytes_values: tuple[int, ...],
    repeats: int,
    legacy_limit: int,
    legacy_segment_bytes: int,
):
    rows = []
    reference = None

    compact_specs = {}
    for modulus, ps in WHEELS.items():
        start = time.perf_counter()
        spec = make_compact_wheel(ps)
        build_seconds = time.perf_counter() - start
        compact_specs[modulus] = spec
        runtime_table_bytes = sys.getsizeof(spec.residues) + sys.getsizeof(spec.residue_rank)
        print(
            f"wheel={modulus:5d} phi={spec.phi:5d} density={spec.density:.9f} "
            f"compact_payload={spec.table_payload_bytes:,} B "
            f"compact_runtime_arrays={runtime_table_bytes:,} B "
            f"build={build_seconds:.6f}s"
        )

    print("\nCompact cache-size sweep")
    for modulus, ps in WHEELS.items():
        spec = compact_specs[modulus]
        for segment_bytes in segment_bytes_values:
            segment_candidates = segment_bytes * 8
            result, stats, timings = _time_compact(
                limit, spec, segment_candidates, repeats
            )
            if reference is None:
                reference = result
            elif result != reference:
                raise AssertionError(
                    f"wheel {modulus} result differs from benchmark reference"
                )

            row = {
                "limit": limit,
                "wheel": modulus,
                "phi": spec.phi,
                "density": spec.density,
                "segment_target_bytes": segment_bytes,
                "segment_candidates": segment_candidates,
                "actual_max_segment_bytes": stats.max_segment_storage_bytes,
                "segments": stats.segments_processed,
                "represented_candidates": stats.represented_candidates,
                "strike_attempts": stats.strike_attempts,
                "newly_removed": stats.newly_removed,
                "wheel_table_payload_bytes": stats.wheel_table_payload_bytes,
                "base_storage_bytes": stats.base_storage_bytes,
                "algorithmic_working_set_bytes": stats.algorithmic_working_set_bytes,
                "repeats": repeats,
                "time_median_seconds": statistics.median(timings),
                "time_min_seconds": min(timings),
                "time_max_seconds": max(timings),
                "count": result[0],
                "checksum": result[1],
                "last": result[2],
            }
            rows.append(row)
            print(
                f"wheel={modulus:5d} seg={segment_bytes // 1024:4d} KiB "
                f"segments={stats.segments_processed:4d} "
                f"working={stats.algorithmic_working_set_bytes / 1024:8.1f} KiB "
                f"median={row['time_median_seconds']:9.4f}s "
                f"min={row['time_min_seconds']:9.4f}s"
            )

    print("\nBest compact segment size per wheel (median time)")
    best_rows = []
    for modulus in WHEELS:
        candidates = [r for r in rows if r["wheel"] == modulus]
        best = min(candidates, key=lambda r: r["time_median_seconds"])
        best_rows.append(best)
        print(
            f"wheel={modulus:5d}: {best['segment_target_bytes'] // 1024:4d} KiB, "
            f"median={best['time_median_seconds']:.6f}s, "
            f"working={best['algorithmic_working_set_bytes'] / 1024:.1f} KiB"
        )

    print("\nLegacy vs compact table representation at a smaller comparison limit")
    legacy_rows = []
    for modulus, ps in WHEELS.items():
        legacy_spec = make_wheel(ps)
        seen = set()
        legacy_table_bytes = (
            _deep_size(legacy_spec.residues, seen)
            + _deep_size(legacy_spec.residue_rank, seen)
        )
        compact_spec = compact_specs[modulus]
        compact_runtime_bytes = (
            sys.getsizeof(compact_spec.residues)
            + sys.getsizeof(compact_spec.residue_rank)
        )
        segment_candidates = legacy_segment_bytes * 8

        gc.collect()
        start = time.perf_counter()
        old_result = consume_primes_primorial_wheel(
            legacy_limit, ps, segment_candidates
        )
        legacy_seconds = time.perf_counter() - start

        gc.collect()
        start = time.perf_counter()
        new_result = consume_primes_compact_spec(
            legacy_limit, compact_spec, segment_candidates
        )
        compact_seconds = time.perf_counter() - start

        if old_result[:3] != new_result[:3]:
            raise AssertionError(f"legacy/compact mismatch for wheel {modulus}")

        legacy_row = {
            "limit": legacy_limit,
            "wheel": modulus,
            "segment_bytes": legacy_segment_bytes,
            "legacy_table_deep_bytes": legacy_table_bytes,
            "compact_table_payload_bytes": compact_spec.table_payload_bytes,
            "compact_table_runtime_bytes": compact_runtime_bytes,
            "table_deep_to_payload_ratio": legacy_table_bytes / compact_spec.table_payload_bytes,
            "legacy_seconds": legacy_seconds,
            "compact_seconds": compact_seconds,
            "speedup_legacy_over_compact": legacy_seconds / compact_seconds,
        }
        legacy_rows.append(legacy_row)
        print(
            f"wheel={modulus:5d}: table {legacy_table_bytes:,} -> "
            f"{compact_spec.table_payload_bytes:,} B payload "
            f"({legacy_table_bytes / compact_spec.table_payload_bytes:.2f}x smaller); "
            f"time {legacy_seconds:.4f}s -> {compact_seconds:.4f}s "
            f"({legacy_seconds / compact_seconds:.3f}x)"
        )

    return rows, best_rows, legacy_rows


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("limit", type=int, nargs="?", default=100_000_000)
    parser.add_argument(
        "--segment-bytes",
        default=",".join(str(x) for x in DEFAULT_SEGMENT_BYTES),
        help="comma-separated target bitset sizes in bytes",
    )
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--legacy-limit", type=int, default=10_000_000)
    parser.add_argument("--legacy-segment-bytes", type=int, default=4096)
    parser.add_argument("--csv", type=Path, default=Path("compact-cache-sweep.csv"))
    parser.add_argument("--best-csv", type=Path, default=Path("compact-cache-best.csv"))
    parser.add_argument("--legacy-csv", type=Path, default=Path("compact-vs-legacy.csv"))
    args = parser.parse_args()

    if args.limit < 2 or args.legacy_limit < 2:
        raise ValueError("limits must be >= 2")
    if args.repeats < 1:
        raise ValueError("repeats must be >= 1")
    segment_bytes_values = tuple(
        int(x.strip()) for x in args.segment_bytes.split(",") if x.strip()
    )
    if not segment_bytes_values or any(x < 1 for x in segment_bytes_values):
        raise ValueError("segment byte values must be positive")

    rows, best_rows, legacy_rows = run_benchmark(
        args.limit,
        segment_bytes_values,
        args.repeats,
        args.legacy_limit,
        args.legacy_segment_bytes,
    )
    write_csv(args.csv, rows)
    write_csv(args.best_csv, best_rows)
    write_csv(args.legacy_csv, legacy_rows)


if __name__ == "__main__":
    main()

"""تحقق مستقل وموسع مناسب للتشغيل داخل GitHub Actions."""

from __future__ import annotations

import argparse
from math import isqrt

from optimized_method import primes_optimized
from original_method import primes_original
from retain_prime_compact_method import primes_retain_compact
from retain_prime_method import primes_retain
from retain_prime_packed_output_method import (
    iter_primes_retain_packed,
    primes_retain_packed,
)
from segmented_method import iter_primes_segmented, primes_segmented_packed
from trace_method import trace_primes
from wheel30_segmented_method import (
    iter_primes_wheel30_segmented,
    primes_wheel30_segmented_packed,
)


def reference_primes(limit: int) -> list[int]:
    result: list[int] = []
    for n in range(2, limit + 1):
        prime = True
        for d in range(2, isqrt(n) + 1):
            if n % d == 0:
                prime = False
                break
        if prime:
            result.append(n)
    return result


def verify(limit: int) -> None:
    expected = reference_primes(limit)

    results = {
        "original_method": primes_original(limit)[0],
        "retain_prime_method": primes_retain(limit)[0],
        "retain_prime_compact_method": primes_retain_compact(limit)[0],
        "retain_prime_packed_output_method": list(primes_retain_packed(limit)[0]),
        "retain_prime_stream": list(iter_primes_retain_packed(limit)),
        "segmented_stream": list(iter_primes_segmented(limit, 257)),
        "segmented_packed": list(primes_segmented_packed(limit, 257)[0]),
        "wheel30_segmented_stream": list(iter_primes_wheel30_segmented(limit, 257)),
        "wheel30_segmented_packed": list(
            primes_wheel30_segmented_packed(limit, 257)[0]
        ),
        "optimized_method": primes_optimized(limit)[0],
        "trace_method": trace_primes(limit)[0],
    }

    for name, actual in results.items():
        if actual != expected:
            raise AssertionError(
                f"{name} differs from independent reference at N={limit}: "
                f"expected {len(expected)} primes, got {len(actual)}"
            )

    print(f"Verified N={limit:,}")
    print(f"Prime count: {len(expected):,}")
    print("All project methods match the independent reference.")


def main() -> int:
    parser = argparse.ArgumentParser(description="Independent verification for CI")
    parser.add_argument("limit", type=int, nargs="?", default=4_999)
    args = parser.parse_args()

    if args.limit < 0:
        parser.error("limit must be >= 0")

    verify(args.limit)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

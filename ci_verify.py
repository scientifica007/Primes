"""تحقق مستقل وموسع مناسب للتشغيل داخل GitHub Actions.

يقارن تطبيقات المشروع مع مرجع بسيط يعتمد القسمة التجريبية حتى الجذر التربيعي.
لا يُستخدم هذا المرجع للأداء؛ دوره التحقق من صحة الناتج فقط.
"""

from __future__ import annotations

import argparse
from math import isqrt

from optimized_method import primes_optimized
from original_method import primes_original
from retain_prime_compact_method import primes_retain_compact
from retain_prime_method import primes_retain
from trace_method import trace_primes


def reference_primes(limit: int) -> list[int]:
    """إرجاع الأعداد الأولية <= limit بطريقة مستقلة وبسيطة."""
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

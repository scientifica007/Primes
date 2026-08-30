"""نسخة محسّنة من الفكرة نفسها للبحث عن الأعداد الأولية.

التحسينات الرئيسة:
- لا ننشئ قائمة الأعداد الطبيعية ولا قوائم منفصلة لكل المضاعفات.
- نستخدم مصفوفة منطقية للمرشحين.
- عند معالجة عدد أولي p نبدأ شطب مضاعفاته من p^2، لأن المضاعفات الأصغر
  سبق شطبها بواسطة عوامل أولية أصغر.
- نستمر حتى sqrt(N)، ثم تكون جميع القيم غير المشطوبة أولية.

هذه الصياغة هي عمليًا غربال إراتوستينس بكفاءة أعلى من النسخة المرجعية.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt
from typing import List, Tuple


@dataclass
class Stats:
    limit: int
    strike_attempts: int = 0
    newly_removed: int = 0
    processed_primes: int = 0


def primes_optimized(limit: int) -> Tuple[List[int], Stats]:
    """إرجاع جميع الأعداد الأولية <= limit باستخدام غربال محسّن."""
    if limit < 0:
        raise ValueError("limit must be >= 0")

    stats = Stats(limit=limit)
    if limit < 2:
        return [], stats

    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"

    # معالجة 2 منفردًا ثم الاقتصار على المرشحين الفرديين يقلل العمل.
    for multiple in range(4, limit + 1, 2):
        stats.strike_attempts += 1
        if is_prime[multiple]:
            is_prime[multiple] = 0
            stats.newly_removed += 1
    stats.processed_primes += 1

    root = isqrt(limit)
    for p in range(3, root + 1, 2):
        if not is_prime[p]:
            continue

        stats.processed_primes += 1
        # البداية من p^2 هي التحسين الحاسم.
        step = 2 * p
        for multiple in range(p * p, limit + 1, step):
            stats.strike_attempts += 1
            if is_prime[multiple]:
                is_prime[multiple] = 0
                stats.newly_removed += 1

    primes = [n for n in range(2, limit + 1) if is_prime[n]]
    return primes, stats


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Optimized prime-search experiment")
    parser.add_argument("limit", type=int, nargs="?", default=100)
    args = parser.parse_args()

    primes, stats = primes_optimized(args.limit)
    print(primes)
    print(stats)

"""نسخة مضغوطة من منهج إبقاء العدد الأولي وحذف المركبات فقط.

الهدف هو حل مشكلة الذاكرة في retain_prime_method.py دون تغيير الفكرة الرياضية:
- العدد الأولي يبقى في فضاء العمل ولا يُحذف.
- نحذف مضاعفاته المركبة بدءًا من p^2.
- العدد الباقي التالي هو العدد الأولي التالي.

بدل set، نمثل الأعداد الفردية فقط ببت واحد لكل مرشح داخل bytearray.
العدد 2 يعالج منفردًا، لأن كل عدد أولي آخر فردي.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple


@dataclass
class Stats:
    limit: int
    represented_odd_candidates: int = 0
    storage_bytes: int = 0
    processed_primes: int = 0
    strike_attempts: int = 0
    newly_removed: int = 0
    next_candidate_tests: int = 0


def primes_retain_compact(limit: int) -> Tuple[List[int], Stats]:
    """إرجاع جميع الأعداد الأولية <= limit باستخدام تمثيل بتّي للأعداد الفردية."""
    if limit < 0:
        raise ValueError("limit must be >= 0")

    stats = Stats(limit=limit)
    if limit < 2:
        return [], stats
    if limit == 2:
        return [2], stats

    # الفهرس i يمثل العدد الفردي n = 2*i + 3.
    odd_count = (limit - 1) // 2
    byte_count = (odd_count + 7) // 8
    alive = bytearray(b"\xff") * byte_count

    stats.represented_odd_candidates = odd_count
    stats.storage_bytes = byte_count

    p = 3
    p_index = 0

    while p * p <= limit:
        stats.processed_primes += 1

        # نمر فقط على المضاعفات الفردية: p^2, p^2 + 2p, ...
        for multiple in range(p * p, limit + 1, 2 * p):
            stats.strike_attempts += 1
            index = (multiple - 3) // 2
            byte_index = index >> 3
            mask = 1 << (index & 7)
            if alive[byte_index] & mask:
                alive[byte_index] &= ~mask
                stats.newly_removed += 1

        # أصغر مرشح فردي باقٍ بعد p هو العدد الأولي التالي.
        next_index = p_index + 1
        while next_index < odd_count:
            stats.next_candidate_tests += 1
            byte_index = next_index >> 3
            mask = 1 << (next_index & 7)
            if alive[byte_index] & mask:
                p_index = next_index
                p = 2 * next_index + 3
                break
            next_index += 1
        else:
            break

    primes: List[int] = [2]
    for index in range(odd_count):
        byte_index = index >> 3
        mask = 1 << (index & 7)
        if alive[byte_index] & mask:
            primes.append(2 * index + 3)

    return primes, stats


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Compact retained-prime search using one bit per odd candidate"
    )
    parser.add_argument("limit", type=int, nargs="?", default=100)
    args = parser.parse_args()

    primes, stats = primes_retain_compact(args.limit)
    print(primes)
    print(stats)

"""نسخة تجريبية تُبقي العدد الأولي وتحذف مضاعفاته المركبة فقط.

الفكرة:
- نبدأ بمجموعة الأعداد من 2 إلى N.
- نأخذ أصغر عنصر باقٍ p.
- نبقي p داخل المجموعة.
- نحذف مضاعفات p المركبة، ونبدأ من p^2 لأن كل مضاعف أصغر سبق أن
  كان قابلًا للحذف بواسطة عامل أولي أصغر.
- العدد الباقي التالي بعد p هو العدد الأولي التالي.
- عندما يصبح p^2 > N تكون كل العناصر الباقية أعدادًا أولية.

هذه النسخة تستخدم set في Python لتجسيد فكرة "الطرح مع إبقاء العدد الأولي"
بصورة مباشرة، ولذلك تختلف خصائص الذاكرة عن النسخة المحسنة المعتمدة على bytearray.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple


@dataclass
class Stats:
    limit: int
    generated_candidates: int = 0
    processed_primes: int = 0
    strike_attempts: int = 0
    newly_removed: int = 0
    next_candidate_tests: int = 0


def primes_retain(limit: int) -> Tuple[List[int], Stats]:
    """إرجاع جميع الأعداد الأولية <= limit مع إبقاء كل عدد أولي في مجموعة العمل."""
    if limit < 0:
        raise ValueError("limit must be >= 0")

    stats = Stats(limit=limit)
    if limit < 2:
        return [], stats

    survivors = set(range(2, limit + 1))
    stats.generated_candidates = len(survivors)

    p = 2
    while p * p <= limit:
        stats.processed_primes += 1

        # p نفسه لا يُحذف. نبدأ من p^2 لأن المضاعفات الأصغر سبق حذفها.
        for multiple in range(p * p, limit + 1, p):
            stats.strike_attempts += 1
            if multiple in survivors:
                survivors.remove(multiple)
                stats.newly_removed += 1

        # أصغر عنصر باقٍ بعد p هو العدد الأولي التالي.
        q = p + 1
        while q <= limit:
            stats.next_candidate_tests += 1
            if q in survivors:
                p = q
                break
            q += 1
        else:
            break

    return sorted(survivors), stats


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Prime-search experiment that keeps each prime and removes only its composite multiples"
    )
    parser.add_argument("limit", type=int, nargs="?", default=100)
    args = parser.parse_args()

    primes, stats = primes_retain(args.limit)
    print(primes)
    print(stats)

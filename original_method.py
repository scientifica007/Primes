"""التطبيق المرجعي للفكرة الأصلية للبحث التدريجي عن الأعداد الأولية.

الفكرة:
1. إنشاء قائمة الأعداد الطبيعية حتى حد N.
2. حذف مضاعفات 2، ثم مضاعفات 3.
3. بعد معالجة عدد أولي p، كل مرشح باقٍ يحقق p < n < p^2
   يُصنَّف أوليًا، مع تجنب إعادة إدراج ما سبق تصنيفه.
4. الانتقال إلى العدد الأولي التالي المكتشف وتكرار العملية.

هذه النسخة مقصودة للوضوح والتجربة، وليست مصممة لتحقيق أفضل أداء ممكن.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple


@dataclass
class Stats:
    limit: int
    generated_naturals: int = 0
    generated_multiples: int = 0
    removal_tests: int = 0
    removed_values: int = 0
    classification_tests: int = 0


def _subtract_multiples(
    candidates: List[int], p: int, limit: int, stats: Stats
) -> List[int]:
    """ينشئ قائمة مضاعفات p ثم يطرحها من قائمة المرشحين حرفيًا تقريبًا."""
    multiples = list(range(p, limit + 1, p))
    stats.generated_multiples += len(multiples)
    multiples_set = set(multiples)

    result: List[int] = []
    for n in candidates:
        stats.removal_tests += 1
        if n in multiples_set:
            stats.removed_values += 1
        else:
            result.append(n)
    return result


def primes_original(limit: int) -> Tuple[List[int], Stats]:
    """إرجاع الأعداد الأولية <= limit وفق المنهج الأصلي التجريبي."""
    if limit < 0:
        raise ValueError("limit must be >= 0")

    stats = Stats(limit=limit)
    if limit < 2:
        return [], stats

    naturals = list(range(1, limit + 1))
    stats.generated_naturals = len(naturals)

    primes: List[int] = [2]
    prime_set = {2}

    # القائمة أ: طرح الأعداد الزوجية من الأعداد الطبيعية.
    candidates = _subtract_multiples(naturals, 2, limit, stats)
    # العدد 1 ليس أوليًا ولا مركبًا؛ نستبعده من فضاء المرشحين.
    candidates = [n for n in candidates if n != 1]

    if limit < 3:
        return primes, stats

    primes.append(3)
    prime_set.add(3)

    # القائمة ب: طرح مضاعفات 3 من القائمة أ.
    candidates = _subtract_multiples(candidates, 3, limit, stats)

    # بعد معالجة 3: الأعداد الباقية بين 3 و9 أولية.
    classified_upper = 3
    upper = min(limit, 3 * 3 - 1)
    for n in candidates:
        if n > 3 and n <= upper:
            stats.classification_tests += 1
            if n not in prime_set:
                primes.append(n)
                prime_set.add(n)
    classified_upper = upper

    # نبدأ بالعدد الأولي التالي بعد 3، أي 5 عندما يكون ضمن المجال.
    process_index = 2

    while classified_upper < limit:
        # إذا لم يكن لدينا بعدُ عدد أولي تالٍ، فأصغر مرشح غير مصنف هو أولي.
        if process_index >= len(primes):
            next_candidates = [n for n in candidates if n > classified_upper]
            if not next_candidates:
                break
            q = min(next_candidates)
            primes.append(q)
            prime_set.add(q)

        p = primes[process_index]
        candidates = _subtract_multiples(candidates, p, limit, stats)

        upper = min(limit, p * p - 1)
        if upper <= classified_upper:
            process_index += 1
            continue

        for n in candidates:
            if n > p and classified_upper < n <= upper:
                stats.classification_tests += 1
                if n not in prime_set:
                    primes.append(n)
                    prime_set.add(n)

        classified_upper = upper
        process_index += 1

    primes = sorted(p for p in prime_set if p <= limit)
    return primes, stats


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Original incremental prime-search experiment")
    parser.add_argument("limit", type=int, nargs="?", default=100)
    args = parser.parse_args()

    primes, stats = primes_original(args.limit)
    print(primes)
    print(stats)

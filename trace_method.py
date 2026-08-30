"""عرض تعليمي خطوة بخطوة لتشكل قوائم تجربة البحث عن الأعداد الأولية.

هذا الملف لا يهدف إلى السرعة، بل إلى جعل كل مرحلة قابلة للمشاهدة:
- قائمة الأعداد الطبيعية.
- قائمة الأعداد الزوجية.
- القائمة أ الناتجة عن طرح الزوجيات.
- القائمة ب بعد طرح مضاعفات 3.
- القوائم ج، د، هـ... بعد طرح مضاعفات الأعداد الأولية التالية.
- الأعداد المحذوفة فعليًا في كل مرحلة.
- الأعداد الأولية الجديدة المصنفة في كل مرحلة.
- قائمة الأعداد الأولية المتراكمة.

تشغيل مثال:
    python trace_method.py 100
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import List


@dataclass
class Stage:
    label: str
    prime: int
    square: int
    multiples: List[int]
    removed: List[int]
    remaining: List[int]
    new_primes: List[int]
    accumulated_primes: List[int]


def _fmt(values: List[int]) -> str:
    return "[" + ", ".join(map(str, values)) + "]"


def trace_primes(limit: int) -> tuple[List[int], List[Stage], List[int], List[int], List[int]]:
    """يبني التتبع الكامل للطريقة الأصلية حتى ``limit``.

    يعاد:
        primes: قائمة الأعداد الأولية المكتشفة.
        stages: مراحل الطرح والتصنيف.
        naturals: الأعداد الطبيعية من 1 إلى limit.
        evens: الأعداد الزوجية حتى limit.
        list_a: القائمة أ = naturals - evens، وتحتوي 1 كما في التعريف الأصلي.

    ملاحظة: العدد 1 يظهر في القوائم الناتجة لأنه ليس من مضاعفات أي عدد أولي،
    لكنه لا يدخل مطلقًا في قائمة الأعداد الأولية.
    """
    if limit < 0:
        raise ValueError("limit must be >= 0")

    naturals = list(range(1, limit + 1))
    evens = list(range(2, limit + 1, 2))
    even_set = set(evens)
    list_a = [n for n in naturals if n not in even_set]

    primes: List[int] = []
    prime_set: set[int] = set()
    stages: List[Stage] = []

    if limit >= 2:
        primes.append(2)
        prime_set.add(2)

    if limit < 3:
        return primes, stages, naturals, evens, list_a

    primes.append(3)
    prime_set.add(3)

    # القائمة ب: طرح مضاعفات 3 من القائمة أ.
    candidates = list_a[:]
    multiples_3 = list(range(3, limit + 1, 3))
    multiples_set = set(multiples_3)
    removed = [n for n in candidates if n in multiples_set]
    candidates = [n for n in candidates if n not in multiples_set]

    upper = min(limit, 3 * 3 - 1)
    new_primes = [n for n in candidates if 3 < n <= upper]
    for n in new_primes:
        if n not in prime_set:
            primes.append(n)
            prime_set.add(n)

    stages.append(
        Stage(
            label="ب",
            prime=3,
            square=9,
            multiples=multiples_3,
            removed=removed,
            remaining=candidates[:],
            new_primes=new_primes,
            accumulated_primes=sorted(prime_set),
        )
    )

    classified_upper = upper
    process_index = 2
    labels = ["ج", "د", "هـ", "و", "ز", "ح", "ط", "ي", "ك", "ل", "م", "ن"]
    label_index = 0

    while classified_upper < limit:
        # إذا لم يكن العدد الأولي التالي قد صُنّف بعد، فأصغر مرشح غير مصنف هو التالي.
        if process_index >= len(primes):
            future = [n for n in candidates if n > classified_upper]
            if not future:
                break
            q = min(future)
            primes.append(q)
            prime_set.add(q)

        p = primes[process_index]
        multiples = list(range(p, limit + 1, p))
        multiples_set = set(multiples)
        removed = [n for n in candidates if n in multiples_set]
        candidates = [n for n in candidates if n not in multiples_set]

        upper = min(limit, p * p - 1)
        new_primes = [
            n
            for n in candidates
            if n > p and classified_upper < n <= upper and n not in prime_set
        ]
        for n in new_primes:
            primes.append(n)
            prime_set.add(n)

        label = labels[label_index] if label_index < len(labels) else f"مرحلة-{label_index + 3}"
        stages.append(
            Stage(
                label=label,
                prime=p,
                square=p * p,
                multiples=multiples,
                removed=removed,
                remaining=candidates[:],
                new_primes=new_primes,
                accumulated_primes=sorted(prime_set),
            )
        )

        classified_upper = upper
        process_index += 1
        label_index += 1

    return sorted(n for n in prime_set if n <= limit), stages, naturals, evens, list_a


def print_trace(limit: int) -> None:
    primes, stages, naturals, evens, list_a = trace_primes(limit)

    print(f"الحد الأعلى N = {limit}")
    print()
    print("قائمة الأعداد الطبيعية:")
    print(_fmt(naturals))
    print()
    print("قائمة الأعداد الزوجية:")
    print(_fmt(evens))
    print()
    print("القائمة أ = الأعداد الطبيعية - الأعداد الزوجية:")
    print(_fmt(list_a))
    print("ملاحظة: 1 يبقى ظاهرًا في القوائم، لكنه لا يصنف عددًا أوليًا.")
    print()
    print("قائمة الأعداد الأولية الابتدائية:")
    print(_fmt([p for p in [2, 3] if p <= limit]))

    for stage in stages:
        print()
        print("=" * 72)
        print(f"القائمة {stage.label}: معالجة العدد الأولي p = {stage.prime}")
        print(f"p² = {stage.square}")
        print(f"مضاعفات {stage.prime} حتى {limit}:")
        print(_fmt(stage.multiples))
        print("العناصر التي حذفت فعليًا من القائمة السابقة:")
        print(_fmt(stage.removed))
        print(f"القائمة {stage.label} بعد الطرح:")
        print(_fmt(stage.remaining))
        print("الأعداد الأولية الجديدة المكتشفة في هذه المرحلة:")
        print(_fmt(stage.new_primes))
        print("قائمة الأعداد الأولية المتراكمة:")
        print(_fmt(stage.accumulated_primes))

    print()
    print("=" * 72)
    print(f"القائمة النهائية للأعداد الأولية <= {limit}:")
    print(_fmt(primes))
    print(f"عددها = {len(primes)}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Step-by-step trace of the prime-list experiment")
    parser.add_argument("limit", type=int, nargs="?", default=100)
    args = parser.parse_args()
    if args.limit < 0:
        parser.error("limit must be >= 0")
    print_trace(args.limit)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

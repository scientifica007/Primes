#!/usr/bin/env python3
"""Generate a lossless normalized dataset for prime intervals around prime squares.

For each of the first 1000 primes p_n, the requested intervals are:
  1) primes strictly between p_n and p_n^2
  2) primes strictly between p_n^2 and p_{n+1}^2

Writing the first list literally inside every one of the 1000 rows would repeat
1,178,547,264 prime values. Instead, this script stores every prime up to
p_1001^2 once in a master CSV and stores exact 1-based start/end indices for
both requested lists in a compact 1000-row index CSV. Thus every list is
recoverable exactly without duplication.
"""

from bisect import bisect_left, bisect_right
import csv
from pathlib import Path

OUT = Path(__file__).resolve().parent
LIMIT = 7927 ** 2  # p_1001^2


def primes_up_to(n: int) -> list[int]:
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    r = int(n ** 0.5)
    for p in range(2, r + 1):
        if sieve[p]:
            start = p * p
            count = ((n - start) // p) + 1
            sieve[start:n + 1:p] = b"\x00" * count
    return [i for i, flag in enumerate(sieve) if flag]


def main() -> None:
    primes = primes_up_to(LIMIT)
    first_1001 = primes[:1001]
    assert first_1001[999] == 7919
    assert first_1001[1000] == 7927

    master = OUT / "primes_up_to_1001st_prime_square.csv"
    with master.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["prime_index_1_based", "prime"])
        for idx, q in enumerate(primes, start=1):
            w.writerow([idx, q])

    index_path = OUT / "first_1000_prime_square_interval_index.csv"
    headers = [
        "الترتيب_n", "p_n", "p_n^2", "p_التالي", "p_التالي^2",
        "فهرس_أول_أولي_بين_p_و_p2", "فهرس_آخر_أولي_بين_p_و_p2",
        "عدد_الأوليات_بين_p_و_p2", "أول_أولي_بين_p_و_p2", "آخر_أولي_بين_p_و_p2",
        "فهرس_أول_أولي_بين_p2_ومربع_التالي", "فهرس_آخر_أولي_بين_p2_ومربع_التالي",
        "عدد_الأوليات_بين_p2_ومربع_التالي", "أول_أولي_بين_p2_ومربع_التالي",
        "آخر_أولي_بين_p2_ومربع_التالي",
    ]

    repeated_first = 0
    second_total = 0
    rows = []
    for i in range(1000):
        p = first_1001[i]
        pn = first_1001[i + 1]
        p2, pn2 = p * p, pn * pn

        a = bisect_right(primes, p)
        b = bisect_left(primes, p2)
        c = bisect_right(primes, p2)
        d = bisect_left(primes, pn2)
        n1, n2 = b - a, d - c
        repeated_first += n1
        second_total += n2

        rows.append([
            i + 1, p, p2, pn, pn2,
            a + 1 if n1 else "", b if n1 else "", n1,
            primes[a] if n1 else "", primes[b - 1] if n1 else "",
            c + 1 if n2 else "", d if n2 else "", n2,
            primes[c] if n2 else "", primes[d - 1] if n2 else "",
        ])

    with index_path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow(headers)
        w.writerows(rows)

    preview = OUT / "first_20_prime_square_interval_lists_preview.csv"
    with preview.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "الترتيب_n", "p_n", "p_n^2", "قائمة_الأوليات_بين_p_و_p2",
            "p_التالي", "p_التالي^2", "قائمة_الأوليات_بين_p2_ومربع_التالي",
        ])
        for i in range(20):
            p, pn = first_1001[i], first_1001[i + 1]
            p2, pn2 = p * p, pn * pn
            a, b = bisect_right(primes, p), bisect_left(primes, p2)
            c, d = bisect_right(primes, p2), bisect_left(primes, pn2)
            w.writerow([
                i + 1, p, p2, "، ".join(map(str, primes[a:b])),
                pn, pn2, "، ".join(map(str, primes[c:d])),
            ])

    note = OUT / "PRIME_SQUARE_INTERVAL_DATASET.md"
    note.write_text(f"""# بيانات الفترات حول مربعات الأعداد الأولية\n\nتغطي البيانات أول 1000 عدد أولي، من 2 إلى 7919، وتستخدم العدد الأولي رقم 1001 وهو 7927 لتحديد نهاية الفترة الأخيرة.\n\n## لماذا استعملنا تمثيلا مفهرسا؟\n\nلو كتبنا قائمة الأوليات بين `p_n` و `p_n^2` حرفيا داخل كل صف، لتكرر في الجدول **{repeated_first:,}** عدد أولي. هذا تضخيم هائل لأن الفترات متداخلة بقوة.\n\nلذلك تحفظ كل الأوليات حتى `7927^2 = {LIMIT:,}` مرة واحدة في:\n\n`primes_up_to_1001st_prime_square.csv`\n\nثم يحفظ الجدول ذو 1000 صف في:\n\n`first_1000_prime_square_interval_index.csv`\n\nفهرس البداية والنهاية (1-based) لكل من القائمتين المطلوبتين. هذا تمثيل **دقيق وفاقده صفر**: كل قائمة يمكن استخراجها كما هي من الملف الرئيسي دون تكرار البيانات.\n\nعدد الأوليات في الملف الرئيسي: **{len(primes):,}**.\n\nومجموع الأوليات في الفترات غير المتداخلة بين `p_n^2` و `p_(n+1)^2` هو **{second_total:,}**.\n\nأضفنا أيضا ملف معاينة لأول 20 صفا بالقوائم مكتوبة حرفيا:\n\n`first_20_prime_square_interval_lists_preview.csv`\n\n## تعريف الفترات\n\nكلتا الفترتين **مفتوحتان** عند الطرفين:\n\n- `p_n < q < p_n^2`\n- `p_n^2 < q < p_(n+1)^2`\n\nحيث `q` عدد أولي.\n""", encoding="utf-8")

    print("Generated:", master.name, index_path.name, preview.name, note.name)
    print("Master prime count:", len(primes))
    print("Repeated literal first-interval entries:", repeated_first)
    print("Square-to-next-square total:", second_total)


if __name__ == "__main__":
    main()

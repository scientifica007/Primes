"""ضغط فضاء المرشحين وقائمة الناتج معًا مع إبقاء العدد الأولي.

هذه النسخة تبني فضاء المرشحين ببت واحد لكل عدد فردي، ثم لا تحول الناتج إلى
list من كائنات int في Python. بدل ذلك تحفظ الأعداد الأولية في array من أعداد
صحيحة غير موقعة ذات 32 بت عندما يكفي المجال، وتنتقل إلى 64 بت عند الحاجة.

كما توفر iter_primes_retain_packed واجهة streaming لا تخزن قائمة الناتج أصلًا.
"""

from __future__ import annotations

from array import array
from dataclasses import dataclass
from typing import Iterator, Tuple


@dataclass
class Stats:
    limit: int
    represented_odd_candidates: int = 0
    candidate_storage_bytes: int = 0
    processed_primes: int = 0
    strike_attempts: int = 0
    newly_removed: int = 0
    next_candidate_tests: int = 0
    output_typecode: str = ""
    output_itemsize: int = 0
    output_storage_bytes: int = 0


def _build_survivor_bits(limit: int, stats: Stats) -> Tuple[bytearray, int]:
    """بناء بتات الأعداد الفردية الباقية بعد حذف المركبات."""
    odd_count = (limit - 1) // 2
    byte_count = (odd_count + 7) // 8
    alive = bytearray(b"\xff") * byte_count

    stats.represented_odd_candidates = odd_count
    stats.candidate_storage_bytes = byte_count

    p = 3
    p_index = 0

    while p * p <= limit:
        stats.processed_primes += 1

        for multiple in range(p * p, limit + 1, 2 * p):
            stats.strike_attempts += 1
            index = (multiple - 3) // 2
            byte_index = index >> 3
            mask = 1 << (index & 7)
            if alive[byte_index] & mask:
                alive[byte_index] &= ~mask
                stats.newly_removed += 1

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

    return alive, odd_count


def _iter_survivors(limit: int, alive: bytearray, odd_count: int) -> Iterator[int]:
    if limit >= 2:
        yield 2

    for index in range(odd_count):
        byte_index = index >> 3
        mask = 1 << (index & 7)
        if alive[byte_index] & mask:
            n = 2 * index + 3
            if n <= limit:
                yield n


def _output_typecode(limit: int) -> str:
    """اختيار أصغر نوع unsigned متاح يكفي لقيمة limit."""
    uint_itemsize = array("I").itemsize
    uint_max = (1 << (8 * uint_itemsize)) - 1
    if limit <= uint_max:
        return "I"

    ull_itemsize = array("Q").itemsize
    ull_max = (1 << (8 * ull_itemsize)) - 1
    if limit <= ull_max:
        return "Q"

    raise OverflowError("limit is too large for packed unsigned integer output")


def primes_retain_packed(limit: int) -> Tuple[array, Stats]:
    """إرجاع الأوليات <= limit داخل array مضغوطة بدل list من int."""
    if limit < 0:
        raise ValueError("limit must be >= 0")

    stats = Stats(limit=limit)
    typecode = _output_typecode(limit)
    primes = array(typecode)

    stats.output_typecode = typecode
    stats.output_itemsize = primes.itemsize

    if limit < 2:
        return primes, stats

    if limit == 2:
        primes.append(2)
        stats.output_storage_bytes = len(primes) * primes.itemsize
        return primes, stats

    alive, odd_count = _build_survivor_bits(limit, stats)
    primes.extend(_iter_survivors(limit, alive, odd_count))
    stats.output_storage_bytes = len(primes) * primes.itemsize
    return primes, stats


def iter_primes_retain_packed(limit: int) -> Iterator[int]:
    """Streaming: توليد الأوليات واحدةً تلو الأخرى دون تخزين قائمة الناتج."""
    if limit < 0:
        raise ValueError("limit must be >= 0")
    if limit < 2:
        return
    if limit == 2:
        yield 2
        return

    stats = Stats(limit=limit)
    alive, odd_count = _build_survivor_bits(limit, stats)
    yield from _iter_survivors(limit, alive, odd_count)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Retained-prime search with bit-packed candidates and packed output"
    )
    parser.add_argument("limit", type=int, nargs="?", default=100)
    parser.add_argument(
        "--stream",
        action="store_true",
        help="Print primes from the streaming interface instead of the packed array",
    )
    args = parser.parse_args()

    if args.stream:
        print(list(iter_primes_retain_packed(args.limit)))
    else:
        primes, stats = primes_retain_packed(args.limit)
        print(list(primes))
        print(stats)

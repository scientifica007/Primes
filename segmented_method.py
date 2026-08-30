"""غربلة مجزأة بتمثيل بتّي مع إبقاء الأعداد الأولية.

الهدف هو عدم تخصيص بت لكل عدد فردي حتى N دفعةً واحدة. بدل ذلك نقسم المجال
إلى مقاطع ثابتة السعة، نمثل كل مقطع ببت واحد لكل عدد فردي، نحذف مركباته
باستخدام الأعداد الأولية الأساسية حتى sqrt(N)، ثم نخرج الأوليات ونتخلص من
المقطع قبل الانتقال إلى المقطع التالي.

هذا يجعل ذاكرة فضاء العمل مرتبطة بحجم المقطع أساسًا، مع كلفة إضافية صغيرة
لقائمة الأعداد الأولية الأساسية حتى sqrt(N).
"""

from __future__ import annotations

from array import array
from dataclasses import dataclass
from math import isqrt
from typing import Iterator, Tuple


DEFAULT_SEGMENT_ODDS = 32_768


@dataclass
class Stats:
    limit: int
    segment_odd_capacity: int
    base_limit: int = 0
    base_prime_count: int = 0
    base_candidate_storage_bytes: int = 0
    base_output_storage_bytes: int = 0
    segments_processed: int = 0
    max_segment_odd_count: int = 0
    max_segment_storage_bytes: int = 0
    strike_attempts: int = 0
    newly_removed: int = 0
    yielded_primes: int = 0
    output_typecode: str = ""
    output_itemsize: int = 0
    output_storage_bytes: int = 0


def _is_alive(bits: bytearray, index: int) -> bool:
    return bool(bits[index >> 3] & (1 << (index & 7)))


def _clear(bits: bytearray, index: int) -> bool:
    """مسح البت وإرجاع True فقط إذا كان حيًا قبل المسح."""
    byte_index = index >> 3
    mask = 1 << (index & 7)
    if bits[byte_index] & mask:
        bits[byte_index] &= ~mask
        return True
    return False


def _typecode_for_value(value: int) -> str:
    uint_itemsize = array("I").itemsize
    if value <= (1 << (8 * uint_itemsize)) - 1:
        return "I"
    ull_itemsize = array("Q").itemsize
    if value <= (1 << (8 * ull_itemsize)) - 1:
        return "Q"
    raise OverflowError("value is too large for packed unsigned integer storage")


def _base_primes(root: int) -> Tuple[array, int]:
    """بناء الأوليات <= root بتمثيل بتّي صغير للأعداد الفردية."""
    typecode = _typecode_for_value(max(root, 2))
    result = array(typecode)
    if root < 2:
        return result, 0

    result.append(2)
    if root < 3:
        return result, 0

    odd_count = (root - 1) // 2  # 3,5,7,...
    byte_count = (odd_count + 7) // 8
    alive = bytearray(b"\xff") * byte_count

    p = 3
    p_index = 0
    while p * p <= root:
        for multiple in range(p * p, root + 1, 2 * p):
            index = (multiple - 3) // 2
            _clear(alive, index)

        next_index = p_index + 1
        while next_index < odd_count:
            if _is_alive(alive, next_index):
                p_index = next_index
                p = 2 * next_index + 3
                break
            next_index += 1
        else:
            break

    for index in range(odd_count):
        if _is_alive(alive, index):
            n = 2 * index + 3
            if n <= root:
                result.append(n)

    return result, byte_count


def _iter_segmented(limit: int, segment_odd_count: int, stats: Stats) -> Iterator[int]:
    if limit >= 2:
        stats.yielded_primes += 1
        yield 2
    if limit < 3:
        return

    root = isqrt(limit)
    base_primes, base_candidate_bytes = _base_primes(root)
    stats.base_limit = root
    stats.base_prime_count = len(base_primes)
    stats.base_candidate_storage_bytes = base_candidate_bytes
    stats.base_output_storage_bytes = len(base_primes) * base_primes.itemsize

    low = 3
    while low <= limit:
        # low و high فرديان، والمقطع يمثل الأعداد low, low+2, ... high.
        high = min(limit, low + 2 * (segment_odd_count - 1))
        if high % 2 == 0:
            high -= 1
        if high < low:
            break

        count = ((high - low) // 2) + 1
        byte_count = (count + 7) // 8
        alive = bytearray(b"\xff") * byte_count

        stats.segments_processed += 1
        stats.max_segment_odd_count = max(stats.max_segment_odd_count, count)
        stats.max_segment_storage_bytes = max(stats.max_segment_storage_bytes, byte_count)

        for p in base_primes:
            if p == 2:
                continue
            if p * p > high:
                break

            start = max(p * p, ((low + p - 1) // p) * p)
            if start % 2 == 0:
                start += p

            for multiple in range(start, high + 1, 2 * p):
                stats.strike_attempts += 1
                index = (multiple - low) // 2
                if _clear(alive, index):
                    stats.newly_removed += 1

        for index in range(count):
            if _is_alive(alive, index):
                n = low + 2 * index
                if n <= limit:
                    stats.yielded_primes += 1
                    yield n

        low = high + 2


def iter_primes_segmented(
    limit: int, segment_odd_count: int = DEFAULT_SEGMENT_ODDS
) -> Iterator[int]:
    """توليد الأوليات <= limit مقطعًا بعد مقطع دون حفظ قائمة الناتج."""
    if limit < 0:
        raise ValueError("limit must be >= 0")
    if segment_odd_count < 1:
        raise ValueError("segment_odd_count must be >= 1")

    stats = Stats(limit=limit, segment_odd_capacity=segment_odd_count)
    yield from _iter_segmented(limit, segment_odd_count, stats)


def consume_primes_segmented(
    limit: int, segment_odd_count: int = DEFAULT_SEGMENT_ODDS
) -> Tuple[int, int, int, Stats]:
    """استهلاك streaming وإرجاع العدد والمجموع وآخر أولي والإحصاءات."""
    if limit < 0:
        raise ValueError("limit must be >= 0")
    if segment_odd_count < 1:
        raise ValueError("segment_odd_count must be >= 1")

    stats = Stats(limit=limit, segment_odd_capacity=segment_odd_count)
    count = 0
    checksum = 0
    last = 0
    for prime in _iter_segmented(limit, segment_odd_count, stats):
        count += 1
        checksum += prime
        last = prime
    return count, checksum, last, stats


def primes_segmented_packed(
    limit: int, segment_odd_count: int = DEFAULT_SEGMENT_ODDS
) -> Tuple[array, Stats]:
    """نسخة مجزأة تحفظ الناتج النهائي في array مضغوطة عند الحاجة إليه كاملًا."""
    if limit < 0:
        raise ValueError("limit must be >= 0")
    if segment_odd_count < 1:
        raise ValueError("segment_odd_count must be >= 1")

    typecode = _typecode_for_value(max(limit, 2))
    output = array(typecode)
    stats = Stats(limit=limit, segment_odd_capacity=segment_odd_count)
    stats.output_typecode = typecode
    stats.output_itemsize = output.itemsize

    output.extend(_iter_segmented(limit, segment_odd_count, stats))
    stats.output_storage_bytes = len(output) * output.itemsize
    return output, stats


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Segmented retained-prime search")
    parser.add_argument("limit", type=int, nargs="?", default=100)
    parser.add_argument(
        "--segment-odds",
        type=int,
        default=DEFAULT_SEGMENT_ODDS,
        help="Number of odd candidates represented per segment",
    )
    parser.add_argument(
        "--packed",
        action="store_true",
        help="Store all output in a packed array instead of pure streaming",
    )
    args = parser.parse_args()

    if args.packed:
        primes, stats = primes_segmented_packed(args.limit, args.segment_odds)
        print(list(primes))
        print(stats)
    else:
        count, checksum, last, stats = consume_primes_segmented(
            args.limit, args.segment_odds
        )
        print(f"count={count}")
        print(f"checksum={checksum}")
        print(f"last={last}")
        print(stats)

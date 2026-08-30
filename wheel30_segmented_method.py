"""غربلة مجزأة بعجلة 30 مع تمثيل بتّي وstreaming.

الفكرة: بدل تمثيل كل الأعداد الفردية، نمثل فقط الأعداد الموافقة للبواقي
1, 7, 11, 13, 17, 19, 23, 29 modulo 30. هذه هي الأعداد غير القابلة
للقسمة على 2 أو 3 أو 5، ولذلك تحتوي جميع الأعداد الأولية الأكبر من 5.

يقسم المجال إلى مقاطع ثابتة بعدد المرشحين (وليس بعدد الأعداد الصحيحة)،
ويستخدم أوليات الأساس حتى sqrt(N) لشطب المركبات. يبدأ شطب p من p^2،
ويبقى p نفسه ضمن المرشحين.
"""

from __future__ import annotations

from array import array
from dataclasses import dataclass
from math import isqrt
from typing import Iterator, Tuple


WHEEL = 30
RESIDUES = (1, 7, 11, 13, 17, 19, 23, 29)
RESIDUE_RANK = [-1] * WHEEL
for _rank, _residue in enumerate(RESIDUES):
    RESIDUE_RANK[_residue] = _rank

DEFAULT_SEGMENT_CANDIDATES = 32_768


@dataclass
class Stats:
    limit: int
    segment_candidate_capacity: int
    base_limit: int = 0
    base_prime_count: int = 0
    base_candidate_storage_bytes: int = 0
    base_output_storage_bytes: int = 0
    segments_processed: int = 0
    represented_candidates: int = 0
    max_segment_candidate_count: int = 0
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


def _candidate_value(index: int) -> int:
    """قيمة المرشح ذي الفهرس العالمي ضمن تسلسل Wheel-30؛ الفهرس 0 هو 1."""
    block, rank = divmod(index, len(RESIDUES))
    return block * WHEEL + RESIDUES[rank]


def _candidate_index(value: int) -> int:
    """فهرس قيمة معروفة بأنها coprime مع 30."""
    block, residue = divmod(value, WHEEL)
    rank = RESIDUE_RANK[residue]
    if rank < 0:
        raise ValueError(f"{value} is not a Wheel-30 candidate")
    return block * len(RESIDUES) + rank


def _first_candidate_index_at_least(value: int) -> int:
    block, residue = divmod(value, WHEEL)
    for rank, candidate_residue in enumerate(RESIDUES):
        if candidate_residue >= residue:
            return block * len(RESIDUES) + rank
    return (block + 1) * len(RESIDUES)


def _last_candidate_index_at_most(value: int) -> int:
    block, residue = divmod(value, WHEEL)
    for rank in range(len(RESIDUES) - 1, -1, -1):
        if RESIDUES[rank] <= residue:
            return block * len(RESIDUES) + rank
    return (block - 1) * len(RESIDUES) + (len(RESIDUES) - 1)


def _base_primes(root: int) -> Tuple[array, int]:
    """أوليات <= root باستخدام bitset فردي صغير؛ هذا ليس فضاء البحث الرئيسي."""
    typecode = _typecode_for_value(max(root, 2))
    result = array(typecode)
    if root < 2:
        return result, 0

    result.append(2)
    if root < 3:
        return result, 0

    odd_count = (root - 1) // 2
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


def _iter_wheel30_segmented(
    limit: int, segment_candidate_count: int, stats: Stats
) -> Iterator[int]:
    # 2 و3 و5 لا تظهر في عجلة 30، لذلك نخرجها صراحة.
    for prime in (2, 3, 5):
        if prime <= limit:
            stats.yielded_primes += 1
            yield prime

    if limit < 7:
        return

    root = isqrt(limit)
    base_primes, base_candidate_bytes = _base_primes(root)
    stats.base_limit = root
    stats.base_prime_count = len(base_primes)
    stats.base_candidate_storage_bytes = base_candidate_bytes
    stats.base_output_storage_bytes = len(base_primes) * base_primes.itemsize

    # الفهرس 0 يمثل العدد 1؛ أول مرشح أولي فعلي في العجلة هو 7 عند الفهرس 1.
    start_index = 1
    last_index = _last_candidate_index_at_most(limit)

    while start_index <= last_index:
        end_index = min(last_index, start_index + segment_candidate_count - 1)
        count = end_index - start_index + 1
        low = _candidate_value(start_index)
        high = _candidate_value(end_index)

        byte_count = (count + 7) // 8
        alive = bytearray(b"\xff") * byte_count

        stats.segments_processed += 1
        stats.represented_candidates += count
        stats.max_segment_candidate_count = max(stats.max_segment_candidate_count, count)
        stats.max_segment_storage_bytes = max(stats.max_segment_storage_bytes, byte_count)

        for p in base_primes:
            if p < 7:
                continue
            if p * p > high:
                break

            # m = p*k. لكي يكون m ضمن Wheel-30، يجب أن يكون k أيضًا coprime مع 30.
            k_min = max(p, (low + p - 1) // p)
            k_index = _first_candidate_index_at_least(k_min)

            while True:
                k = _candidate_value(k_index)
                multiple = p * k
                if multiple > high:
                    break

                stats.strike_attempts += 1
                multiple_index = _candidate_index(multiple)
                local_index = multiple_index - start_index
                if _clear(alive, local_index):
                    stats.newly_removed += 1
                k_index += 1

        for local_index in range(count):
            if _is_alive(alive, local_index):
                n = _candidate_value(start_index + local_index)
                if n <= limit:
                    stats.yielded_primes += 1
                    yield n

        start_index = end_index + 1


def iter_primes_wheel30_segmented(
    limit: int, segment_candidate_count: int = DEFAULT_SEGMENT_CANDIDATES
) -> Iterator[int]:
    """توليد أوليات <= limit عبر Wheel-30 مجزأة دون حفظ الناتج كاملًا."""
    if limit < 0:
        raise ValueError("limit must be >= 0")
    if segment_candidate_count < 1:
        raise ValueError("segment_candidate_count must be >= 1")

    stats = Stats(limit=limit, segment_candidate_capacity=segment_candidate_count)
    yield from _iter_wheel30_segmented(limit, segment_candidate_count, stats)


def consume_primes_wheel30_segmented(
    limit: int, segment_candidate_count: int = DEFAULT_SEGMENT_CANDIDATES
) -> Tuple[int, int, int, Stats]:
    """استهلاك streaming وإرجاع العدد والمجموع وآخر أولي والإحصاءات."""
    if limit < 0:
        raise ValueError("limit must be >= 0")
    if segment_candidate_count < 1:
        raise ValueError("segment_candidate_count must be >= 1")

    stats = Stats(limit=limit, segment_candidate_capacity=segment_candidate_count)
    count = 0
    checksum = 0
    last = 0
    for prime in _iter_wheel30_segmented(limit, segment_candidate_count, stats):
        count += 1
        checksum += prime
        last = prime
    return count, checksum, last, stats


def primes_wheel30_segmented_packed(
    limit: int, segment_candidate_count: int = DEFAULT_SEGMENT_CANDIDATES
) -> Tuple[array, Stats]:
    """نسخة Wheel-30 المجزأة مع حفظ الناتج النهائي في array مضغوطة."""
    if limit < 0:
        raise ValueError("limit must be >= 0")
    if segment_candidate_count < 1:
        raise ValueError("segment_candidate_count must be >= 1")

    typecode = _typecode_for_value(max(limit, 2))
    output = array(typecode)
    stats = Stats(limit=limit, segment_candidate_capacity=segment_candidate_count)
    stats.output_typecode = typecode
    stats.output_itemsize = output.itemsize

    output.extend(_iter_wheel30_segmented(limit, segment_candidate_count, stats))
    stats.output_storage_bytes = len(output) * output.itemsize
    return output, stats


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Wheel-30 segmented retained-prime search")
    parser.add_argument("limit", type=int, nargs="?", default=100)
    parser.add_argument(
        "--segment-candidates",
        type=int,
        default=DEFAULT_SEGMENT_CANDIDATES,
        help="Wheel-30 candidates represented per segment",
    )
    parser.add_argument(
        "--packed",
        action="store_true",
        help="Store the complete output in a packed array instead of streaming",
    )
    args = parser.parse_args()

    if args.packed:
        primes, stats = primes_wheel30_segmented_packed(
            args.limit, args.segment_candidates
        )
        print(list(primes))
        print(stats)
    else:
        count, checksum, last, stats = consume_primes_wheel30_segmented(
            args.limit, args.segment_candidates
        )
        print(f"count={count}, checksum={checksum}, last={last}")
        print(stats)

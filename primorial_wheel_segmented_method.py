"""غربلة مجزأة عامة بعجلة primorial.

تسمح باختبار عجلات مثل 30، 210، 2310، 30030 مع نفس منطق التمثيل:
- لا نمثل إلا البواقي المتباينة أوليًا مع M.
- بت واحد لكل مرشح ممثل.
- مقاطع ثابتة في عدد المرشحين.
- إخراج streaming.
"""
from __future__ import annotations

from array import array
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from math import gcd, isqrt, prod
from typing import Iterator, Sequence, Tuple

from segmented_method import _base_primes

DEFAULT_SEGMENT_CANDIDATES = 32_768


@dataclass(frozen=True)
class WheelSpec:
    primes: tuple[int, ...]
    modulus: int
    residues: tuple[int, ...]
    residue_rank: tuple[int, ...]


@dataclass
class Stats:
    limit: int
    wheel_modulus: int
    wheel_phi: int
    wheel_density: float
    wheel_primes: tuple[int, ...]
    segment_candidate_capacity: int
    segments_processed: int = 0
    represented_candidates: int = 0
    max_segment_storage_bytes: int = 0
    strike_attempts: int = 0
    newly_removed: int = 0
    yielded_primes: int = 0
    base_prime_count: int = 0
    base_storage_bytes: int = 0


def make_wheel(wheel_primes: Sequence[int]) -> WheelSpec:
    ps = tuple(int(p) for p in wheel_primes)
    if not ps or ps[0] != 2:
        raise ValueError("wheel_primes must start with 2")
    if tuple(sorted(set(ps))) != ps:
        raise ValueError("wheel_primes must be sorted distinct primes")
    modulus = prod(ps)
    residues = tuple(n for n in range(1, modulus) if gcd(n, modulus) == 1)
    rank = [-1] * modulus
    for i, r in enumerate(residues):
        rank[r] = i
    return WheelSpec(ps, modulus, residues, tuple(rank))


def _is_alive(bits: bytearray, index: int) -> bool:
    return bool(bits[index >> 3] & (1 << (index & 7)))


def _clear(bits: bytearray, index: int) -> bool:
    byte_index = index >> 3
    mask = 1 << (index & 7)
    if bits[byte_index] & mask:
        bits[byte_index] &= ~mask
        return True
    return False


def _candidate_value(spec: WheelSpec, index: int) -> int:
    block, rank = divmod(index, len(spec.residues))
    return block * spec.modulus + spec.residues[rank]


def _candidate_index(spec: WheelSpec, value: int) -> int:
    block, residue = divmod(value, spec.modulus)
    rank = spec.residue_rank[residue]
    if rank < 0:
        raise ValueError(f"{value} is not a candidate for wheel {spec.modulus}")
    return block * len(spec.residues) + rank


def _first_candidate_index_at_least(spec: WheelSpec, value: int) -> int:
    block, residue = divmod(value, spec.modulus)
    rank = bisect_left(spec.residues, residue)
    if rank < len(spec.residues):
        return block * len(spec.residues) + rank
    return (block + 1) * len(spec.residues)


def _last_candidate_index_at_most(spec: WheelSpec, value: int) -> int:
    block, residue = divmod(value, spec.modulus)
    rank = bisect_right(spec.residues, residue) - 1
    if rank >= 0:
        return block * len(spec.residues) + rank
    return (block - 1) * len(spec.residues) + len(spec.residues) - 1


def _iter_primes(
    limit: int,
    spec: WheelSpec,
    segment_candidate_count: int,
    stats: Stats,
) -> Iterator[int]:
    for p in spec.primes:
        if p <= limit:
            stats.yielded_primes += 1
            yield p

    first_after_wheel = spec.primes[-1] + 1
    if limit < first_after_wheel:
        return

    root = isqrt(limit)
    base_primes, base_candidate_bytes = _base_primes(root)
    stats.base_prime_count = len(base_primes)
    stats.base_storage_bytes = base_candidate_bytes + len(base_primes) * base_primes.itemsize

    start_index = _first_candidate_index_at_least(spec, first_after_wheel)
    last_index = _last_candidate_index_at_most(spec, limit)

    while start_index <= last_index:
        end_index = min(last_index, start_index + segment_candidate_count - 1)
        count = end_index - start_index + 1
        low = _candidate_value(spec, start_index)
        high = _candidate_value(spec, end_index)
        byte_count = (count + 7) // 8
        alive = bytearray(b"\xff") * byte_count

        stats.segments_processed += 1
        stats.represented_candidates += count
        stats.max_segment_storage_bytes = max(stats.max_segment_storage_bytes, byte_count)

        for p in base_primes:
            if p <= spec.primes[-1]:
                continue
            if p * p > high:
                break

            k_min = max(p, (low + p - 1) // p)
            k_index = _first_candidate_index_at_least(spec, k_min)
            while True:
                k = _candidate_value(spec, k_index)
                multiple = p * k
                if multiple > high:
                    break
                stats.strike_attempts += 1
                local = _candidate_index(spec, multiple) - start_index
                if _clear(alive, local):
                    stats.newly_removed += 1
                k_index += 1

        for local in range(count):
            if _is_alive(alive, local):
                n = _candidate_value(spec, start_index + local)
                if n <= limit:
                    stats.yielded_primes += 1
                    yield n
        start_index = end_index + 1


def consume_primes_primorial_wheel(
    limit: int,
    wheel_primes: Sequence[int],
    segment_candidate_count: int = DEFAULT_SEGMENT_CANDIDATES,
) -> Tuple[int, int, int, Stats]:
    if limit < 0:
        raise ValueError("limit must be >= 0")
    if segment_candidate_count < 1:
        raise ValueError("segment_candidate_count must be >= 1")
    spec = make_wheel(wheel_primes)
    stats = Stats(
        limit=limit,
        wheel_modulus=spec.modulus,
        wheel_phi=len(spec.residues),
        wheel_density=len(spec.residues) / spec.modulus,
        wheel_primes=spec.primes,
        segment_candidate_capacity=segment_candidate_count,
    )
    count = checksum = last = 0
    for p in _iter_primes(limit, spec, segment_candidate_count, stats):
        count += 1
        checksum += p
        last = p
    return count, checksum, last, stats


def iter_primes_primorial_wheel(
    limit: int,
    wheel_primes: Sequence[int],
    segment_candidate_count: int = DEFAULT_SEGMENT_CANDIDATES,
) -> Iterator[int]:
    if limit < 0:
        raise ValueError("limit must be >= 0")
    spec = make_wheel(wheel_primes)
    stats = Stats(
        limit=limit,
        wheel_modulus=spec.modulus,
        wheel_phi=len(spec.residues),
        wheel_density=len(spec.residues) / spec.modulus,
        wheel_primes=spec.primes,
        segment_candidate_capacity=segment_candidate_count,
    )
    yield from _iter_primes(limit, spec, segment_candidate_count, stats)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("limit", type=int, nargs="?", default=1_000_000)
    parser.add_argument("--wheel-primes", default="2,3,5,7,11")
    parser.add_argument("--segment-candidates", type=int, default=DEFAULT_SEGMENT_CANDIDATES)
    args = parser.parse_args()
    wheel_primes = tuple(int(x) for x in args.wheel_primes.split(",") if x)
    count, checksum, last, stats = consume_primes_primorial_wheel(
        args.limit, wheel_primes, args.segment_candidates
    )
    print(f"wheel={stats.wheel_modulus}, phi={stats.wheel_phi}, density={stats.wheel_density:.9f}")
    print(f"count={count:,}, last={last:,}, checksum={checksum:,}")
    print(stats)

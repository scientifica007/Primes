"""Cache-aware segmented prime sieve using compact primorial-wheel tables.

This is an experimental companion to ``primorial_wheel_segmented_method.py``.
The mathematical sieve is the same, but the wheel representation is packed:

* residues use unsigned 8/16/32-bit arrays when possible;
* residue -> rank uses signed 8/16/32-bit arrays with -1 as the sentinel;
* WheelSpec objects can be reused across many segment-size experiments;
* the hot loops keep wheel fields in local variables and advance candidate
  block/rank state incrementally instead of repeatedly decoding every index.

The goal is to separate the value of a larger wheel from accidental Python
object overhead and to make cache-size experiments meaningful.
"""
from __future__ import annotations

from array import array
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from functools import lru_cache
from math import gcd, isqrt, prod
from typing import Iterator, Sequence, Tuple

from segmented_method import _base_primes

DEFAULT_SEGMENT_CANDIDATES = 32_768


@dataclass(frozen=True)
class CompactWheelSpec:
    primes: tuple[int, ...]
    modulus: int
    phi: int
    residues: array
    residue_rank: array
    residues_typecode: str
    rank_typecode: str
    residues_payload_bytes: int
    rank_payload_bytes: int

    @property
    def table_payload_bytes(self) -> int:
        return self.residues_payload_bytes + self.rank_payload_bytes

    @property
    def density(self) -> float:
        return self.phi / self.modulus


@dataclass
class Stats:
    limit: int
    wheel_modulus: int
    wheel_phi: int
    wheel_density: float
    wheel_primes: tuple[int, ...]
    segment_candidate_capacity: int
    wheel_residue_bytes: int = 0
    wheel_rank_bytes: int = 0
    wheel_table_payload_bytes: int = 0
    segments_processed: int = 0
    represented_candidates: int = 0
    max_segment_storage_bytes: int = 0
    strike_attempts: int = 0
    newly_removed: int = 0
    yielded_primes: int = 0
    base_prime_count: int = 0
    base_storage_bytes: int = 0

    @property
    def algorithmic_working_set_bytes(self) -> int:
        """Main packed tables + largest sieve segment + base-prime storage.

        This deliberately excludes Python object headers and interpreter state.
        It is a stable algorithmic payload useful for cross-wheel/cache studies.
        """
        return (
            self.wheel_table_payload_bytes
            + self.max_segment_storage_bytes
            + self.base_storage_bytes
        )


def _unsigned_typecode(max_value: int) -> str:
    if max_value <= 0xFF:
        return "B"
    if max_value <= 0xFFFF:
        return "H"
    if max_value <= 0xFFFFFFFF:
        return "I"
    return "Q"


def _signed_typecode(max_rank: int) -> str:
    if max_rank <= 0x7F:
        return "b"
    if max_rank <= 0x7FFF:
        return "h"
    if max_rank <= 0x7FFFFFFF:
        return "i"
    return "q"


def _validate_wheel_primes(wheel_primes: Sequence[int]) -> tuple[int, ...]:
    ps = tuple(int(p) for p in wheel_primes)
    if not ps or ps[0] != 2:
        raise ValueError("wheel_primes must start with 2")
    if tuple(sorted(set(ps))) != ps:
        raise ValueError("wheel_primes must be sorted distinct primes")
    return ps


@lru_cache(maxsize=32)
def _make_compact_wheel_cached(ps: tuple[int, ...]) -> CompactWheelSpec:
    modulus = prod(ps)
    residue_values = [n for n in range(1, modulus) if gcd(n, modulus) == 1]
    phi = len(residue_values)

    residue_code = _unsigned_typecode(modulus - 1)
    rank_code = _signed_typecode(phi - 1)
    residues = array(residue_code, residue_values)
    rank = array(rank_code, [-1]) * modulus
    for i, r in enumerate(residue_values):
        rank[r] = i

    return CompactWheelSpec(
        primes=ps,
        modulus=modulus,
        phi=phi,
        residues=residues,
        residue_rank=rank,
        residues_typecode=residue_code,
        rank_typecode=rank_code,
        residues_payload_bytes=len(residues) * residues.itemsize,
        rank_payload_bytes=len(rank) * rank.itemsize,
    )


def make_compact_wheel(wheel_primes: Sequence[int]) -> CompactWheelSpec:
    """Build or reuse a compact wheel specification."""
    return _make_compact_wheel_cached(_validate_wheel_primes(wheel_primes))


def _is_alive(bits: bytearray, index: int) -> bool:
    return bool(bits[index >> 3] & (1 << (index & 7)))


def _clear(bits: bytearray, index: int) -> bool:
    byte_index = index >> 3
    mask = 1 << (index & 7)
    if bits[byte_index] & mask:
        bits[byte_index] &= ~mask
        return True
    return False


def _first_candidate_index_at_least(spec: CompactWheelSpec, value: int) -> int:
    block, residue = divmod(value, spec.modulus)
    r = bisect_left(spec.residues, residue)
    if r < spec.phi:
        return block * spec.phi + r
    return (block + 1) * spec.phi


def _last_candidate_index_at_most(spec: CompactWheelSpec, value: int) -> int:
    block, residue = divmod(value, spec.modulus)
    r = bisect_right(spec.residues, residue) - 1
    if r >= 0:
        return block * spec.phi + r
    return (block - 1) * spec.phi + spec.phi - 1


def _iter_primes_with_spec(
    limit: int,
    spec: CompactWheelSpec,
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

    modulus = spec.modulus
    phi = spec.phi
    residues = spec.residues
    rank_map = spec.residue_rank
    last_wheel_prime = spec.primes[-1]

    start_index = _first_candidate_index_at_least(spec, first_after_wheel)
    last_index = _last_candidate_index_at_most(spec, limit)

    while start_index <= last_index:
        end_index = min(last_index, start_index + segment_candidate_count - 1)
        count = end_index - start_index + 1

        low_block, low_rank = divmod(start_index, phi)
        high_block, high_rank = divmod(end_index, phi)
        low = low_block * modulus + residues[low_rank]
        high = high_block * modulus + residues[high_rank]

        byte_count = (count + 7) // 8
        alive = bytearray(b"\xff") * byte_count
        stats.segments_processed += 1
        stats.represented_candidates += count
        stats.max_segment_storage_bytes = max(stats.max_segment_storage_bytes, byte_count)

        for p in base_primes:
            if p <= last_wheel_prime:
                continue
            if p * p > high:
                break

            k_min = max(p, (low + p - 1) // p)
            k_index = _first_candidate_index_at_least(spec, k_min)
            k_block, k_rank = divmod(k_index, phi)

            while True:
                k = k_block * modulus + residues[k_rank]
                multiple = p * k
                if multiple > high:
                    break

                stats.strike_attempts += 1
                m_block, m_residue = divmod(multiple, modulus)
                m_rank = rank_map[m_residue]
                # p and k are both coprime to the wheel modulus, so this rank
                # must exist. Keeping the guard makes corruption obvious.
                if m_rank < 0:
                    raise AssertionError("wheel candidate multiplication lost coprimality")
                local = m_block * phi + m_rank - start_index
                if _clear(alive, local):
                    stats.newly_removed += 1

                k_rank += 1
                if k_rank == phi:
                    k_rank = 0
                    k_block += 1

        out_block, out_rank = divmod(start_index, phi)
        for local in range(count):
            if _is_alive(alive, local):
                n = out_block * modulus + residues[out_rank]
                if n <= limit:
                    stats.yielded_primes += 1
                    yield n
            out_rank += 1
            if out_rank == phi:
                out_rank = 0
                out_block += 1

        start_index = end_index + 1


def _new_stats(limit: int, spec: CompactWheelSpec, segment_candidate_count: int) -> Stats:
    return Stats(
        limit=limit,
        wheel_modulus=spec.modulus,
        wheel_phi=spec.phi,
        wheel_density=spec.density,
        wheel_primes=spec.primes,
        segment_candidate_capacity=segment_candidate_count,
        wheel_residue_bytes=spec.residues_payload_bytes,
        wheel_rank_bytes=spec.rank_payload_bytes,
        wheel_table_payload_bytes=spec.table_payload_bytes,
    )


def consume_primes_compact_spec(
    limit: int,
    spec: CompactWheelSpec,
    segment_candidate_count: int = DEFAULT_SEGMENT_CANDIDATES,
) -> Tuple[int, int, int, Stats]:
    if limit < 0:
        raise ValueError("limit must be >= 0")
    if segment_candidate_count < 1:
        raise ValueError("segment_candidate_count must be >= 1")
    stats = _new_stats(limit, spec, segment_candidate_count)
    count = checksum = last = 0
    for p in _iter_primes_with_spec(limit, spec, segment_candidate_count, stats):
        count += 1
        checksum += p
        last = p
    return count, checksum, last, stats


def consume_primes_compact_wheel(
    limit: int,
    wheel_primes: Sequence[int],
    segment_candidate_count: int = DEFAULT_SEGMENT_CANDIDATES,
) -> Tuple[int, int, int, Stats]:
    return consume_primes_compact_spec(
        limit,
        make_compact_wheel(wheel_primes),
        segment_candidate_count,
    )


def iter_primes_compact_wheel(
    limit: int,
    wheel_primes: Sequence[int],
    segment_candidate_count: int = DEFAULT_SEGMENT_CANDIDATES,
) -> Iterator[int]:
    if limit < 0:
        raise ValueError("limit must be >= 0")
    if segment_candidate_count < 1:
        raise ValueError("segment_candidate_count must be >= 1")
    spec = make_compact_wheel(wheel_primes)
    stats = _new_stats(limit, spec, segment_candidate_count)
    yield from _iter_primes_with_spec(limit, spec, segment_candidate_count, stats)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Compact cache-aware primorial wheel sieve")
    parser.add_argument("limit", type=int, nargs="?", default=1_000_000)
    parser.add_argument("--wheel-primes", default="2,3,5,7")
    parser.add_argument("--segment-candidates", type=int, default=DEFAULT_SEGMENT_CANDIDATES)
    args = parser.parse_args()
    wheel_primes = tuple(int(x) for x in args.wheel_primes.split(",") if x)
    spec = make_compact_wheel(wheel_primes)
    count, checksum, last, stats = consume_primes_compact_spec(
        args.limit, spec, args.segment_candidates
    )
    print(
        f"wheel={spec.modulus}, phi={spec.phi}, density={spec.density:.9f}, "
        f"table_payload={spec.table_payload_bytes:,} B"
    )
    print(f"count={count:,}, last={last:,}, checksum={checksum:,}")
    print(stats)

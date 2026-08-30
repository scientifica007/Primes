"""اختبارات صحة للنسخ وتمثيلات الناتج والتتبع التعليمي."""

import unittest

from optimized_method import primes_optimized
from original_method import primes_original
from retain_prime_compact_method import primes_retain_compact
from retain_prime_method import primes_retain
from retain_prime_packed_output_method import (
    iter_primes_retain_packed,
    primes_retain_packed,
)
from segmented_method import (
    consume_primes_segmented,
    iter_primes_segmented,
    primes_segmented_packed,
)
from trace_method import trace_primes
from wheel30_segmented_method import (
    consume_primes_wheel30_segmented,
    iter_primes_wheel30_segmented,
    primes_wheel30_segmented_packed,
)


EXPECTED_100 = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97,
]


class PrimeSearchTests(unittest.TestCase):
    def test_known_primes_to_100(self):
        original, _ = primes_original(100)
        retained, _ = primes_retain(100)
        compact, _ = primes_retain_compact(100)
        packed, _ = primes_retain_packed(100)
        streamed = list(iter_primes_retain_packed(100))
        segmented = list(iter_primes_segmented(100, segment_odd_count=8))
        segmented_packed, _ = primes_segmented_packed(100, segment_odd_count=8)
        wheel30 = list(iter_primes_wheel30_segmented(100, segment_candidate_count=5))
        wheel30_packed, _ = primes_wheel30_segmented_packed(
            100, segment_candidate_count=5
        )
        optimized, _ = primes_optimized(100)
        traced, stages, _, _, _ = trace_primes(100)

        self.assertEqual(original, EXPECTED_100)
        self.assertEqual(retained, EXPECTED_100)
        self.assertEqual(compact, EXPECTED_100)
        self.assertEqual(list(packed), EXPECTED_100)
        self.assertEqual(streamed, EXPECTED_100)
        self.assertEqual(segmented, EXPECTED_100)
        self.assertEqual(list(segmented_packed), EXPECTED_100)
        self.assertEqual(wheel30, EXPECTED_100)
        self.assertEqual(list(wheel30_packed), EXPECTED_100)
        self.assertEqual(optimized, EXPECTED_100)
        self.assertEqual(traced, EXPECTED_100)
        self.assertEqual([stage.label for stage in stages], ["ب", "ج", "د", "هـ"])
        self.assertEqual(stages[0].new_primes, [5, 7])
        self.assertEqual(stages[1].new_primes, [11, 13, 17, 19, 23])
        self.assertEqual(stages[2].new_primes, [29, 31, 37, 41, 43, 47])
        self.assertEqual(stages[3].new_primes, [53, 59, 61, 67, 71, 73, 79, 83, 89, 97])

    def test_small_limits(self):
        expected = {
            0: [],
            1: [],
            2: [2],
            3: [2, 3],
            4: [2, 3],
            5: [2, 3, 5],
            8: [2, 3, 5, 7],
            9: [2, 3, 5, 7],
            10: [2, 3, 5, 7],
        }
        for limit, primes in expected.items():
            with self.subTest(limit=limit):
                self.assertEqual(primes_original(limit)[0], primes)
                self.assertEqual(primes_retain(limit)[0], primes)
                self.assertEqual(primes_retain_compact(limit)[0], primes)
                self.assertEqual(list(primes_retain_packed(limit)[0]), primes)
                self.assertEqual(list(iter_primes_retain_packed(limit)), primes)
                self.assertEqual(list(iter_primes_segmented(limit, 3)), primes)
                self.assertEqual(list(primes_segmented_packed(limit, 3)[0]), primes)
                self.assertEqual(list(iter_primes_wheel30_segmented(limit, 3)), primes)
                self.assertEqual(
                    list(primes_wheel30_segmented_packed(limit, 3)[0]), primes
                )
                self.assertEqual(primes_optimized(limit)[0], primes)
                self.assertEqual(trace_primes(limit)[0], primes)

    def test_methods_agree_over_range(self):
        for limit in range(0, 501):
            with self.subTest(limit=limit):
                original = primes_original(limit)[0]
                self.assertEqual(original, primes_retain(limit)[0])
                self.assertEqual(original, primes_retain_compact(limit)[0])
                self.assertEqual(original, list(primes_retain_packed(limit)[0]))
                self.assertEqual(original, list(iter_primes_retain_packed(limit)))
                self.assertEqual(original, list(iter_primes_segmented(limit, 17)))
                self.assertEqual(
                    original, list(iter_primes_wheel30_segmented(limit, 11))
                )
                self.assertEqual(original, primes_optimized(limit)[0])
                self.assertEqual(original, trace_primes(limit)[0])

    def test_retain_method_removes_only_composites(self):
        primes, stats = primes_retain(100)
        self.assertEqual(primes, EXPECTED_100)
        self.assertEqual(stats.newly_removed, 74)
        self.assertEqual(stats.processed_primes, 4)
        for composite in [4, 9, 25, 49]:
            self.assertNotIn(composite, primes)
        for prime in [2, 3, 5, 7]:
            self.assertIn(prime, primes)

    def test_compact_method_uses_bit_storage(self):
        primes, stats = primes_retain_compact(100)
        self.assertEqual(primes, EXPECTED_100)
        self.assertEqual(stats.represented_odd_candidates, 49)
        self.assertEqual(stats.storage_bytes, 7)
        self.assertEqual(stats.processed_primes, 3)
        self.assertEqual(stats.newly_removed, 25)

    def test_packed_output_compresses_result(self):
        primes, stats = primes_retain_packed(100)
        self.assertEqual(list(primes), EXPECTED_100)
        self.assertEqual(stats.represented_odd_candidates, 49)
        self.assertEqual(stats.candidate_storage_bytes, 7)
        self.assertIn(stats.output_typecode, ("I", "Q"))
        self.assertEqual(stats.output_itemsize, primes.itemsize)
        self.assertEqual(stats.output_storage_bytes, len(primes) * primes.itemsize)
        self.assertLess(stats.output_storage_bytes, 25 * 28)

    def test_streaming_output_matches_without_result_container(self):
        self.assertEqual(list(iter_primes_retain_packed(100)), EXPECTED_100)

    def test_segmented_stream_uses_bounded_segments(self):
        count, checksum, last, stats = consume_primes_segmented(
            100, segment_odd_count=8
        )
        self.assertEqual(count, len(EXPECTED_100))
        self.assertEqual(checksum, sum(EXPECTED_100))
        self.assertEqual(last, 97)
        self.assertEqual(stats.base_limit, 10)
        self.assertEqual(stats.base_prime_count, 4)
        self.assertEqual(stats.segments_processed, 7)
        self.assertLessEqual(stats.max_segment_odd_count, 8)
        self.assertLessEqual(stats.max_segment_storage_bytes, 1)
        self.assertEqual(stats.newly_removed, 25)
        self.assertEqual(stats.yielded_primes, 25)

    def test_wheel30_represents_only_coprime_candidates(self):
        count, checksum, last, stats = consume_primes_wheel30_segmented(
            100, segment_candidate_count=8
        )
        self.assertEqual(count, len(EXPECTED_100))
        self.assertEqual(checksum, sum(EXPECTED_100))
        self.assertEqual(last, 97)
        self.assertEqual(stats.base_limit, 10)
        self.assertEqual(stats.base_prime_count, 4)
        self.assertEqual(stats.represented_candidates, 25)
        self.assertEqual(stats.segments_processed, 4)
        self.assertLessEqual(stats.max_segment_candidate_count, 8)
        self.assertLessEqual(stats.max_segment_storage_bytes, 1)
        self.assertEqual(stats.strike_attempts, 3)
        self.assertEqual(stats.newly_removed, 3)
        self.assertEqual(stats.yielded_primes, 25)

    def test_trace_preserves_one_but_never_classifies_it(self):
        primes, stages, _, _, list_a = trace_primes(100)
        self.assertIn(1, list_a)
        self.assertNotIn(1, primes)
        self.assertTrue(all(1 in stage.remaining for stage in stages))
        self.assertTrue(all(1 not in stage.new_primes for stage in stages))

    def test_negative_limit_rejected(self):
        with self.assertRaises(ValueError):
            primes_original(-1)
        with self.assertRaises(ValueError):
            primes_retain(-1)
        with self.assertRaises(ValueError):
            primes_retain_compact(-1)
        with self.assertRaises(ValueError):
            primes_retain_packed(-1)
        with self.assertRaises(ValueError):
            list(iter_primes_retain_packed(-1))
        with self.assertRaises(ValueError):
            list(iter_primes_segmented(-1))
        with self.assertRaises(ValueError):
            consume_primes_segmented(-1)
        with self.assertRaises(ValueError):
            primes_segmented_packed(-1)
        with self.assertRaises(ValueError):
            list(iter_primes_wheel30_segmented(-1))
        with self.assertRaises(ValueError):
            consume_primes_wheel30_segmented(-1)
        with self.assertRaises(ValueError):
            primes_wheel30_segmented_packed(-1)
        with self.assertRaises(ValueError):
            primes_optimized(-1)
        with self.assertRaises(ValueError):
            trace_primes(-1)

    def test_invalid_segment_size_rejected(self):
        with self.assertRaises(ValueError):
            list(iter_primes_segmented(100, 0))
        with self.assertRaises(ValueError):
            consume_primes_segmented(100, 0)
        with self.assertRaises(ValueError):
            primes_segmented_packed(100, 0)
        with self.assertRaises(ValueError):
            list(iter_primes_wheel30_segmented(100, 0))
        with self.assertRaises(ValueError):
            consume_primes_wheel30_segmented(100, 0)
        with self.assertRaises(ValueError):
            primes_wheel30_segmented_packed(100, 0)


if __name__ == "__main__":
    unittest.main()

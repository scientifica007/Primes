"""اختبارات صحة للنسخ والتتبع التعليمي."""

import unittest

from optimized_method import primes_optimized
from original_method import primes_original
from retain_prime_compact_method import primes_retain_compact
from retain_prime_method import primes_retain
from trace_method import trace_primes


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
        optimized, _ = primes_optimized(100)
        traced, stages, _, _, _ = trace_primes(100)
        self.assertEqual(original, EXPECTED_100)
        self.assertEqual(retained, EXPECTED_100)
        self.assertEqual(compact, EXPECTED_100)
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
                self.assertEqual(primes_optimized(limit)[0], primes)
                self.assertEqual(trace_primes(limit)[0], primes)

    def test_methods_agree_over_range(self):
        for limit in range(0, 501):
            with self.subTest(limit=limit):
                original = primes_original(limit)[0]
                self.assertEqual(original, primes_retain(limit)[0])
                self.assertEqual(original, primes_retain_compact(limit)[0])
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
        for composite in [9, 25, 49, 77, 91]:
            self.assertNotIn(composite, primes)
        for prime in [2, 3, 5, 7]:
            self.assertIn(prime, primes)

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
            primes_optimized(-1)
        with self.assertRaises(ValueError):
            trace_primes(-1)


if __name__ == "__main__":
    unittest.main()

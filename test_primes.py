"""اختبارات صحة للنسختين."""

import unittest

from optimized_method import primes_optimized
from original_method import primes_original


EXPECTED_100 = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97,
]


class PrimeSearchTests(unittest.TestCase):
    def test_known_primes_to_100(self):
        original, _ = primes_original(100)
        optimized, _ = primes_optimized(100)
        self.assertEqual(original, EXPECTED_100)
        self.assertEqual(optimized, EXPECTED_100)

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
                self.assertEqual(primes_optimized(limit)[0], primes)

    def test_methods_agree_over_range(self):
        for limit in range(0, 501):
            with self.subTest(limit=limit):
                self.assertEqual(
                    primes_original(limit)[0],
                    primes_optimized(limit)[0],
                )

    def test_negative_limit_rejected(self):
        with self.assertRaises(ValueError):
            primes_original(-1)
        with self.assertRaises(ValueError):
            primes_optimized(-1)


if __name__ == "__main__":
    unittest.main()

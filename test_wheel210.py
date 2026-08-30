"""اختبارات مستقلة لـ Wheel-210 ومقارنتها بـ Wheel-30."""

import unittest

from original_method import primes_original
from wheel30_segmented_method import iter_primes_wheel30_segmented
from wheel210_segmented_method import (
    consume_primes_wheel210_segmented,
    iter_primes_wheel210_segmented,
    primes_wheel210_segmented_packed,
)


EXPECTED_100 = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97,
]


class Wheel210Tests(unittest.TestCase):
    def test_known_primes_to_100(self):
        self.assertEqual(
            list(iter_primes_wheel210_segmented(100, 5)), EXPECTED_100
        )
        packed, _ = primes_wheel210_segmented_packed(100, 5)
        self.assertEqual(list(packed), EXPECTED_100)

    def test_agrees_over_range(self):
        for limit in range(0, 501):
            with self.subTest(limit=limit):
                reference = primes_original(limit)[0]
                self.assertEqual(
                    list(iter_primes_wheel210_segmented(limit, 7)), reference
                )
                self.assertEqual(
                    list(iter_primes_wheel30_segmented(limit, 7)), reference
                )

    def test_wheel210_removes_seven_from_representation(self):
        count, checksum, last, stats = consume_primes_wheel210_segmented(
            100, segment_candidate_count=8
        )
        self.assertEqual(count, 25)
        self.assertEqual(checksum, sum(EXPECTED_100))
        self.assertEqual(last, 97)
        # الأعداد الممثلة من 11 إلى 100 وغير القابلة للقسمة على 2,3,5,7.
        self.assertEqual(stats.represented_candidates, 21)
        self.assertEqual(stats.newly_removed, 0)
        self.assertEqual(stats.strike_attempts, 0)
        self.assertEqual(stats.yielded_primes, 25)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            list(iter_primes_wheel210_segmented(-1))
        with self.assertRaises(ValueError):
            list(iter_primes_wheel210_segmented(100, 0))
        with self.assertRaises(ValueError):
            consume_primes_wheel210_segmented(100, 0)
        with self.assertRaises(ValueError):
            primes_wheel210_segmented_packed(100, 0)


if __name__ == "__main__":
    unittest.main()

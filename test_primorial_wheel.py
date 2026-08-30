import unittest

from original_method import primes_original
from primorial_wheel_segmented_method import (
    consume_primes_primorial_wheel,
    iter_primes_primorial_wheel,
    make_wheel,
)


class PrimorialWheelTests(unittest.TestCase):
    def test_wheel_parameters(self):
        cases = {
            (2, 3): (6, 2),
            (2, 3, 5): (30, 8),
            (2, 3, 5, 7): (210, 48),
            (2, 3, 5, 7, 11): (2310, 480),
            (2, 3, 5, 7, 11, 13): (30030, 5760),
        }
        for ps, expected in cases.items():
            with self.subTest(ps=ps):
                spec = make_wheel(ps)
                self.assertEqual((spec.modulus, len(spec.residues)), expected)

    def test_wheels_match_reference_to_1000(self):
        expected = primes_original(1000)[0]
        wheels = [
            (2, 3),
            (2, 3, 5),
            (2, 3, 5, 7),
            (2, 3, 5, 7, 11),
            (2, 3, 5, 7, 11, 13),
        ]
        for ps in wheels:
            with self.subTest(ps=ps):
                self.assertEqual(list(iter_primes_primorial_wheel(1000, ps, 37)), expected)

    def test_stats(self):
        count, checksum, last, stats = consume_primes_primorial_wheel(
            1000, (2, 3, 5, 7, 11), 64
        )
        expected = primes_original(1000)[0]
        self.assertEqual(count, len(expected))
        self.assertEqual(checksum, sum(expected))
        self.assertEqual(last, expected[-1])
        self.assertEqual(stats.wheel_modulus, 2310)
        self.assertEqual(stats.wheel_phi, 480)
        self.assertLessEqual(stats.max_segment_storage_bytes, 8)


if __name__ == "__main__":
    unittest.main()

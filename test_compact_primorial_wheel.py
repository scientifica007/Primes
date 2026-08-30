import unittest

from original_method import primes_original
from compact_primorial_wheel_method import (
    consume_primes_compact_spec,
    iter_primes_compact_wheel,
    make_compact_wheel,
)


WHEELS = [
    (2, 3, 5),
    (2, 3, 5, 7),
    (2, 3, 5, 7, 11),
    (2, 3, 5, 7, 11, 13),
]


class CompactPrimorialWheelTests(unittest.TestCase):
    def test_compact_table_parameters(self):
        expected = {
            (2, 3, 5): (30, 8),
            (2, 3, 5, 7): (210, 48),
            (2, 3, 5, 7, 11): (2310, 480),
            (2, 3, 5, 7, 11, 13): (30030, 5760),
        }
        for ps, values in expected.items():
            with self.subTest(ps=ps):
                spec = make_compact_wheel(ps)
                self.assertEqual((spec.modulus, spec.phi), values)
                self.assertEqual(spec.table_payload_bytes,
                                 spec.residues_payload_bytes + spec.rank_payload_bytes)

    def test_30030_table_is_compact(self):
        spec = make_compact_wheel((2, 3, 5, 7, 11, 13))
        self.assertEqual(spec.residues_typecode, "H")
        self.assertEqual(spec.rank_typecode, "h")
        self.assertEqual(spec.residues_payload_bytes, 5760 * 2)
        self.assertEqual(spec.rank_payload_bytes, 30030 * 2)
        self.assertEqual(spec.table_payload_bytes, 71580)

    def test_wheels_match_reference_across_small_segments(self):
        expected = primes_original(5000)[0]
        for ps in WHEELS:
            for segment_candidates in (1, 7, 37, 256, 4096):
                with self.subTest(ps=ps, segment=segment_candidates):
                    actual = list(iter_primes_compact_wheel(5000, ps, segment_candidates))
                    self.assertEqual(actual, expected)

    def test_reused_spec_and_stats(self):
        spec = make_compact_wheel((2, 3, 5, 7))
        count, checksum, last, stats = consume_primes_compact_spec(10000, spec, 32768)
        expected = primes_original(10000)[0]
        self.assertEqual(count, len(expected))
        self.assertEqual(checksum, sum(expected))
        self.assertEqual(last, expected[-1])
        self.assertEqual(stats.wheel_modulus, 210)
        self.assertEqual(stats.wheel_table_payload_bytes, spec.table_payload_bytes)
        self.assertLessEqual(stats.max_segment_storage_bytes, 4096)
        self.assertGreater(stats.algorithmic_working_set_bytes, 0)

    def test_edge_cases(self):
        for limit in range(0, 20):
            expected = primes_original(limit)[0]
            for ps in WHEELS:
                with self.subTest(limit=limit, ps=ps):
                    self.assertEqual(list(iter_primes_compact_wheel(limit, ps, 8)), expected)


if __name__ == "__main__":
    unittest.main()

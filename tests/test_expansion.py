# -*- coding: utf-8 -*-
"""动态膨胀模块单元测试（0.14/0.15 宇宙八相变）. """
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest

from cosmogeo import expansion
from cosmogeo.constants import EPS0


class TestExpansion(unittest.TestCase):
    def test_epsilon(self):
        self.assertAlmostEqual(expansion.epsilon(0), EPS0, delta=1e-9)
        self.assertLess(expansion.epsilon(5), expansion.epsilon(0), "ε 应随 s 衰减")

    def test_w_eff(self):
        w = expansion.w_eff(0)
        self.assertAlmostEqual(w, -1 - EPS0, delta=1e-9)
        # 几乎 de Sitter：N₅→N₇ 区间（s>0.4）内 |w+1| < 0.0044（0.14 推论 0.14.2.02）
        self.assertLess(abs(expansion.w_eff(1.0) + 1), 0.0044)

    def test_hubble_hz_zero(self):
        self.assertAlmostEqual(expansion.hubble_hz(0), 68.9, delta=0.5)

    def test_hubble_increasing(self):
        self.assertGreater(expansion.hubble_hz(3), expansion.hubble_hz(0))

    def test_distance_modulus(self):
        # 标准宇宙学对照：z=0.1 → μ≈38.2、z=1.0 → μ≈44.0
        self.assertAlmostEqual(expansion.distance_modulus(0.1), 38.2, delta=0.5)
        self.assertAlmostEqual(expansion.distance_modulus(1.0), 44.0, delta=0.5)

    def test_luminosity_distance(self):
        self.assertTrue(400 < expansion.luminosity_distance(0.1) < 500)

    def test_redshift_drift(self):
        rd = expansion.redshift_drift()
        self.assertLess(rd["magnitude_per_year"], 1e-10)

    def test_sn_deviation(self):
        self.assertEqual(expansion.sn_deviation(1.0), 0.0, "z≤1.5 无偏离")
        self.assertAlmostEqual(expansion.sn_deviation(2.0), 0.01, "z>1.5 更亮 0.01 mag")

    def test_phase_labels(self):
        self.assertEqual(expansion.phase_label(-5), "相区I 剧烈展开")
        self.assertIn("收缩", expansion.phase_label(0))

    def test_inflection_count(self):
        self.assertEqual(len(expansion.inflection_points()), 5)


if __name__ == "__main__":
    unittest.main()

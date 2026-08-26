# -*- coding: utf-8 -*-
"""8.10 太阳系 + 8.12 卫星系统单元测试."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest

from cosmogeo import solar, satellite
from cosmogeo.constants import S_E


class TestSolar(unittest.TestCase):
    def test_r0(self):
        self.assertAlmostEqual(solar.r0(), 0.3087, delta=1e-4)

    def test_rn_geometric(self):
        self.assertAlmostEqual(solar.r_n(4), 0.3087 * 16, delta=0.01)  # 木星壳层 4.94
        self.assertAlmostEqual(solar.r_n(7), 0.3087 * 128, delta=0.5)  # 39.5

    def test_planet_table_outer(self):
        table = {name: (r_geo, r_obs) for name, _, r_geo, r_obs, _ in solar.planet_table()}
        # 外层行星 5% 内（除海王星 31% 见第11章）
        self.assertAlmostEqual(table["土星"][0] / table["土星"][1], 1.031, delta=0.01)
        self.assertAlmostEqual(table["天王星"][0] / table["天王星"][1], 1.031, delta=0.01)
        # 冥王星 0.09%
        self.assertAlmostEqual(table["冥王星"][0] / table["冥王星"][1], 1.0009, delta=0.005)

    def test_orbital_fine_structure(self):
        fs0 = solar.orbital_fine_structure(0)
        self.assertAlmostEqual(fs0["S_orbital"], S_E, delta=0.1)
        self.assertAlmostEqual(fs0["alpha_orbital"], 1 / S_E, delta=1e-8)
        fs4 = solar.orbital_fine_structure(4)
        self.assertAlmostEqual(fs4["S_orbital"], S_E * 16, delta=0.1)

    def test_etno_attractor(self):
        self.assertAlmostEqual(solar.etno_attractor(), 247.28, delta=0.01)

    def test_spectral_gap(self):
        self.assertEqual(solar.spectral_gap_stability()["mu2"], 5.18)


class TestSatellite(unittest.TestCase):
    def test_period_ratio(self):
        theory = satellite.period_ratio_theory()
        self.assertAlmostEqual(theory, 5.0 / S_E, delta=1e-9)
        # 偏差 < 0.04%（文章自报，实测 ~0.013%）
        dev = satellite.period_ratio_deviation()
        self.assertLess(dev, 0.0004, f"周期比偏差应 <0.04%，实际 {dev*100:.3f}%")

    def test_mass_ratio(self):
        self.assertAlmostEqual(satellite.mass_ratio_theory(), 81.0, delta=0.1)
        self.assertLess(satellite.mass_ratio_deviation(), 0.07, "质量比偏差 <7%")

    def test_tidal_lock(self):
        self.assertAlmostEqual(satellite.tidal_lock_period_day(), 27.32, delta=0.01)

    def test_recession(self):
        self.assertAlmostEqual(satellite.recession_rate_cm_year(), 3.8, delta=0.1)

    def test_moon_distance(self):
        self.assertAlmostEqual(satellite.earth_moon_distance_km(), 384400.0, delta=100)


if __name__ == "__main__":
    unittest.main()

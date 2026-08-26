# -*- coding: utf-8 -*-
"""静态→动态转化测试（0.7.0）."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest

from cosmogeo import hubble, satellite, solar, cmb, df
from cosmogeo.constants import A0

M_MW = 1.2e41
R_SUN = 8.2 * 3.0857e19


class TestCosmicAge(unittest.TestCase):
    def test_age_universe(self):
        # 几何论无暗物质 → 年龄偏大（可检验差异）：14-20 Gyr
        age = hubble.age_universe_gyr()
        self.assertTrue(14 < age < 20, f"几何论宇宙年龄 {age:.1f} Gyr")

    def test_age_recombination(self):
        # t(z=1100) 应 ~0.1-1 Myr（量级）
        t_rec = hubble.cosmic_age(1100) * 1e3  # Gyr → Myr（×1e3）
        self.assertTrue(0.2 < t_rec < 2.0, f"复合时代 {t_rec:.2f} Myr（标准 0.38）")

    def test_age_monotonic(self):
        self.assertGreater(hubble.cosmic_age(0), hubble.cosmic_age(2))


class TestMoonEvolution(unittest.TestCase):
    def test_distance_past(self):
        d = satellite.moon_distance_past(1000)  # 10 亿年前
        self.assertAlmostEqual(d, 346400, delta=100)

    def test_distance_future(self):
        d = satellite.moon_distance_evolution(-1e8)  # 未来 1 亿年
        self.assertGreater(d, satellite.EARTH_MOON_DISTANCE_KM)

    def test_tidal_lock_epoch(self):
        t = satellite.tidal_lock_epoch()["epoch_years_ago"]
        self.assertTrue(5e9 < t < 2e10, f"潮汐锁定时代 {t/1e9:.1f} e9 年")


class TestMercuryPrecession(unittest.TestCase):
    def test_per_century(self):
        self.assertAlmostEqual(solar.mercury_precession_per_century(), 42.98, delta=0.01)

    def test_angle_evolution(self):
        self.assertAlmostEqual(solar.mercury_precession_angle(10), 429.8, delta=0.1)

    def test_orbits_per_century(self):
        self.assertAlmostEqual(solar.mercury_orbits_per_century(), 415.2, delta=1.0)


class TestEarlyUniverse(unittest.TestCase):
    def test_temperature_redshift(self):
        self.assertAlmostEqual(cmb.temperature_at_redshift(1100), 3000.2, delta=1.0)

    def test_temperature_time(self):
        t_early = cmb.temperature_at_time(1e6 * 3.15576e7)
        self.assertTrue(100 < t_early < 1000, f"T(1e6年)={t_early:.0f} K（辐射主导简化）")

    def test_recombination(self):
        ri = cmb.recombination_info()
        self.assertAlmostEqual(ri["z_star"], 1100)
        self.assertTrue(5 < ri["t_star_kyr"] < 50, f"t_*={ri['t_star_kyr']:.0f} kyr（简化模型）")

    def test_bbn(self):
        b = cmb.bbn_parameters()
        self.assertAlmostEqual(b["eta"], 6.1e-10)
        self.assertEqual(b["N_nu"], 3)


class TestOrbitSampling(unittest.TestCase):
    def test_sample_orbit(self):
        s = df.sample_orbit(R_SUN, M_MW)
        self.assertTrue(50e3 < s["v_sampled"] < 300e3, f"采样 v={s['v_sampled']/1e3:.0f} km/s")
        self.assertGreater(len(s["orbit"]["r"]), 50, "轨道应被积分")


if __name__ == "__main__":
    unittest.main()

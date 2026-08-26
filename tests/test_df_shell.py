# -*- coding: utf-8 -*-
"""8.11 壳层/DF 扩展单元测试（对照 galpy df 模块）. """
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
import math

from cosmogeo import galaxy, distribution, df
from cosmogeo.constants import A0

KPC = 3.0857e19
M_MW = 1.2e41
R_SUN = 8.2 * KPC


class TestShell(unittest.TestCase):
    def test_delta_n(self):
        # 8.11 §4.4 数值表验证
        self.assertAlmostEqual(galaxy.delta_n(3), 0.154, delta=0.001)
        self.assertAlmostEqual(galaxy.delta_n(4), 0.100, delta=0.001)
        self.assertAlmostEqual(galaxy.delta_n(5), 0.069, delta=0.001)

    def test_shell_radius(self):
        self.assertAlmostEqual(galaxy.shell_radius(1), 8.96, delta=0.01)  # 内层重构区表值
        self.assertAlmostEqual(galaxy.shell_radius(2), 11.7, delta=0.1)
        self.assertAlmostEqual(galaxy.shell_radius(3), 12.2, delta=0.1)
        self.assertAlmostEqual(galaxy.shell_radius(4), 12.1, delta=0.1)
        # 饱和：n 大 → R_c ≈ 11.7
        self.assertAlmostEqual(galaxy.shell_radius(20), 11.7, delta=0.2)

    def test_fill_factor(self):
        self.assertAlmostEqual(galaxy.shell_fill_factor(), 1.5)

    def test_transition_radius(self):
        self.assertAlmostEqual(galaxy.transition_radius_kpc(), 5.85, delta=0.05)

    def test_bulge_sigma(self):
        self.assertAlmostEqual(galaxy.sigma_bulge() / 1e3, 100.0)


class TestSigmaProfile(unittest.TestCase):
    def test_bulge_region(self):
        s = distribution.sigma_profile(0.5 * KPC, M_MW)
        self.assertAlmostEqual(s, 100.0e3)  # 核球区 σ≈100 km/s

    def test_disk_region(self):
        s = distribution.sigma_profile(R_SUN, M_MW)
        v = __import__('cosmogeo.potential', fromlist=['v_circ']).v_circ(R_SUN, M_MW)
        self.assertAlmostEqual(s, v / math.sqrt(2))

    def test_halo_locked_sigma(self):
        # 晕区（>12 kpc）速度锁定：v_c⁴ = G_eff·M·a_0 常数 → σ 常数（8.11 §7.1）
        s1 = distribution.sigma_profile(15 * KPC, M_MW)
        s2 = distribution.sigma_profile(30 * KPC, M_MW)
        self.assertAlmostEqual(s1, s2, delta=1e-6, msg="晕区色散应锁定为常数（v_c⁴ 常数）")


class TestOsipkovMerritt(unittest.TestCase):
    def test_beta_limits(self):
        b_in = distribution.osipkov_merritt_beta(1e3 * KPC * 0, M_MW)  # r→0
        self.assertAlmostEqual(distribution.osipkov_merritt_beta(1e-6, m=M_MW), 0.0, delta=1e-6)
        self.assertGreater(distribution.osipkov_merritt_beta(30 * KPC, m=M_MW), 0.8,
                           "大半径 β → 1（切向各向异性）")

    def test_beta_half(self):
        # r = r_a 处 β = 1/2
        r_m = distribution.potential_permeation_radius(M_MW)
        self.assertAlmostEqual(distribution.osipkov_merritt_beta(r_m, r_a=r_m), 0.5, delta=1e-6)


class TestDF(unittest.TestCase):
    def test_gaussian_normalized(self):
        # 高斯 DF 在 σ 内积分 ≈ 1（数值）
        sig = 100e3
        n = 500
        total = 0.0
        for i in range(n):
            v = sig * 6 * i / n
            total += df.gaussian_velocity_df(v, sig) * 4 * math.pi * v * v * (sig * 6 / n)
        self.assertAlmostEqual(total, 1.0, delta=0.01)

    def test_escape_truncation(self):
        v_esc_r = __import__('cosmogeo.potential', fromlist=['v_esc']).v_esc(R_SUN, M_MW)
        self.assertEqual(df.velocity_distribution(v_esc_r * 1.1, R_SUN, M_MW), 0.0)

    def test_peak_velocity(self):
        m = df.maxwell_from_sigma(100e3)
        self.assertAlmostEqual(m["v_peak"], math.sqrt(2) * 100e3)
        self.assertAlmostEqual(m["v_rms"], math.sqrt(3) * 100e3)

    def test_circular_pdf(self):
        c = df.circular_velocity_pdf(R_SUN, M_MW)
        self.assertTrue(150e3 < c["mean_v"] < 250e3, f"圆速度均值应在 150-250 km/s")


if __name__ == "__main__":
    unittest.main()

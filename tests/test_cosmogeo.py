# -*- coding: utf-8 -*-
"""cosmogeo 单元测试（对照文章数值）."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
import math

from cosmogeo import rotation, hubble, cmb, lensing, constants


class TestHubble(unittest.TestCase):
    def test_h0_geo(self):
        h = hubble.h0_geo()
        self.assertAlmostEqual(h, 68.9, delta=1.0, msg="H_0^geo ≈ 68.9 km/s/Mpc")

    def test_h0_between_observations(self):
        h = hubble.h0_geo()
        self.assertTrue(constants.H0_PLANCK < h < constants.H0_SH0ES,
                        f"H_0^geo={h} 应介于 Planck 67.4 与 SH0ES 73.0 之间")

    def test_lambda_res(self):
        self.assertAlmostEqual(constants.LAMBDA_RES, 1.66e-52, delta=0.1e-52,
                               msg="Λ_res ≈ 1.66e-52 m⁻²")

    def test_tau_m(self):
        self.assertAlmostEqual(hubble.tau_m_years() / 1e12, 2.46, delta=0.2,
                               msg="τ_M ≈ 2.46e12 年")

    def test_four_stage(self):
        epochs = hubble.four_stage_epochs()
        self.assertAlmostEqual(epochs["formation_days"], 7.24, delta=0.1)
        self.assertAlmostEqual(epochs["gap_days"], 50.7, delta=0.5, msg="间隙阶段 ≈ 50.7 日")
        self.assertAlmostEqual(epochs["maha_cycle_years"] / 1e9, 11.41, delta=0.5,
                               msg="大周期 ≈ 114.1 亿年")


class TestRotation(unittest.TestCase):
    def test_g_eff(self):
        self.assertAlmostEqual(rotation.g_eff() / constants.G, 1.2106, delta=1e-3)

    def test_newton_regime(self):
        a = rotation.modified_accel(1e-7)  # a_N ≫ a_0
        self.assertAlmostEqual(a / 1e-7, 1.0, delta=0.01)

    def test_deep_mond_regime(self):
        a = rotation.modified_accel(1e-13)  # a_N ≪ a_0
        expected = math.sqrt(1e-13 * constants.A0)
        self.assertAlmostEqual(a / expected, 1.0, delta=0.01)

    def test_flat_limit(self):
        v = rotation.flat_limit_velocity(1.0e42)
        expected = (rotation.g_eff() * 1.0e42 * constants.A0) ** 0.25
        self.assertAlmostEqual(v, expected)


class TestCMB(unittest.TestCase):
    def test_first_peak(self):
        self.assertEqual(cmb.first_peak_l(), 220.0)

    def test_ns(self):
        self.assertAlmostEqual(cmb.n_s(), 0.965, delta=1e-6)

    def test_r_zero(self):
        self.assertEqual(cmb.tensor_to_scalar_ratio(), 0.0)

    def test_sound_speed(self):
        cs = cmb.sound_speed(0.62)
        expected = constants.C / math.sqrt(3 * 1.62)
        self.assertAlmostEqual(cs, expected)


class TestLensing(unittest.TestCase):
    def test_deflection(self):
        alpha = lensing.deflection_angle(1e41, 5e20)
        base = 4 * constants.G * 1e41 / (constants.C ** 2 * 5e20)
        self.assertAlmostEqual(alpha, base * (1 + constants.SIGMA_C))

    def test_apparent_dm(self):
        self.assertAlmostEqual(lensing.apparent_dm_mass(1e41), 1e41 * constants.SIGMA_C)


if __name__ == "__main__":
    unittest.main()

# -*- coding: utf-8 -*-
"""aainv 完整化 + 盘 DF 精确化测试（0.10.0）."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
import math

from cosmogeo import actionangle, df, potential
from cosmogeo.constants import A0

KPC = 3.0857e19
M_MW = 1.2e41
R_SUN = 8.2 * KPC


class TestFrequencies(unittest.TestCase):
    def test_frequency_ratios(self):
        f = actionangle.frequencies(R_SUN, M_MW)
        # 渗透势过渡区：κ/Ω ∈ (1, √2)（Kepler=1、对数势=√2）
        ratio = f["Omega_R"] / f["Omega_phi"]
        self.assertTrue(1.0 < ratio < math.sqrt(2) + 0.01,
                        f"κ/Ω={ratio:.3f} 应介于 1 与 √2")
        # ν = 3Ω（薄盘近似）
        self.assertAlmostEqual(f["Omega_z"], 3 * f["Omega_phi"], delta=1e-9)


class TestActionAngles(unittest.TestCase):
    def test_circular(self):
        v0 = potential.v_circ(R_SUN, M_MW)
        aa = actionangle.action_angles(R_SUN, 0, v0, 0, M_MW)
        self.assertAlmostEqual(aa["JR"], 0.0, delta=aa["Jphi"] * 1e-6, msg="圆轨道 JR≈0")
        self.assertAlmostEqual(aa["Jphi"], R_SUN * v0)
        self.assertAlmostEqual(aa["Jz"], 0.0, delta=aa["Jphi"] * 1e-6)

    def test_eccentric(self):
        v0 = potential.v_circ(R_SUN, M_MW)
        aa = actionangle.action_angles(R_SUN, 0.3 * v0, 0.7 * v0, 0, M_MW)
        self.assertGreater(aa["JR"], 0)
        # J_φ 守恒 = r·vφ
        self.assertAlmostEqual(aa["Jphi"], R_SUN * 0.7 * v0)

    def test_vertical(self):
        v0 = potential.v_circ(R_SUN, M_MW)
        nu = actionangle.vertical_frequency(R_SUN, M_MW)
        aa = actionangle.action_angles(R_SUN, 0, v0, 20e3, M_MW, z=0.2 * KPC)
        e_z = 0.5 * (20e3) ** 2 + 0.5 * nu * nu * (0.2 * KPC) ** 2
        self.assertAlmostEqual(aa["Jz"], e_z / nu, delta=e_z / nu * 1e-9, msg="谐波 J_z=E_z/ν")


class TestQuasiIsothermalDF(unittest.TestCase):
    def test_monotonic(self):
        v0 = potential.v_circ(R_SUN, M_MW)
        f_circ = df.quasiisothermal_df(R_SUN, 0, v0, 0, M_MW)
        f_ecc = df.quasiisothermal_df(R_SUN, 0.3 * v0, 0.7 * v0, 0, M_MW)
        f_hot = df.quasiisothermal_df(R_SUN, 0.3 * v0, 0.7 * v0, 50e3, M_MW, z=0.3 * KPC)
        self.assertGreater(f_circ, f_ecc, "圆轨道 DF 应最大")
        self.assertGreater(f_ecc, f_hot, "热粒子 DF 应更小")

    def test_positive(self):
        v0 = potential.v_circ(R_SUN, M_MW)
        f = df.quasiisothermal_df(R_SUN, 0, v0, 0, M_MW)
        self.assertGreater(f, 0)

    def test_sampling(self):
        s = df.disk_df_sample(R_SUN, M_MW, n=10, seed=3)
        self.assertEqual(s["n"], 10)
        self.assertEqual(len(s["samples"]), 10)
        v0 = potential.v_circ(R_SUN, M_MW)
        for p in s["samples"]:
            self.assertAlmostEqual(p["r"], R_SUN)
            self.assertTrue(abs(p["vphi"]) < 3 * v0, "采样切向速度应在合理范围")


if __name__ == "__main__":
    unittest.main()

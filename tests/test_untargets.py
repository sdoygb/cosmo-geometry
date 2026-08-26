# -*- coding: utf-8 -*-
"""未达目标实现测试：矮星系 core/太阳锁定/径向作用量."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
import math

from cosmogeo import dwarf, galaxy, orbit, potential
from cosmogeo.constants import A0, G

KPC = 3.0857e19
M_MW = 1.2e41
M_DWARF = 1e9 * 1.989e30   # 10^9 M_sun
R_SUN = 8.2 * KPC


class TestDwarf(unittest.TestCase):
    def test_core_flat(self):
        # 矮星系深 MOND 区（r > r_M≈1.17 kpc）：v 常数（core 型）
        v1 = dwarf.dwarf_velocity(10.0 * KPC, M_DWARF)
        v2 = dwarf.dwarf_velocity(100.0 * KPC, M_DWARF)
        self.assertAlmostEqual(v1 / v2, 1.0, delta=0.02, msg="矮星系深 MOND 区 v 应平坦")

    def test_cusp_core_signature(self):
        # 中心区（0.05–0.3 kpc ≪ r_s=1 kpc）：NFW cusp v ∝ √r 上升（slope>0），
        # 几何论无暗物质中心牛顿下降（slope<0）——cusp vs core 的可检验区分
        sig = dwarf.cusp_core_signature(M_DWARF, r_in=0.05 * KPC, r_out=0.3 * KPC,
                                        r_s=1 * KPC, m_halo=3e10 * 1.989e30)
        self.assertGreater(sig["slope_nfw"], 0.3, "NFW 中心 cusp 斜率应 >0.3（上升）")
        self.assertLess(sig["slope_geo"], 0.0, "几何论中心斜率应 <0（无 cusp，下降/平坦）")

    def test_deep_mond_radius(self):
        r_dm = dwarf.deep_mond_radius(M_DWARF)
        self.assertGreater(r_dm, 0)


class TestSolarLock(unittest.TestCase):
    def test_lock_ratio(self):
        sl = galaxy.solar_lock_position()
        self.assertAlmostEqual(sl["ratio"], 1.0, delta=0.05,
                               msg=f"r_crit/r_sun 应 ≈1（{sl['ratio']:.3f}）")
        self.assertAlmostEqual(sl["r_crit_kpc"], 8.3, delta=0.5)


class TestRadialAction(unittest.TestCase):
    def test_circular_action_zero(self):
        v0 = potential.v_circ(R_SUN, M_MW)
        L = R_SUN * v0
        E = -abs(potential.potential(R_SUN, M_MW)) + 0.5 * v0 ** 2
        Jr = orbit.radial_action(E, L, M_MW)
        self.assertLess(Jr, L * 1e-6, "圆轨道 J_r 应 ≈0")

    def test_eccentric_action_positive(self):
        v0 = potential.v_circ(R_SUN, M_MW)
        L2 = 0.7 * R_SUN * v0
        E2 = -abs(potential.potential(R_SUN, M_MW)) + 0.5 * (0.3 * v0) ** 2 + 0.5 * (L2 / R_SUN) ** 2
        Jr = orbit.radial_action(E2, L2, M_MW)
        self.assertGreater(Jr, 0, "偏心轨道 J_r 应 >0")
        # 偏心越大 J_r 越大
        L3 = 0.3 * R_SUN * v0
        E3 = -abs(potential.potential(R_SUN, M_MW)) + 0.5 * (0.7 * v0) ** 2 + 0.5 * (L3 / R_SUN) ** 2
        Jr3 = orbit.radial_action(E3, L3, M_MW)
        self.assertGreater(Jr3, Jr, "更偏心轨道 J_r 应更大")


if __name__ == "__main__":
    unittest.main()

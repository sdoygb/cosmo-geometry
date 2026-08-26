# -*- coding: utf-8 -*-
"""galpy 功能覆盖测试：势家族/质量分布/轨道诊断（0.8.0）."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
import math

from cosmogeo import potential, orbit
from cosmogeo.constants import A0, G

KPC = 3.0857e19
M_MW = 1.2e41
R_SUN = 8.2 * KPC
YEAR = 3.15576e7


class TestMassProfile(unittest.TestCase):
    def test_mass_within_increasing(self):
        m1 = potential.mass_within(1 * KPC, M_MW)
        m8 = potential.mass_within(8.2 * KPC, M_MW)
        m30 = potential.mass_within(30 * KPC, M_MW)
        self.assertLess(m1, m8, "M(<r) 应随 r 增长")
        self.assertLess(m8, m30)
        # 量级：10^10-10^11 M_sun
        self.assertTrue(1e10 < m8 / 1.989e30 < 1e12, f"M(<8.2kpc)={m8/1.989e30:.1e} M_sun")


class TestPotentials(unittest.TestCase):
    def test_kepler(self):
        phi = potential.kepler_potential(R_SUN, M_MW)
        self.assertLess(phi, 0)

    def test_logarithmic_flat(self):
        # 对数势 Φ=v_c²·ln(r)：v_circ² = r·dΦ/dr = v_c² 常数（平坦旋转曲线）
        v_c = 200e3
        # 数值导数验证：v_circ(r) = v_c 对任意 r
        for r_kpc in (1, 10, 50):
            r = r_kpc * KPC
            eps = r * 1e-4
            dphi_dr = ((potential.logarithmic_potential(r + eps, v_c)
                        - potential.logarithmic_potential(r - eps, v_c)) / (2 * eps))
            v_circ = (r * dphi_dr) ** 0.5
            self.assertAlmostEqual(v_circ, v_c, delta=1.0,
                                   msg=f"对数势 v_circ({r_kpc} kpc) 应=v_c")

    def test_plummer(self):
        phi = potential.plummer_potential(R_SUN, M_MW, b=1 * KPC)
        self.assertLess(phi, 0)
        # 中心密度
        rho = potential.plummer_density(0.1 * KPC, M_MW, 1 * KPC)
        self.assertGreater(rho, 0)


class TestOrbitDiagnostics(unittest.TestCase):
    def _integrate(self, vr_frac, vt_frac, t_frac=2.5):
        v0 = potential.v_circ(R_SUN, M_MW)
        t_total = orbit.orbital_period(R_SUN, M_MW) * YEAR * t_frac
        dt = t_total / 4000
        return orbit.integrate_orbit(R_SUN, vr_frac * v0, vt_frac * v0, M_MW,
                                     t_total=t_total, dt=dt)

    def test_circular(self):
        orb = self._integrate(0.0, 1.0)
        d = orbit.orbit_diagnostics(orb)
        self.assertLess(d["ecc"], 0.05, "圆轨道 e≈0")

    def test_eccentric(self):
        orb = self._integrate(0.4, 0.8)
        d = orbit.orbit_diagnostics(orb)
        self.assertGreater(d["ecc"], 0.1, "偏心轨道 e>0.1")
        self.assertLess(d["ecc"], 0.9)

    def test_classify(self):
        self.assertEqual(orbit.classify_orbit(0.01), "圆轨道")
        self.assertEqual(orbit.classify_orbit(0.5), "椭圆轨道")
        self.assertEqual(orbit.classify_orbit(None), "逃逸/未束缚")

    def test_j_phi(self):
        v0 = potential.v_circ(R_SUN, M_MW)
        self.assertAlmostEqual(orbit.j_phi(R_SUN, v0), R_SUN * v0)

    def test_vertical_action(self):
        Jz = orbit.vertical_action_approx(R_SUN, M_MW, 0.1 * KPC, 30e3)
        self.assertGreater(Jz, 0)


if __name__ == "__main__":
    unittest.main()

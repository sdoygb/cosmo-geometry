# -*- coding: utf-8 -*-
"""3D 轨道 + 恒星流测试（0.9.0）."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
import math

from cosmogeo import orbit, potential, stream
from cosmogeo.constants import A0

KPC = 3.0857e19
M_MW = 1.2e41
R_SUN = 8.2 * KPC
YEAR = 3.15576e7


class TestOrbit3D(unittest.TestCase):
    def _orb3d(self, z0=0.3 * KPC, vz0=30e3, t_frac=0.5):
        v0 = potential.v_circ(R_SUN, M_MW)
        t_total = orbit.orbital_period(R_SUN, M_MW) * YEAR * t_frac
        dt = t_total / 1500
        return orbit.integrate_orbit_3d(R_SUN, z0, 0.0, v0, vz0, M_MW,
                                        t_total=t_total, dt=dt)

    def test_radial_circular(self):
        d = orbit.orbit_diagnostics_3d(self._orb3d())
        self.assertLess(d["ecc"], 0.05, "3D 圆轨道径向 e≈0")

    def test_vertical_oscillation(self):
        orb = self._orb3d()
        zs = orb["z"]
        self.assertGreater(max(zs), 0, "应有垂直振荡")
        self.assertAlmostEqual(max(zs), -min(zs), delta=max(zs) * 0.05,
                               msg="垂直振荡应对称")

    def test_zmax(self):
        d = orbit.orbit_diagnostics_3d(self._orb3d())
        self.assertGreater(d["zmax"], 0)
        self.assertLess(d["zmax"], 2 * KPC, "zmax 应有限（束缚）")

    def test_angular_momentum(self):
        orb = self._orb3d()
        # vphi·R 守恒
        L0 = R_SUN * orb["vphi"][0]
        L_mid = orb["R"][len(orb["R"]) // 2] * orb["vphi"][len(orb["R"]) // 2]
        self.assertAlmostEqual(L0, L_mid, delta=L0 * 1e-6)


class TestStream(unittest.TestCase):
    def test_analytic_length(self):
        L = stream.stream_analytic_length(1e9, 10e3)
        self.assertAlmostEqual(L, 10e3 * 1e9 * YEAR / KPC, delta=0.01)

    def test_angle_spread_formula(self):
        an = stream.stream_angle_spread(R_SUN, M_MW, 1e8)
        self.assertGreater(an["angle_span_deg"], 0)
        self.assertGreater(an["length_kpc"], 0)

    def test_numeric_stream(self):
        res = stream.evolve_stream(R_SUN, M_MW, t_years=2.2e7,
                                   n_particles=20, sigma_v_frac=0.01, steps=80)
        self.assertGreater(res["length_kpc"], 0)
        self.assertGreater(res["width_kpc"], 0)
        self.assertEqual(len(res["angles_deg"]), 20)

    def test_numeric_matches_analytic_order(self):
        # 小弥散线性极限：数值/解析 同数量级（0.5-5 倍）
        an = stream.stream_angle_spread(R_SUN, M_MW, t_years=2.2e7, sigma_v_frac=0.01)
        res = stream.evolve_stream(R_SUN, M_MW, t_years=2.2e7,
                                   n_particles=30, sigma_v_frac=0.01, steps=100)
        angles = res["angles_deg"]
        span = max(angles) - min(angles)
        ratio = span / an["angle_span_deg"]
        self.assertTrue(0.5 < ratio < 5.0,
                        f"数值/解析={ratio:.2f}（应同数量级）")

    def test_stream_grows_with_time(self):
        short = stream.evolve_stream(R_SUN, M_MW, t_years=2.2e7,
                                     n_particles=20, sigma_v_frac=0.01, steps=80)
        long = stream.evolve_stream(R_SUN, M_MW, t_years=4.4e7,
                                    n_particles=20, sigma_v_frac=0.01, steps=80)
        a_short = max(short["angles_deg"]) - min(short["angles_deg"])
        a_long = max(long["angles_deg"]) - min(long["angles_deg"])
        self.assertGreater(a_long, a_short, "流随时间展开")


if __name__ == "__main__":
    unittest.main()

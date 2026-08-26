# -*- coding: utf-8 -*-
"""galpy 对照模块单元测试（galaxy/potential/distribution/orbit）."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import unittest
import math

from cosmogeo import galaxy, potential, distribution, orbit
from cosmogeo.constants import A0, G, YEAR_S

M_MW = 1.2e41          # kg：银河系重子质量量级
R_SUN = 8.2e3 * 3.0857e16  # m：太阳位置 8.2 kpc


class TestGalaxy(unittest.TestCase):
    def test_v_flat(self):
        self.assertAlmostEqual(galaxy.v_flat_km_s(), 235.0, delta=1.0)

    def test_scale_radius(self):
        self.assertAlmostEqual(galaxy.scale_radius_kpc(), 11.7, delta=0.1)

    def test_shell_parameters(self):
        self.assertEqual(galaxy.shell_exponent(), 2.0)
        self.assertAlmostEqual(galaxy.shell_correction(), 0.25)

    def test_critical_radius(self):
        r = galaxy.critical_radius(M_MW)
        self.assertAlmostEqual(r / 3.0857e19, 8.3, delta=1.0, msg="r_crit ≈ 8.3 kpc")

    def test_baryonic_mass(self):
        m = galaxy.baryonic_mass_for_v_flat()
        self.assertTrue(5e10 < m / 1.989e30 < 3e11, f"M 应在银河系量级 1e11±，实际 {m/1.989e30:.1e}")


class TestPotential(unittest.TestCase):
    def test_v_circ_solar(self):
        v = potential.v_circ(R_SUN, M_MW)
        self.assertTrue(200e3 < v < 260e3, f"太阳位置 v_circ={v/1e3:.0f} km/s 应在 200-260")

    def test_permeation_radius(self):
        r_m = potential.permeation_radius(M_MW)
        self.assertTrue(7e19 < r_m < 3.2e20, f"r_M 应在 ~9 kpc（{r_m/3.0857e19:.1f}）")

    def test_v_esc_bounded(self):
        v = potential.v_esc(R_SUN, M_MW, r_boundary=100 * 3.0857e19)
        self.assertTrue(300e3 < v < 600e3, f"v_esc 应在 300-600 km/s，实际 {v/1e3:.0f}")

    def test_potential_negative(self):
        self.assertLess(potential.potential(R_SUN, M_MW), 0)

    def test_circular_freq(self):
        omega = potential.circular_angular_freq(R_SUN, M_MW)
        T = 2 * math.pi / omega / YEAR_S
        self.assertTrue(200e6 < T < 260e6, f"轨道周期 ≈ 230 Myr，实际 {T/1e6:.0f}")


class TestDistribution(unittest.TestCase):
    def test_sigma_iso_relation(self):
        s = distribution.sigma_iso(R_SUN, M_MW)
        v = potential.v_circ(R_SUN, M_MW)
        self.assertAlmostEqual(s, v / math.sqrt(2))

    def test_beta(self):
        self.assertAlmostEqual(distribution.anisotropy_beta(0.75), 0.25)
        self.assertAlmostEqual(distribution.anisotropy_beta(1.0), 0.0)

    def test_sigma_t_relation(self):
        s_r = distribution.sigma_r(R_SUN, M_MW)
        s_t = distribution.sigma_t(R_SUN, M_MW, beta=0.25)
        self.assertAlmostEqual(s_t, s_r * math.sqrt(0.75))


class TestOrbit(unittest.TestCase):
    def test_circular_action(self):
        L = orbit.angular_momentum(R_SUN, M_MW)
        self.assertAlmostEqual(orbit.circular_action(R_SUN, M_MW), L)

    def test_circular_orbit_stable(self):
        vt = potential.v_circ(R_SUN, M_MW)
        orb = orbit.integrate_orbit(R_SUN, 0.0, vt, M_MW, t_total=2 * YEAR_S, dt=0.01 * YEAR_S)
        ratio = max(orb["r"]) / min(orb["r"])
        self.assertLess(ratio, 1.01, f"圆轨道应稳定（r 波动 {ratio:.5f}）")

    def test_escape_condition(self):
        vt = potential.v_circ(R_SUN, M_MW)
        self.assertFalse(orbit.escape_condition(R_SUN, 0.0, vt, M_MW))
        v_esc = potential.v_esc(R_SUN, M_MW)
        self.assertTrue(orbit.escape_condition(R_SUN, 0.0, v_esc * 1.05, M_MW))


if __name__ == "__main__":
    unittest.main()

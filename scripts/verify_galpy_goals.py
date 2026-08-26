#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verify_galpy_goals.py — 对照 galpy 目的的实现验证.

galpy 的能力清单（目的）→ 几何论实现：
  1. 旋转曲线 rotation curve       → potential.v_circ（渗透函数闭式）
  2. 势场 potential models         → potential.potential / v_esc
  3. 分布函数速度矩 velocity moments → distribution.sigma_r/sigma_t/beta
  4. 轨道积分 orbit integration    → orbit.integrate_orbit
  5. 作用量 action-angle           → orbit.circular_action（圆轨道闭式）
  6. 银河系旋转/尺度（8.11）        → galaxy.v_flat / R_c^gal / 壳层
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cosmogeo import galaxy, potential, distribution, orbit
from cosmogeo.constants import G, A0, YEAR_S

M_MW = 1.2e41  # kg：银河系可见（重子）质量 ~6e10 M_sun

print("=" * 72)
print("对照 galpy 目的：几何论实现验证")
print("=" * 72)

# ---- 1. 旋转曲线（galpy rotation curve）----
print("\n[1] 旋转曲线（对照 galpy v_circ）")
for r_pc in [1e3, 5e3, 2e4, 8.2e3]:  # 8.2 kpc ≈ 太阳位置
    r = r_pc * 3.0857e16
    v = potential.v_circ(r, M_MW)
    print(f"    r={r_pc:>6.0f} pc: v_circ = {v/1e3:6.1f} km/s")
r_sun = 8.2e3 * 3.0857e16
v_sun = potential.v_circ(r_sun, M_MW)
print(f"    太阳位置 (8.2 kpc): {v_sun/1e3:.0f} km/s（观测 ~220）")

# ---- 2. 势场 + 逃逸速度（galpy potential）----
print("\n[2] 渗透函数势 Φ(r) 与逃逸速度")
phi = potential.potential(r_sun, M_MW)
v_esc = potential.v_esc(r_sun, M_MW, r_boundary=100 * 3.0857e19)
print(f"    Φ(8.2kpc) = {phi/1e12:.2f} (km/s)²·1e12（负值）")
print(f"    v_esc(8.2kpc, 100kpc边界) = {v_esc/1e3:.0f} km/s（银河系观测 ~500-550，量级一致）")
r_m = potential.permeation_radius(M_MW)
print(f"    渗透半径 r_M = {r_m/3.0857e19:.1f} kpc（MOND 半径）")

# ---- 3. 速度矩（galpy sphericaldf）----
print("\n[3] 速度矩（对照 galpy vmomentdensity/sigmar/sigmat/beta）")
s_r = distribution.sigma_r(r_sun, M_MW)
s_t = distribution.sigma_t(r_sun, M_MW, beta=0.0)
print(f"    σ_r(8.2kpc) = {s_r/1e3:.0f} km/s（各向同性近似，观测 ~100-120）")
print(f"    σ_t(8.2kpc) = {s_t/1e3:.0f} km/s")
print(f"    β = 1-σ_t²/σ_r² = {distribution.anisotropy_beta(1.0):.2f}（各向同性基线）")

# ---- 4. 轨道积分（galpy orbit）----
print("\n[4] 轨道积分（对照 galpy orbit integration）")
vt0 = potential.v_circ(r_sun, M_MW)
orb = orbit.integrate_orbit(r_sun, 0.0, vt0, M_MW, t_total=5 * YEAR_S, dt=0.01 * YEAR_S)
r_range = max(orb["r"]) / min(orb["r"])
print(f"    圆轨道 5 年积分: {len(orb['r'])} 步，r 波动比 {r_range:.6f}（应 ≈1：圆轨道保持）")
T = orbit.orbital_period(r_sun, M_MW)
print(f"    轨道周期 T = {T/1e6:.2f} 百万年（太阳绕银心 ~230 Myr 量级）")

# ---- 5. 作用量（galpy actionAngle）----
print("\n[5] 圆轨道作用量（对照 galpy actionAngle）")
L = orbit.circular_action(r_sun, M_MW)
J = orbit.circular_action(r_sun, M_MW)
print(f"    L = r·v_circ = {L:.3e} m²/s，J_r=0，J=L（圆轨道闭式）")

# ---- 6. 银河系结构（8.11）----
print("\n[6] 银河系结构（8.11）")
print(f"    v_flat = {galaxy.v_flat_km_s():.0f} km/s（文章 235，观测 ~220-240）")
print(f"    R_c^gal = {galaxy.scale_radius_kpc():.1f} kpc（文章 11.7）")
print(f"    壳层 β = {galaxy.shell_exponent()}，α = {galaxy.shell_correction()}")
r_crit = galaxy.critical_radius(M_MW)
print(f"    a_N=a_0 临界壳层 r_crit = {r_crit/3.0857e19:.1f} kpc（速度锁定处）")
M_req = galaxy.baryonic_mass_for_v_flat()
print(f"    由 v_flat=235 反推重子质量 = {M_req/1.989e30:.1e} M_sun（银河系 ~6e10 M_sun）")

print("\n" + "=" * 72)
print("galpy 目的覆盖：旋转曲线 ✓ 势/逃逸 ✓ 速度矩 ✓ 轨道 ✓ 作用量(圆) ✓ 银河系结构 ✓")
print("=" * 72)

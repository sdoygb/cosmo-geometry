#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""quickstart.py — 快速上手：观测锚定 → 动态演化 + 静态计算.

实用主义路线：不做第一性原理推导，全部用现有观测数据锚定，
然后直接算动态（时间演化）与静态（闭式量）。
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cosmogeo import anchors, hubble, satellite, solar, galaxy, expansion, cmb, orbit, potential

print("=" * 72)
print("cosmo-geometry 快速上手（观测锚定版）")
print("=" * 72)

print("\n[1] 观测锚定一览")
for k, v in anchors.all_anchors().items():
    print(f"    {k}: {v}")

print("\n[2] 静态计算（闭式，观测锚定输入）")
print(f"    地月周期比 5α = {satellite.period_ratio_theory():.6f}（观测 {satellite.period_ratio_observed():.6f}）")
print(f"    水星进动 = {solar.mercury_precession_per_century():.1f} 角秒/世纪（观测 {anchors.MERCURY_PRECESSION_ARCSEC_CENTURY}）")
print(f"    银河系 v_flat = {galaxy.v_flat_km_s():.0f} km/s（观测 {anchors.V_FLAT_OBS_KM_S}）")
print(f"    太阳位置锁定 r_crit = {galaxy.solar_lock_position()['r_crit_kpc']:.1f} kpc（观测 r_sun = {anchors.SUN_R_KPC}）")

print("\n[3] 动态演化（时间依赖）")
print(f"    宇宙年龄 t(0) = {hubble.age_universe_gyr():.1f} Gyr（观测 H0={anchors.H0_OBS_KM_S_MPC} 对应 ~14）")
print(f"    地月 10 亿年前 = {satellite.moon_distance_past(1000):.0f} km")
print(f"    水星未来 10 世纪进动 = {solar.mercury_precession_angle(10):.1f} 角秒")
print(f"    H(z=1) = {expansion.hubble_hz(1.0):.1f} km/s/Mpc")
print(f"    T(z=1100) = {cmb.temperature_at_redshift(1100):.0f} K（复合）")

print("\n[4] 轨道（动态积分，观测锚定）")
v0 = potential.v_circ(anchors.SUN_R_KPC * 3.0857e19, 1.2e41)
print(f"    太阳位置圆速度 = {v0/1e3:.0f} km/s（观测 {anchors.V_SUN_KM_S}）")
T = orbit.orbital_period(anchors.SUN_R_KPC * 3.0857e19, 1.2e41)
print(f"    轨道周期 = {T/1e6:.0f} 百万年（观测 ~230）")

print("\n" + "=" * 72)
print("所有计算仅需观测锚定（anchors 模块），无需理论闭环。")
print("=" * 72)

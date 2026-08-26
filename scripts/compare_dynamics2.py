#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""compare_dynamics2.py — 抽查两个动态演化：逃逸速度 + 垂直振荡.

几何论（渗透势，无暗物质）vs galpy MWPotential2014（含暗物质）：
  1. 逃逸速度 v_esc(8.2 kpc)：动态轨道理论量（galpy vesc 函数）
  2. 3D 垂直振荡：z₀=0.5 kpc 扰动轨道的 z(t)——垂直周期/频率对比
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.sans-serif"] = ["PingFang SC", "Heiti TC", "Arial Unicode MS", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

import astropy.units as u
from galpy.potential import MWPotential2014, vcirc as gvcirc, vesc as gvesc
from galpy.orbit import Orbit

from cosmogeo import potential as cpot
from cosmogeo import orbit as corbit

KPC = 3.0857e19
YEAR = 3.15576e7
G_RO, G_VO = 8.0, 220.0
R_SUN_KPC = 8.2

# ---- 几何论质量校准（与 compare_galpy.py 一致）----
def v_geo(r_kpc, M):
    return cpot.v_circ(r_kpc * KPC, M) / 1e3
target = gvcirc(MWPotential2014, R_SUN_KPC / G_RO) * G_VO
lo, hi = 1e40, 1e42
for _ in range(60):
    mid = (lo + hi) / 2
    if v_geo(R_SUN_KPC, mid) > target:
        hi = mid
    else:
        lo = mid
M_CAL = (lo + hi) / 2
print(f"几何论校准 M = {M_CAL:.3e} kg（v(8.2kpc)={v_geo(R_SUN_KPC, M_CAL):.1f} km/s）")

# ============ 1. 逃逸速度对比 ============
print("\n[1] 逃逸速度 v_esc(8.2 kpc)")
# galpy（物理单位 vesc：返回无量纲需 ×vo）
v_esc_galpy = gvesc(MWPotential2014, R_SUN_KPC / G_RO) * G_VO  # km/s
print(f"  galpy v_esc(8.2kpc) = {v_esc_galpy:.0f} km/s")
# 几何论（100 kpc 有限边界——银河系边界逃逸）
v_esc_geo_inf = cpot.v_esc(R_SUN_KPC * KPC, M_CAL) / 1e3
v_esc_geo_100 = cpot.v_esc(R_SUN_KPC * KPC, M_CAL, r_boundary=100 * KPC) / 1e3
print(f"  几何论 v_esc(8.2kpc, 无穷远) = {v_esc_geo_inf:.0f} km/s")
print(f"  几何论 v_esc(8.2kpc, 100kpc边界) = {v_esc_geo_100:.0f} km/s")
print(f"  观测（Gaia RVS）：~500-550 km/s")

# ============ 2. 垂直振荡对比 ============
print("\n[2] 3D 垂直振荡（z₀=0.5 kpc，vz=0，圆轨道背景）")
Z0_KPC = 0.5

# galpy 3D 轨道（物理单位）
vT_circ = gvcirc(MWPotential2014, R_SUN_KPC / G_RO) * G_VO
o3 = Orbit([R_SUN_KPC * u.kpc, 0 * u.km / u.s, vT_circ * u.km / u.s,
            Z0_KPC * u.kpc, 0 * u.km / u.s, 0 * u.rad])
ts3 = np.linspace(0, 0.3, 1000) * u.Gyr  # 0.3 Gyr（覆盖垂直振荡多周期）
o3.integrate(ts3, MWPotential2014)
z_galpy = np.array(o3.z(ts3))  # kpc

# 几何论 3D 轨道
v0 = cpot.v_circ(R_SUN_KPC * KPC, M_CAL)
o_geo3 = corbit.integrate_orbit_3d(R_SUN_KPC * KPC, Z0_KPC * KPC, 0.0, v0, 0.0, M_CAL,
                                   t_total=0.3e9 * YEAR, dt=1e5 * YEAR)
z_geo = np.array(o_geo3["z"]) / KPC

# 垂直振荡周期检测（z 过零或极值）
def vertical_period(z_arr, dt_years):
    """z(t) 相邻同方向过零的时间间隔 → 垂直周期."""
    crossing = []
    for i in range(1, len(z_arr)):
        if z_arr[i-1] * z_arr[i] < 0:  # 过零
            crossing.append(i)
    if len(crossing) < 2:
        return None
    periods = []
    for k in range(1, len(crossing)):
        d = (crossing[k] - crossing[k-1]) * dt_years
        if d > 0:
            periods.append(d)
    if not periods:
        return None
    return 2.0 * np.mean(periods)  # 相邻过零间隔 ×2 = 全周期

T_z_galpy = vertical_period(z_galpy, 0.3e9 / 1000)
T_z_geo = vertical_period(z_geo, 1e5)
print(f"  galpy 垂直振荡周期 = {T_z_galpy/1e6:.1f} Myr（若测到）")
print(f"  几何论垂直振荡周期 = {T_z_geo/1e6:.1f} Myr（ν_z=3Ω 理论：{2*np.pi/(3*8.976e-16)/YEAR/1e6:.1f} Myr）")
# 垂直幅度
print(f"  galpy z 范围 = [{z_galpy.min():.3f}, {z_galpy.max():.3f}] kpc（初始 0.5）")
print(f"  几何论 z 范围 = [{z_geo.min():.3f}, {z_geo.max():.3f}] kpc（初始 0.5）")

# 图：z(t) 对比
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
t_galpy = ts3.to_value(u.Gyr)
t_geo = np.linspace(0, 0.3, len(z_geo))
ax[0].plot(t_galpy, z_galpy, "b-", lw=1.2)
ax[0].set_title(f"galpy 垂直振荡（周期 {T_z_galpy/1e6 if T_z_galpy else 0:.0f} Myr）")
ax[0].set_xlabel("t (Gyr)"); ax[0].set_ylabel("z (kpc)"); ax[0].set_ylim(-1, 1)
ax[1].plot(t_geo, z_geo, "r-", lw=1.2)
ax[1].set_title(f"几何论垂直振荡（周期 {T_z_geo/1e6 if T_z_geo else 0:.0f} Myr，ν=3Ω）")
ax[1].set_xlabel("t (Gyr)"); ax[1].set_ylabel("z (kpc)"); ax[1].set_ylim(-1, 1)
plt.tight_layout()
plt.savefig("vertical_compare.png", dpi=120)
print("图 vertical_compare.png 已保存")

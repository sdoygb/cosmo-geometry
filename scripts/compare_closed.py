#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""compare_closed.py — 三项闭合后的最终对比图（几何论 vs galpy vs 观测）.

面板：
  1. 旋转曲线：单点模型 / 分布模型（闭合）/ galpy MWPotential2014
  2. v_esc(8.2kpc)：galpy / 几何论分布（观测标定边界）/ 观测带 500-550
  3. 垂直周期：galpy / 几何论盘势（闭合）/ 观测带 80-100 Myr
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.sans-serif"] = ["PingFang SC", "Heiti TC", "Arial Unicode MS", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

from galpy.potential import MWPotential2014, vcirc as gvcirc, vesc as gvesc

from cosmogeo import rotation, actionangle
from cosmogeo.constants import A0, YEAR_S

KPC = 3.0857e19
G_RO, G_VO = 8.0, 220.0
R_SUN_KPC = 8.2

# ---- 几何论质量校准 ----
def v_geo(r_kpc, M):
    return rotation.rotation_velocity(r_kpc * KPC, M) / 1e3
target = gvcirc(MWPotential2014, R_SUN_KPC / G_RO) * G_VO
lo, hi = 1e40, 1e42
for _ in range(60):
    mid = (lo + hi) / 2
    if v_geo(R_SUN_KPC, mid) > target:
        hi = mid
    else:
        lo = mid
M_CAL = (lo + hi) / 2
R8 = R_SUN_KPC * KPC

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# ===== 面板 1：旋转曲线 =====
R_kpc = np.linspace(1, 30, 60)
v_galpy = np.array([gvcirc(MWPotential2014, r / G_RO) * G_VO for r in R_kpc])
v_point = np.array([rotation.rotation_velocity(r * KPC, M_CAL) / 1e3 for r in R_kpc])
v_dist = np.array([rotation.rotation_velocity_distributed(r * KPC, M_CAL) / 1e3 for r in R_kpc])

ax = axes[0]
ax.plot(R_kpc, v_galpy, "b-", lw=2, label="galpy MW（含暗物质）")
ax.plot(R_kpc, v_point, "g--", lw=1.5, alpha=0.6, label="几何论单点（闭合前）")
ax.plot(R_kpc, v_dist, "r-", lw=2, label="几何论分布（闭合后）")
ax.axvline(8.2, color="gray", ls=":", lw=1)
ax.set_xlabel("R (kpc)"); ax.set_ylabel("v_circ (km/s)")
ax.set_title("旋转曲线：闭合后 vs galpy")
ax.legend(fontsize=8); ax.set_xlim(0, 30); ax.set_ylim(0, 300)

# ===== 面板 2：v_esc =====
v_esc_galpy = gvesc(MWPotential2014, R_SUN_KPC / G_RO) * G_VO  # km/s
v_esc_geo = rotation.v_esc_distributed(R8, M_CAL, 395 * KPC) / 1e3
v_esc_old100 = __import__("cosmogeo.potential", fromlist=["v_esc"]).v_esc(R8, M_CAL, r_boundary=100 * KPC) / 1e3

ax = axes[1]
bars = ax.bar(["galpy", "几何论分布\n(闭合)", "几何论单点\n(100kpc)"],
              [v_esc_galpy, v_esc_geo, v_esc_old100],
              color=["#4472C4", "#ED7D31", "#A5A5A5"])
ax.axhspan(500, 550, color="green", alpha=0.15, label="观测 500-550")
for b, v in zip(bars, [v_esc_galpy, v_esc_geo, v_esc_old100]):
    ax.text(b.get_x() + b.get_width() / 2, v + 8, f"{v:.0f}", ha="center", fontsize=10)
ax.set_ylabel("v_esc (km/s)")
ax.set_title("逃逸速度 v_esc(8.2 kpc)")
ax.legend(fontsize=8); ax.set_ylim(0, 650)

# ===== 面板 3：垂直周期 =====
nu_geo = actionangle.vertical_frequency(R8, M_CAL)
T_z_geo = 2 * math.pi / nu_geo / YEAR_S / 1e6  # Myr
T_z_galpy = 110.0  # 从 3D 轨道实测（compare_dynamics2.py）
T_z_old3w = 76.5

ax = axes[2]
bars = ax.bar(["galpy\n(实测)", "几何论盘势\n(闭合)", "几何论 3Ω\n(旧)"],
              [T_z_galpy, T_z_geo, T_z_old3w],
              color=["#4472C4", "#ED7D31", "#A5A5A5"])
ax.axhspan(80, 100, color="green", alpha=0.15, label="观测太阳邻域 80-100")
for b, v in zip(bars, [T_z_galpy, T_z_geo, T_z_old3w]):
    ax.text(b.get_x() + b.get_width() / 2, v + 2, f"{v:.0f}", ha="center", fontsize=10)
ax.set_ylabel("垂直周期 (Myr)")
ax.set_title("垂直振荡周期 T_z")
ax.legend(fontsize=8); ax.set_ylim(0, 140)

plt.tight_layout()
plt.savefig("closed_compare.png", dpi=130)
print("closed_compare.png 已保存")

# 数值摘要
print(f"\n旋转曲线（r kpc: galpy / 分布 / 单点）:")
for r in [1, 4, 8.2, 15, 30]:
    i = int((r - 1) / 29 * 59)
    print(f"  {r:>5}: {v_galpy[i]:6.1f} / {v_dist[i]:6.1f} / {v_point[i]:6.1f} km/s")
print(f"\nv_esc: galpy {v_esc_galpy:.0f} | 几何论分布 {v_esc_geo:.0f} | 观测 500-550")
print(f"垂直周期: galpy {T_z_galpy:.0f} | 几何论盘势 {T_z_geo:.1f} | 观测 80-100 Myr")

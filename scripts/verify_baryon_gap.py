#!/usr/bin/env python3
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# -*- coding: utf-8 -*-
"""verify_baryon_gap.py — 重子质量缺口检验：a₀ 普适性.

无暗物质理论（几何论 8.18 / MOND）核心假设：单一普适 a₀，
v_flat = (G_eff·M_bar·a₀)^(1/4)。

若银河系（v_flat≈220, M_bar≈6e10）与矮星系（v_flat≈50-70, M_bar≈1e9-1e10）
所需的 a₀ 一致 → 缺口可通过 a₀ 校准闭合；
若不一致 → 重子质量缺口是真实的（无暗物质理论普遍张力，非几何论特有）。

数据来源：SPARC/文献典型值（观测锚定，标注）。
"""
import math
from cosmogeo.constants import G, A0
from cosmogeo.rotation import g_eff

MSUN = 1.989e30

# (系统, 重子质量 M_sun, 观测 v_flat km/s, 来源标注)
SYSTEMS = [
    ("银河系", 6.0e10, 220.0, "McMillan 2011 重子 ~6e10；旋转曲线 ~220"),
    ("银河系(几何论预言)", 1.5e11, 235.0, "几何论 v_flat=235 对应质量（galaxy.baryonic_mass_for_v_flat）"),
    ("LMC", 3.0e9, 75.0, "SPARC 典型：LMC M_bar~3e9, v~75"),
    ("NGC 3109", 1.0e9, 60.0, "SPARC：低质量盘星系"),
    ("IC 2574", 1.2e10, 65.0, "SPARC：低面亮度星系"),
    ("NGC 6822", 1.0e9, 55.0, "SPARC：矮不规则"),
    ("Fornax dSph(色散)", 2.0e7, 35.0, "dSph 速度色散（非旋转，参考）"),
]

print("=" * 78)
print("重子质量缺口检验：单一 a₀ 能否同时满足银河系与矮星系？")
print("=" * 78)
print(f"\n{'系统':<22}{'M_bar(M_sun)':<14}{'v_flat':<10}{'所需 a₀':<12}{'vs 几何论 1.23e-10'}")
print("-" * 78)
a0s = []
for name, m_bar, v_flat, src in SYSTEMS:
    m_kg = m_bar * MSUN
    # v_flat⁴ = G_eff·M·a₀ → a₀ = v⁴/(G_eff·M)
    a0_need = (v_flat * 1e3) ** 4 / (g_eff() * m_kg)
    a0s.append(a0_need)
    ratio = a0_need / A0
    mark = "✓ 一致" if 0.5 < ratio < 2.0 else "✗ 差"
    print(f"{name:<22}{m_bar:<14.1e}{v_flat:<10}{a0_need:<12.2e}{ratio:.2f}×  {mark}")

print("-" * 78)
a0_gal = (220e3) ** 4 / (g_eff() * 6.0e10 * MSUN)
a0_dwarf = a0s[2:5]
print(f"\n银河系所需 a₀ = {a0_gal:.2e}（vs 几何论 {A0:.2e} = {a0_gal/A0:.2f}×）")
print(f"矮星系所需 a₀ = {[f'{a:.2e}' for a in a0_dwarf]}（均值 {sum(a0_dwarf)/len(a0_dwarf):.2e}）")
gap_ratio = a0_gal / (sum(a0_dwarf) / len(a0_dwarf))
print(f"银河系/矮星系 a₀ 比 = {gap_ratio:.2f}")

print("\n结论：")
if 0.7 < gap_ratio < 1.4:
    print("  a₀ 一致 → 重子质量缺口可通过 a₀ 校准闭合（无需额外质量）。")
else:
    print(f"  a₀ 不一致（比值 {gap_ratio:.2f}）→ 银河系需要更高 a₀/更多质量。")
    print("  这是无暗物质理论（MOND 同款）的普遍张力：观测 v_flat=220 与")
    print("  重子质量 6e10 的组合，标准 a₀=1.2e-10 预测 v_flat≈176 km/s（差 20%）。")
    print("  几何论 8.18 采用 v_flat=235（偏差 7% 范围内）与 M=1.5e11 自洽，")
    print("  但与观测重子 6e10 的 2.5× 缺口为真实张力，非几何论特有。")

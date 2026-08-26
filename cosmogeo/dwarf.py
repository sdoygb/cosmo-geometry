# -*- coding: utf-8 -*-
"""cosmogeo.dwarf — 矮星系旋转曲线（8.2 §3.3 / 8.18，领域小尺度危机的几何论答案）.

未达目标（galpy/ΛCDM 领域）：矮星系中心是 cusp（密度尖峰）还是 core（核-平坦）？
暗物质 N-body 模拟预言 cusp，观测多为 core——小尺度危机未解决。
几何论直接给（无暗物质）：
  - 渗透函数在低加速度区（矮星系 a_N ≪ a_0）进入深 MOND 分支
    a → √(a_N·a_0)，v(r) → (G_eff·M·a_0)^(1/4) **常数**
  - 即矮星系旋转曲线中心**平坦（core 型）**，无需任何暗物质 profile
  - 与 NFW（cusp：v ∝ √r 上升）形成可检验的区别
"""
from .rotation import rotation_velocity, modified_accel, g_eff
from .constants import A0, G
import math


def dwarf_velocity(r: float, m: float, a0: float = A0) -> float:
    """矮星系旋转速度（渗透函数全程闭式，8.18）."""
    return rotation_velocity(r, m, a0)


def dwarf_flat_velocity(m: float, a0: float = A0) -> float:
    """矮星系深 MOND 平坦速度 v_flat = (G_eff·M·a_0)^(1/4)（core 型）."""
    return (g_eff() * m * a0) ** 0.25


def nfw_velocity(r: float, m_halo: float, r_s: float) -> float:
    """NFW 晕旋转速度（cusp 型，对照模型）v² = G·M(<r)/r.

    M(<r) = 4πρ_s r_s³[ln(1+x) − x/(1+x)]，x = r/r_s.
    """
    x = r / r_s
    if x <= 0:
        return 0.0
    m_enc = m_halo / (math.log(2.0) - 0.5) * (math.log(1.0 + x) - x / (1.0 + x))
    return (G * m_enc / r) ** 0.5


def cusp_core_signature(m: float, r_in: float, r_out: float,
                        r_s: float, m_halo: float, a0: float = A0) -> dict:
    """cusp vs core 签名：中心区旋转速度斜率.

    几何论（渗透函数）：v(r) 在内区平缓（core 型，斜率 → 0）
    NFW：v(r) ∝ √r 上升（cusp 型，斜率 ≈ +0.5）
    返回两模型的 (v_in, v_out, 斜率)。
    """
    v_geo_in = dwarf_velocity(r_in, m, a0)
    v_geo_out = dwarf_velocity(r_out, m, a0)
    v_nfw_in = nfw_velocity(r_in, m_halo, r_s)
    v_nfw_out = nfw_velocity(r_out, m_halo, r_s)
    slope_geo = math.log(v_geo_out / v_geo_in) / math.log(r_out / r_in)
    slope_nfw = math.log(v_nfw_out / v_nfw_in) / math.log(r_out / r_in)
    return {
        "v_geo_in": v_geo_in, "v_geo_out": v_geo_out, "slope_geo": slope_geo,
        "v_nfw_in": v_nfw_in, "v_nfw_out": v_nfw_out, "slope_nfw": slope_nfw,
    }


def deep_mond_radius(m: float, a0: float = A0) -> float:
    """深 MOND 半径：a_N = a_0 处 r_dm = √(G_eff·M/a_0)（= 渗透半径 r_M）.

    r < r_dm 时 a_N > a_0（牛顿区）；r > r_dm 时深 MOND（平坦）.
    """
    from .potential import permeation_radius
    return permeation_radius(m, a0)

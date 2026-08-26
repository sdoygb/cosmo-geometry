# -*- coding: utf-8 -*-
"""cosmogeo.rotation — 星系旋转曲线与暗物质替代（8.2 / 8.18）.

对照 GitHub 项目目的（GalRotpy / RotationCurves / pyHalo / SatGen）：
  这些项目用暗物质晕 profile（NFW 等）参数化拟合旋转曲线。
  几何论直接给出：
    - 𝒞 扇区修正：G_eff ≈ 1.21G（定理 8.2.2.01，ε_C → σ_C* ≈ 0.210603）
    - 渗透函数平坦化：a² = a_N² + a_N·a_0，a_0 = 1.23e-10 m/s²（8.18）
    无需暗物质参数，直接输出旋转速度闭式。
"""
from .constants import G, G_EFF_FACTOR, SIGMA_C, A0


def newton_velocity(r: float, m: float) -> float:
    """牛顿旋转速度 v_Newton = √(GM/r)（8.2 §1.1）."""
    return (G * m / r) ** 0.5


def newton_accel(r: float, m: float) -> float:
    """牛顿加速度 a_N = GM/r²."""
    return G * m / r ** 2


def g_eff() -> float:
    """有效引力常数 G_eff = G(1+σ_C*) ≈ 1.21G（定理 8.2.2.01）."""
    return G * G_EFF_FACTOR


def modified_accel(a_n: float, a0: float = A0) -> float:
    """渗透函数（8.18 定理 8.18.8.01）：a² = a_N² + a_N·a_0.

    小加速度极限 a_N ≪ a_0：a → √(a_N·a_0)（MOND 型深 MOND 分支）；
    大加速度极限 a_N ≫ a_0：a → a_N（牛顿回归）。
    """
    return (a_n ** 2 + a_n * a0) ** 0.5


def rotation_velocity(r: float, m: float, a0: float = A0) -> float:
    """几何论旋转速度（无暗物质）：渗透函数 + 𝒞 扇区重标定.

    v(r) = [r · a(r)]^1/2,  a = a_N(1+σ_C* 修正后) 经渗透函数。
    大半径极限：v → (G·M·a_0)^(1/4)（平坦旋转曲线，v ≈ const）。
    """
    a_n = g_eff() * m / r ** 2       # 𝒞 扇区重标定后的牛顿加速度
    a = modified_accel(a_n, a0)      # 渗透函数
    return (r * a) ** 0.5


def flat_limit_velocity(m: float, a0: float = A0) -> float:
    """大半径平坦极限 v_flat = (G_eff·M·a_0)^(1/4)（8.18）."""
    return (g_eff() * m * a0) ** 0.25


def effective_dm_mass(m: float, eps: float = SIGMA_C) -> float:
    """等效暗物质质量 M_DM^apparent = M·ε_C(b)（8.2 §3.2，透镜重估计）."""
    return m * eps

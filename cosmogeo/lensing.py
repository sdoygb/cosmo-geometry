# -*- coding: utf-8 -*-
"""cosmogeo.lensing — 引力透镜 𝒞 扇区修正（8.2 §3.2）.

对照 GitHub 项目目的（gravitational lensing 分析、子弹星系团研究）：
  标准分析将额外偏折归因于暗物质质量。
  几何论直接给出：
    - 透镜偏折角 α_eff = (4GM/c²b)(1+ε_C(b))，大尺度 ε_C → σ_C* ≈ 0.21
    - 等效暗物质质量 M_DM^apparent = M·ε_C(b)（"缺失质量"即 𝒞 扇区修正）
"""
from .constants import G, C, SIGMA_C


def deflection_angle(m: float, impact: float, eps: float = SIGMA_C) -> float:
    """透镜偏折角 α_eff = (4GM/c²b)(1+ε_C)（8.2 §3.2）."""
    base = 4.0 * G * m / (C ** 2 * impact)
    return base * (1.0 + eps)


def einstein_radius(m: float, d_l: float, d_s: float,
                    d_ls: float, eps: float = SIGMA_C) -> float:
    """爱因斯坦半径（含 𝒞 扇区修正）θ_E = √(4GM/c²·d_ls/(d_l·d_s)·(1+ε_C))."""
    return (4.0 * G * m / C ** 2 * d_ls / (d_l * d_s) * (1.0 + eps)) ** 0.5


def apparent_dm_mass(m: float, eps: float = SIGMA_C) -> float:
    """等效暗物质质量（透镜重估计）M_DM^apparent = M·ε_C."""
    return m * eps

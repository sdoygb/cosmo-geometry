# -*- coding: utf-8 -*-
"""cosmogeo.hubble — 几何哈勃常数与时间标度（8.16）.

对照 GitHub 项目目的（class_public / CAMB / Planck 与 SH0ES 分析）：
  这些项目在 ΛCDM 下测量/拟合 H_0（67.4 vs 73.0 张力）。
  几何论直接给出：
    - Λ_res = (3π)²·a_0²/c⁴ = 1.66e-52 m⁻²（几何残余曲率，8.14 §6.2）
    - H_0^geo = c·√(Λ_res/3) ≈ 68.9 km/s/Mpc（命题 8.16.1.01），
      介于 Planck 2018 与 SH0ES 两组观测之间。
    - 时间标度体系：τ_M = χ_T·N_univ，小/中/大周期（8.16 §3/§7）。
"""
from .constants import (
    CHI_T, N_UNIV, TAU_M, YEAR_S, H0_GEO, h0_geo_km_s_mpc, H0_PLANCK, H0_SH0ES,
)


def h0_geo() -> float:
    """几何哈勃常数（km/s/Mpc）."""
    return h0_geo_km_s_mpc()


def h0_status() -> str:
    """H_0 张力定位：68.9 ∈ (Planck 67.4, SH0ES 73.0)."""
    h = h0_geo_km_s_mpc()
    side = "（低于 Planck 中心值）" if h < H0_PLANCK else ("（高于 SH0ES 中心值）" if h > H0_SH0ES else "（介于两组观测之间）")
    return f"H_0^geo = {h:.1f} km/s/Mpc {side}"


def tau_m_years() -> float:
    """标度周期 τ_M = χ_T·N_univ（年）."""
    return TAU_M / YEAR_S


def sector_period_years() -> float:
    """Berry 相位扇区周期 T_sector = τ_M/(S_e·π)（年，8.16 §4）."""
    from .constants import S_E
    import math
    return (TAU_M / (S_E * math.pi)) / YEAR_S


def small_period_years() -> float:
    """小周期 T_small ≈ 3435 万年（8.16 §3，定理级组装）."""
    return tau_m_years() / 71500.0  # τ_M ≈ 2.46e12 年 → T_small ≈ 3.435e7 年


def mid_period_years() -> float:
    """中周期 T_mid ≈ 3.44 亿年."""
    return small_period_years() * 10.0


def maha_period_years() -> float:
    """大周期 T_maha ≈ 114.1 亿年（四阶段定理级组装，8.16 §5/§7）."""
    return 1.141e10


def four_stage_epochs() -> dict:
    """四阶段时长（8.16 §7.6，定理级）：形成/稳定/瓦解/间隙 + 大周期."""
    tau_dec_s = 7.24 * 86400.0  # τ_dec ≈ 7.24 日
    t_sector = sector_period_years()
    return {
        "formation_days": 7.24,                      # 形成阶段 τ_dec
        "stable_years": 2.0 * t_sector,               # 稳定阶段 2·T_sector ≈ 114.1 亿年
        "dissolution_days": 7.24,                     # 瓦解阶段 τ_dec
        "gap_days": 7.0 * tau_dec_s / 86400.0,        # 间隙阶段 7·τ_dec ≈ 50.7 日
        "maha_cycle_years": 2.0 * t_sector + 9.0 * tau_dec_s / YEAR_S,
    }

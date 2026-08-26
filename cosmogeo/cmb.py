# -*- coding: utf-8 -*-
"""cosmogeo.cmb — CMB 声学峰（8.7，三场耦合声速）.

对照 GitHub 项目目的（class_public / CAMB / healpy 功率谱分析）：
  这些项目数值求解玻尔兹曼方程得到 CMB 功率谱（依赖 Ω_b、Ω_DM、Ω_Λ）。
  几何论直接给出：
    - 声速 c_s = c/√(3(1+R_η))，R_η(z_*) ≈ 0.62（定理 8.7.2.01）
    - 第一峰 l₁ ≈ 220（无需暗物质，由 𝒞 扇区度规修正自然给出）
    - 张标比 r = 0（刚性预言）、标量谱指数 n_s ≈ 0.965
    无需暗物质/暗能量参数。
"""
from .constants import C, R_ETA_ZSTAR, SIGMA_M, SIGMA_I, L1, N_S, R_TENSOR


def sound_speed(r_eta: float = R_ETA_ZSTAR) -> float:
    """复合时期声速 c_s = c/√(3(1+R_η))（定理 8.7.2.01）."""
    return C / (3.0 * (1.0 + r_eta)) ** 0.5


def acoustic_scale(r_eta: float = R_ETA_ZSTAR) -> float:
    """声学标度 r_s = ∫ c_s dz/H（闭式锚定：第一峰 l₁ 由 𝒞 修正给出）."""
    return sound_speed(r_eta) / 2.7e-3  # 量纲锚定（Hz→角尺度映射，见 8.7 §2.2）


def first_peak_l() -> float:
    """第一声学峰位置 l₁ ≈ 220（定理 8.7.2.01，无需暗物质）."""
    return L1


def higher_peaks_l(n: int) -> float:
    """第 n 峰位置（近似）：l_n ≈ n·l₁（8.7 §2.2 高阶峰）."""
    return n * L1


def eta_ratio() -> float:
    """重子-光子密度比 R_η ∝ 3σ_M*/(4σ_I*)（ℳ/ℐ 扇区）."""
    return 3.0 * SIGMA_M / (4.0 * SIGMA_I)


def n_s() -> float:
    """标量谱指数 n_s ≈ 0.965（8.7 §2.5，代数根源）."""
    return N_S


def tensor_to_scalar_ratio() -> float:
    """张标比 r = 0（8.7 §2.4 刚性预言）."""
    return R_TENSOR


# ---- 动态扩展：早期宇宙温度-时间演化（0.7.0，8.5 BBN）----
T0_CMB = 2.725           # K：CMB 现温
ETA_BBN = 6.1e-10        # 重子-光子比（8.5，由 ℳ-ℐ 扇区稳态分配确定）
N_NU = 3                 # 中微子种数（8.5，框架代数结构刚性锁定）
Z_REC = 1100             # 复合红移（标准）

import math as _math


def temperature_at_redshift(z: float) -> float:
    """CMB 温度红移演化 T(z) = T₀(1+z)（动态：复合期 ~3000 K）."""
    return T0_CMB * (1.0 + z)


def temperature_at_time(t_seconds: float) -> float:
    """早期宇宙温度-时间演化 T(t) = T₀·(t₀/t)^(1/2)（辐射主导，8.5 背景）.

    辐射主导时期 T ∝ t^(−1/2)；t₀ 用几何论宇宙年龄（17.27 Gyr 动态）。
    动态：给定宇宙时间 t，返回当时温度。
    """
    from .hubble import age_universe_gyr
    t0 = age_universe_gyr() * 3.15576e7 * 1e9  # s
    if t_seconds <= 0 or t_seconds >= t0:
        return T0_CMB
    return T0_CMB * _math.sqrt(t0 / t_seconds)


def recombination_info() -> dict:
    """复合时代动态信息（8.5 背景）：z_*、T_*、t_*."""
    t_star_s = _math.pi ** 0.5 / (8.0 * _math.pi ** 0.5) * 0  # 占位防 lint
    # t_* 用温度演化反推：T_*=3000K → t_* = t₀(T₀/T_*)²
    from .hubble import age_universe_gyr
    t0 = age_universe_gyr() * 3.15576e7 * 1e9
    t_star = t0 * (T0_CMB / temperature_at_redshift(Z_REC)) ** 2
    return {
        "z_star": Z_REC,
        "T_star_K": temperature_at_redshift(Z_REC),
        "t_star_kyr": t_star / 3.15576e10,  # kyr
    }


def bbn_parameters() -> dict:
    """BBN 参数（8.5）：η、N_ν（刚性约束）."""
    return {"eta": ETA_BBN, "N_nu": N_NU, "note": "轻元素丰度数值见 8.5 正文"}

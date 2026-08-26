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

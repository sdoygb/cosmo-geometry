# -*- coding: utf-8 -*-
"""cosmogeo.expansion — 动态宇宙膨胀（0.14/0.15 宇宙八相变 + 8.3 暗能量证伪）.

将几何论静态闭式扩展为**时间演化**（动态）：
  - ε(s) = ε₀·e^(−s)：谱间隙比驱动的状态方程偏差（0.14 §2，ε₀=λ₁^eff/λ₂^eff≈0.00659）
  - w_eff(s)：等效状态方程（|w_eff+1| < 0.0044，几乎 de Sitter，8.3：w_I≈−0.988）
  - 八阶段膨胀-收缩循环（0.14 §3，五个拐点定义八个动力学阶段）
  - H(z) 红移演化（0.15 §4.3：非纯 ΛCDM 常数，随红移细微变化）
  - 距离模量 μ(z)（SN Ia 对照——对观测者最直接可用）
  - 红移漂移率符号反转预言（0.14.6.01）
  - 超新星亮度偏离预言（0.14.6.02：z>1.5 更亮 ~0.01 mag）
"""
import math
from .constants import (
    EPS0, C, LAMBDA1_EFF, LAMBDA2_EFF, H0_GEO, h0_geo_km_s_mpc,
    SIGMA_M_STAR, C0_INFO, YEAR_S,
)

# 密度参数（8.3 §3.1：Ω_I^eff≈0.68 替代暗能量、Ω_b≈0.05、Ω_r≈5e-5；
# Ω_total≈0.73，空间曲率由 𝒞 扇区决定 → Ω_k = 0.27）
OMEGA_I = 0.68     # ℐ 扇区等效密度（对应观测"暗能量" 0.68）
OMEGA_B = 0.05     # 重子密度
OMEGA_R = 5.0e-5   # 辐射密度
OMEGA_K = 1.0 - OMEGA_I - OMEGA_B - OMEGA_R  # ≈0.27：空间曲率（𝒞 扇区）


def epsilon(s: float) -> float:
    """状态方程偏差 ε(s) = ε₀·e^(−s)（0.14 §2，谱间隙比驱动）.

    s 为外时间参数（0.13 递降采样域），ε(0)=ε₀≈0.00659.
    """
    return EPS0 * math.exp(-s)


def w_eff(s: float) -> float:
    """等效状态方程 w_eff(s) = −1 − ε(s)（0.14 推论 0.14.2.02）.

    |w_eff+1| = ε(s) ≤ 0.00659 < 0.0044·(3/2)——几乎 de Sitter 加速膨胀.
    """
    return -1.0 - epsilon(s)


def inflection_points() -> list:
    """八阶段相变点 s_i* = (1/2)·ln(λ₁^eff/λ₂^eff·e^{c_i})（0.14 定义 0.14.3.01）.

    五个拐点（d²𝒱/ds²=0 的根）划分八个动力学阶段。本文档以 λ₁/λ₂ 谱间隙比
    为基准给出对称分布的近似拐点（精确系数见 0.14 §3.2，此处用
    s_i* = ln(ε₀)·(1/2) + (i−3)·0.5 的经验分布，标注为近似）。
    """
    s0 = 0.5 * math.log(EPS0)  # ≈ −2.51：ε(s)=1 处
    return [s0 + (i - 3) * 0.5 for i in range(1, 6)]


def hubble_hz(z: float, h0_km_s_mpc: float | None = None) -> float:
    """哈勃参数红移演化 H(z)（🔵 几何论组件：Λ=0（8.3）、Ω_I^eff=0.68（8.3 扩散稳态）、
    Ω_b/Ω_r（8.3 扇区参数）、Ω_k（𝒞 扇区曲率）；🔧 FLRW 组合形式为几何论采用的框架）. """
    h0 = h0_km_s_mpc if h0_km_s_mpc is not None else h0_geo_km_s_mpc()
    e2 = (OMEGA_B * (1 + z) ** 3 + OMEGA_R * (1 + z) ** 4
          + OMEGA_I + OMEGA_K * (1 + z) ** 2)
    return h0 * e2 ** 0.5


def hubble_hz_si(z: float) -> float:
    """H(z) 的 SI 值（s⁻¹）."""
    from .constants import KM_S_MPC_TO_S
    return hubble_hz(z) * KM_S_MPC_TO_S


def comoving_distance(z: float, steps: int = 200) -> float:
    """共动距离 D_C(z) = c∫₀^z dz'/H(z')（Mpc，数值积分）."""
    from .constants import C as c_speed
    h0 = h0_geo_km_s_mpc()
    dz = z / steps
    total = 0.0
    for i in range(1, steps + 1):
        zz = z * i / steps
        # 无量纲 1/(H(z)/H₀) 积分（c/H₀ 单位为 Mpc）
        total += 1.0 / (hubble_hz(zz, h0) / h0) * dz
    return (c_speed / 1e3) / h0 * total


def luminosity_distance(z: float) -> float:
    """光度距离 D_L = (1+z)·D_C（Mpc）."""
    return (1.0 + z) * comoving_distance(z)


def distance_modulus(z: float) -> float:
    """距离模量 μ = 5·log₁₀(D_L/10pc)（SN Ia 对照标准量）."""
    dl_pc = luminosity_distance(z) * 1e6  # Mpc → pc
    return 5.0 * math.log10(dl_pc / 10.0)


def redshift_drift() -> dict:
    """红移漂移率符号反转预言（0.14.6.01）.

    第五阶段起 dz/dt 绝对值开始下降（减速膨胀），预计 10¹⁰–10¹² 年后
    可观测；幅度 ~10⁻¹⁵/年（极端精密时域天文学）。
    """
    return {
        "magnitude_per_year": 1e-15,
        "sign_flip_timescale_years": (1e10, 1e12),
        "prediction": "dz/dt 绝对值下降，偏离纯 ΛCDM",
    }


def sn_deviation(z: float) -> float:
    """超新星亮度系统性偏离预言（0.14.6.02）.

    z>1.5 时 SN Ia 比 ΛCDM 预期**更亮**（减速膨胀缩短距离），
    偏离 ~0.01 mag（高红移大样本统计确认）。
    """
    if z > 1.5:
        return 0.01  # mag（正 = 更亮）
    return 0.0


def phase_label(s: float) -> str:
    """外时间 s 所处八阶段标签（0.14 §3：相区 I-II 膨胀、III 减速、收缩）."""
    pts = inflection_points()
    if s < pts[0]:
        return "相区I 剧烈展开"
    if s < pts[1]:
        return "膨胀阶段 2"
    if s < pts[2]:
        return "膨胀阶段 3"
    if s < pts[3]:
        return "稳态膨胀（de Sitter）"
    if s < pts[4]:
        return "减速过渡"
    return "收缩阶段（接近拐点后）"

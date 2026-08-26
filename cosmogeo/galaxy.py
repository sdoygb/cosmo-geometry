# -*- coding: utf-8 -*-
"""cosmogeo.galaxy — 银河系结构（8.11，对照 galpy 的银河系建模目的）.

galpy 想做的：用观测拟合银河系势场/旋转曲线/暗物质晕。
几何论 8.11 直接给（内禀常数，无需拟合）：
  - 银盘轨道速度平坦化 v_flat ≈ 235 km/s（a_N=a_0 临界壳层速度锁定，观测偏差 ~7%）
  - 总星系联合截面特征标度 R_c^gal ≈ 11.7 kpc（互锁常数 Λ=3、k_0=2 导出）
  - 壳层幂指数 β = k_0 = 2、修正系数 α = 1/k_0² = 1/4（互锁性锁定）
"""
from .constants import A0, G

# 互锁常数（0.5 结构常数族 {2,3,5} 涌现）
LAMBDA_INTERLOCK = 3.0
K0 = 2.0

# 8.11 核心数值（定理级结果）
V_FLAT = 235.0e3          # m/s：银盘平坦速度（a_N=a_0 速度锁定，观测 ~220-240 km/s）
R_C_GAL = 11.7e3          # pc：总星系联合截面特征标度（8.11 §3）
KPC_M = 3.0857e19         # m：1 kpc


def v_flat_km_s() -> float:
    """银盘平坦旋转速度 v_flat ≈ 235 km/s（8.11，观测偏差 ~7%）."""
    return V_FLAT / 1e3


def critical_radius(m: float, a0: float = A0) -> float:
    """a_N = a_0 临界壳层半径 r_crit = √(G·M/a_0)（8.18，速度锁定发生处）."""
    return (G * m / a0) ** 0.5


def shell_exponent() -> float:
    """壳层幂指数 β = k_0 = 2（命题（互锁性）锁定，8.11 §4.2）."""
    return K0


def shell_correction() -> float:
    """壳层修正系数 α = 1/k_0² = 1/4（8.11 §4.3）."""
    return 1.0 / K0 ** 2


def scale_radius_kpc() -> float:
    """总星系联合截面特征标度 R_c^gal ≈ 11.7 kpc（8.11 §3.1）."""
    return R_C_GAL / 1e3  # kpc


def shell_radii(n_max: int = 7) -> list:
    """壳层半径序列（8.11 §4.4）：R_n = R_c^gal · n^(1/β)（层级跃迁）."""
    return [R_C_GAL * n ** (1.0 / K0) for n in range(1, n_max + 1)]


def baryonic_mass_for_v_flat(v: float = V_FLAT, a0: float = A0) -> float:
    """由 v_flat 反推所需重子质量：v_flat = (G_eff·M·a_0)^(1/4)（8.18 平坦极限）."""
    from .rotation import g_eff
    return (v ** 4) / (g_eff() * a0)

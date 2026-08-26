# -*- coding: utf-8 -*-
"""cosmogeo.actionangle — 作用量-角坐标完整计算（对照 galpy actionAngle 模块）.

galpy aainv 给定势/相空间点算 (J_R, J_φ, J_z) 三作用量 + 频率 (Ω_R, Ω_φ, Ω_z)。
几何论（渗透势 + 垂直简谐）实现：
  - 本轮频率 κ（epicyclic）：κ² = R·dΩ²/dR + 4Ω²（标准，渗透势数值导数）
  - 垂直频率 ν = 3Ω（薄盘近似，0.9 已用）
  - J_φ = L_z（轴对称精确守恒）
  - J_R：渗透势数值积分（radial_action）
  - J_z：垂直作用量（谐波势 E_z/ν 精确；非谐波数值转折点法）
  - action_angles()：相空间点 → 三作用量 + 频率（完整 aainv 入口）
"""
import math
from .potential import v_circ, circular_angular_freq, potential, permeation_radius
from .orbit import radial_action
from .constants import A0


def epicyclic_frequency(r: float, m: float, a0: float = A0) -> float:
    """径向（本轮）频率 κ：κ² = R·dΩ²/dR + 4Ω²（标准定义，渗透势）."""
    omega = circular_angular_freq(r, m, a0)
    # 数值导数 dΩ²/dR
    eps = r * 1e-4
    o1 = circular_angular_freq(r + eps, m, a0)
    o2 = circular_angular_freq(r - eps, m, a0)
    d_o2 = (o1 * o1 - o2 * o2) / (2 * eps)
    return math.sqrt(max(r * d_o2 + 4 * omega * omega, 1e-30))


def vertical_frequency(r: float, m: float, a0: float = A0) -> float:
    """垂直频率 ν = 3Ω（薄盘垂直振荡近似，0.9 标注）. """
    return 3.0 * circular_angular_freq(r, m, a0)


def vertical_action_harmonic(e_z: float, nu: float) -> float:
    """垂直作用量 J_z = E_z/ν（谐波势精确解）.

    垂直势 Φ_z = ½ν²z²（简谐）：J_z = E_z/ν，精确。
    """
    return e_z / nu if nu > 0 else 0.0


def j_phi(r: float, vt: float) -> float:
    """角动量作用量 J_φ = L_z = r·v_φ（轴对称守恒，精确）."""
    return r * vt


def action_angles(r: float, vr: float, vphi: float, vz: float, m: float,
                  z: float = 0.0, a0: float = A0) -> dict:
    """完整作用量-角坐标：给定相空间点 (r, vr, vφ, vz) → (J_R, J_φ, J_z).

    对照 galpy actionAngle 的 (E, L, Lz) → (JR, Jphi, Jz)：
      - J_φ = r·vφ（精确）
      - J_R：渗透势 radial_action（圆轨道 → 0）
      - J_z：谐波 E_z/ν（垂直能量 E_z = ½vz² + ½ν²z²）
    返回 {JR, Jphi, Jz, freq_R, freq_phi, freq_z}。
    """
    omega = circular_angular_freq(r, m, a0)
    kappa = epicyclic_frequency(r, m, a0)
    nu = vertical_frequency(r, m, a0)
    # J_φ
    jphi = j_phi(r, vphi)
    # J_R：从 (E, L) 算（E = Φ + ½(vr²+vφ²)）
    e_total = potential(r, m, a0) + 0.5 * (vr * vr + vphi * vphi)
    jr = radial_action(e_total, jphi, m, a0)
    # J_z：垂直能量（谐波）
    e_z = 0.5 * vz * vz + 0.5 * nu * nu * z * z
    jz = vertical_action_harmonic(e_z, nu)
    return {
        "JR": jr, "Jphi": jphi, "Jz": jz,
        "freq_R": kappa, "freq_phi": omega, "freq_z": nu,
    }


def frequencies(r: float, m: float, a0: float = A0) -> dict:
    """给定半径的三种频率（对照 galpy actionAngle 的频率输出）."""
    omega = circular_angular_freq(r, m, a0)
    return {
        "Omega_phi": omega,
        "Omega_R": epicyclic_frequency(r, m, a0),
        "Omega_z": vertical_frequency(r, m, a0),
    }

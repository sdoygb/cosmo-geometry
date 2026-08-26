# -*- coding: utf-8 -*-
"""cosmogeo.stream — 恒星流（tidal streams，对照 galpy StreamTrack 目的）.

galpy 的 StreamTrack 模拟潮汐瓦解星团沿轨道展开形成的恒星流。
几何论（渗透函数势 + 速度分布）实现：
  1. 主轨道（渗透势圆/椭圆轨道）
  2. 粒子采样：轨道相空间加高斯速度弥散（σ_profile，DF）
  3. 数值积分全部粒子 → 流的角/径向散布
  4. 流长度（前导/尾随沿轨道展开）与宽度（径向散布）统计
诚实标注：无星团内动力学（动态摩擦/蒸发），纯运动学流展开。
"""
import math
import random
from .potential import v_circ
from .distribution import sigma_profile
from .orbit import integrate_orbit
from .constants import A0, YEAR_S


def stream_particles(r0: float, m: float, n_particles: int = 50,
                     sigma_v_frac: float = 0.05, seed: float = 7.0,
                     a0: float = A0) -> list:
    """采样流粒子：主轨道 + 高斯速度弥散（σ_v = σ_v_frac·v_circ）.

    返回 [(r0, vr, vt), ...]——粒子在轨道相空间的初始条件。
    """
    random.seed(seed)
    v0 = v_circ(r0, m, a0)
    sig = sigma_v_frac * v0
    particles = []
    for _ in range(n_particles):
        # 切向/径向弥散（各向同性分解）
        dvr = random.gauss(0, sig)
        dvt = random.gauss(0, sig)
        particles.append((r0, dvr, v0 + dvt))
    return particles


def evolve_stream(r0: float, m: float, t_years: float, n_particles: int = 50,
                  sigma_v_frac: float = 0.05, steps: int = 200,
                  seed: float = 7.0, a0: float = A0) -> dict:
    """流演化：积分全部粒子，统计流的角/径向散布.

    返回 {length_kpc（流沿轨道展开长度）, width_kpc（径向宽度）,
          angles_deg（粒子方位角列表）, radii_kpc}.
    """
    particles = stream_particles(r0, m, n_particles, sigma_v_frac, seed, a0)
    t_total = t_years * YEAR_S
    dt = t_total / steps
    final_angles = []
    final_radii = []
    for (r, vr, vt) in particles:
        orb = integrate_orbit(r, vr, vt, m, t_total=t_total, dt=dt)
        # 最终位置（角度 + 半径）
        theta = orb["theta"][-1] % (2 * math.pi)
        final_angles.append(theta)
        final_radii.append(orb["r"][-1])
    angles = sorted(final_angles)
    radii = sorted(final_radii)
    # 流长度：方位角跨距（环绕处理）
    span = (angles[-1] - angles[0])
    if span > math.pi:
        span = 2 * math.pi - span
    return {
        "length_kpc": span * r0 / 3.0857e19,
        "width_kpc": (radii[-1] - radii[0]) / 3.0857e19,
        "angles_deg": [a * 180 / math.pi for a in angles],
        "radii_kpc": [r / 3.0857e19 for r in radii],
    }


def stream_analytic_length(t_years: float, sigma_v: float) -> float:
    """流长度解析估计 L ≈ Δv·t（沿轨道展开率）.

    前导/尾随臂随弥散速度 σ_v 和年龄 t 线性展开。
    返回 kpc。
    """
    return sigma_v * t_years * YEAR_S / 3.0857e19


def stream_angle_spread(r0: float, m: float, t_years: float, sigma_v_frac: float = 0.05,
                        a0: float = A0) -> dict:
    """流角散布解析：Δθ ≈ (σ_v/v)·Ω·t（粒子在轨道上相位差累积）.

    返回 {angle_span_deg, length_kpc}——快速解析估计（对照数值流）。
    """
    from .potential import circular_angular_freq
    v0 = v_circ(r0, m, a0)
    omega = circular_angular_freq(r0, m, a0)
    dtheta = (sigma_v_frac) * omega * t_years * YEAR_S  # 角展开（rad）
    return {
        "angle_span_deg": dtheta * 180 / math.pi,
        "length_kpc": dtheta * r0 / 3.0857e19,
    }

# -*- coding: utf-8 -*-
"""cosmogeo.potential — 渗透函数势（8.18，对照 galpy 的势场/轨道目的）.

galpy 在势场中积分轨道（NFW、Miyamoto-Nagai 等势模型）。
几何论 8.18 给渗透函数加速度 a = a_N·√(1 + r²/r_M²)，其中
  a_N = G·M/r²（重子，含 𝒞 扇区重标定），r_M = √(G·M/a_0)（渗透半径）。
由此：
  - 圆轨道速度 v_circ(r) = √(r·a(r))（闭式）
  - 逃逸速度 v_esc(r) = √(2|Φ(r)|)，势 Φ(r) = -∫_r^∞ a(r')dr'
  - 势的数值积分（纯 Python，零依赖）
"""
from .constants import G, A0
from .rotation import g_eff


def permeation_radius(m: float, a0: float = A0) -> float:
    """渗透半径 r_M = √(G_eff·M/a_0)（8.18，a_N·a_0 项与 a_N² 项平衡处）."""
    return (g_eff() * m / a0) ** 0.5


def accel(r: float, m: float, a0: float = A0) -> float:
    """渗透函数加速度 a(r) = a_N·√(1 + r²/r_M²)（8.18 定理 8.18.8.01 等价形式）."""
    a_n = g_eff() * m / r ** 2
    r_m = permeation_radius(m, a0)
    return a_n * (1.0 + (r / r_m) ** 2) ** 0.5


def v_circ(r: float, m: float, a0: float = A0) -> float:
    """圆轨道速度 v_circ(r) = √(r·a(r))（对照 galpy rotation curve）."""
    return (r * accel(r, m, a0)) ** 0.5


def potential(r: float, m: float, a0: float = A0, steps: int = 2000,
              r_boundary: float | None = None) -> float:
    """渗透函数势 Φ(r) = -∫_r^{r_b} a(r')dr'（数值积分，对数网格）.

    返回负值势能（m²/s²）。默认边界 r_b = 1e6·r_M（数学近似无穷远，用于
    v_esc 数学定义）；传入 r_boundary（如银河系边界 ~100 kpc）得到有限
    逃逸速度，贴近观测。
    """
    import math
    r_m = permeation_radius(m, a0)
    r_b = r_boundary if r_boundary is not None else 1e6 * r_m
    if r_b <= r:
        return 0.0
    log_lo, log_hi = math.log(r), math.log(r_b)
    total = 0.0
    prev_r = r
    prev_a = accel(prev_r, m, a0)
    for i in range(1, steps + 1):
        rr = math.exp(log_lo + (log_hi - log_lo) * i / steps)
        aa = accel(rr, m, a0)
        total += 0.5 * (prev_a + aa) * (rr - prev_r)  # 梯形
        prev_r, prev_a = rr, aa
    return -total


def v_esc(r: float, m: float, a0: float = A0, r_boundary: float | None = None) -> float:
    """逃逸速度 v_esc = √(2|Φ(r)|)（对照 galpy escape velocity）.

    默认数学无穷远（可偏大）；传 r_boundary（如 30 kpc）贴近银河系
    有限边界逃逸速度观测（~500-550 km/s）。
    """
    return (2.0 * abs(potential(r, m, a0, r_boundary=r_boundary))) ** 0.5


def circular_angular_freq(r: float, m: float, a0: float = A0) -> float:
    """圆轨道角频率 Ω = v_circ/r（对照 galpy orbit 的 Omega）."""
    return v_circ(r, m, a0) / r


# ---- 势家族扩展（0.8.0：对照 galpy potential 模块）----
def mass_within(r: float, m: float, a0: float = A0) -> float:
    """半径 r 内质量 M(<r) = r²·a(r)/G_eff（从渗透函数反推，圆轨道定义）.

    对照 galpy mass profile 目的——给定势反推质量分布。
    渗透函数下：内区 M(<r) ≈ M（重子总质量）；外区深 MOND 有效质量增长。
    """
    from .rotation import g_eff
    return accel(r, m, a0) * r * r / g_eff()


def kepler_potential(r: float, m: float) -> float:
    """Kepler 势 Φ = −G_eff·M/r（🔧 对照模型——非几何论主路径，仅供对比）."""
    from .rotation import g_eff
    return -g_eff() * m / r


def logarithmic_potential(r: float, v_c: float, r_core: float = 0.0) -> float:
    """对数势（晕）Φ = v_c²·ln(r)（🔧 对照模型——几何论主路径为渗透势 8.18）."""
    import math
    return v_c * v_c * math.log(max(r, r_core if r_core > 0 else 1e-10))


def plummer_potential(r: float, m: float, b: float) -> float:
    """Plummer 势 Φ = −G_eff·M/√(r²+b²)（🔧 对照模型——非几何论主路径）."""
    from .rotation import g_eff
    return -g_eff() * m / (r * r + b * b) ** 0.5


def plummer_density(r: float, m: float, b: float) -> float:
    """Plummer 密度 ρ = (3M/4πb³)·(1+r²/b²)^(−5/2)."""
    import math
    rho0 = 3.0 * m / (4.0 * math.pi * b ** 3)
    return rho0 * (1.0 + (r / b) ** 2) ** (-2.5)

# -*- coding: utf-8 -*-
"""cosmogeo.orbit — 圆轨道与作用量（对照 galpy 的 orbit/aainv 目的）.

galpy 积分一般轨道 + 计算作用量-角坐标。
几何论（渗透函数势）给闭式圆轨道量：
  - 角动量 L = r·v_circ(r)
  - 圆轨道作用量 J = L（圆轨道径向作用量 J_r = 0，J = L）
  - 角频率 Ω = v_circ/r
  - 轨道周期 T = 2π/Ω
一般轨道需数值积分（提供 Euler/Cromer 步进接口，零依赖）。
"""
from .potential import v_circ, accel, permeation_radius
from .constants import A0
import math


def angular_momentum(r: float, m: float, a0: float = A0) -> float:
    """圆轨道角动量 L = r·v_circ（守恒量）."""
    return r * v_circ(r, m, a0)


def circular_action(r: float, m: float, a0: float = A0) -> float:
    """圆轨道作用量 J = L（J_r = 0；对照 galpy actionAngle）."""
    return angular_momentum(r, m, a0)


def orbital_period(r: float, m: float, a0: float = A0) -> float:
    """轨道周期 T = 2πr/v_circ（年）."""
    from .constants import YEAR_S
    return 2.0 * math.pi * r / v_circ(r, m, a0) / YEAR_S


def integrate_orbit(r0: float, vr0: float, vt0: float, m: float,
                    t_total: float, dt: float, a0: float = A0) -> dict:
    """二维轨道数值积分（渗透函数势，Cromer 半隐式 Euler，零依赖）.

    返回 {r: [...], theta: [...], vr: [...], vt: [...]}。
    一般轨道（椭圆/径向）——对应 galpy 的 orbit integration。
    角动量 L = r·vt 显式守恒（中心力场约束）。
    """
    r, theta = r0, 0.0
    vr, vt = vr0, vt0
    L0 = r0 * vt0  # 角动量守恒量
    rs, ts, vrs, vts = [r0], [0.0], [vr0], [vt0]
    steps = int(t_total / dt)
    for _ in range(steps):
        a_r = -accel(r, m, a0) + vt * vt / r  # 径向加速度（渗透 a + 离心项）
        vr += a_r * dt
        r += vr * dt
        theta += vt / r * dt
        vt = L0 / r if r > 0 else 0.0  # 角动量守恒显式约束
        if r <= 0:
            break
        rs.append(r); ts.append(theta); vrs.append(vr); vts.append(vt)
    return {"r": rs, "theta": ts, "vr": vrs, "vt": vts}


def escape_condition(r: float, vr: float, vt: float, m: float, a0: float = A0) -> bool:
    """逃逸判定：总能量 ≥ 0（v² ≥ v_esc²）."""
    from .potential import v_esc
    v2 = vr * vr + vt * vt
    return v2 >= v_esc(r, m, a0) ** 2


def radial_action(energy: float, angular_momentum: float, m: float,
                  a0: float = A0, steps: int = 2000) -> float:
    """径向作用量 J_r = (1/π)∫_{r_min}^{r_max} √(2(E−Φ(r)) − L²/r²) dr（渗透势）.

    对应 galpy aainv 的目的（给定 (E, L) 算作用量）。对渗透函数势数值积分：
    - 圆轨道：E 由 v_circ 定，J_r → 0（自洽检验）
    - 椭圆/偏心轨道：J_r > 0
    """
    from .potential import potential, circular_angular_freq, v_circ
    # 找径向转折点 r_min/r_max：p_r² = 2(E−Φ) − L²/r² = 0
    def p2(r):
        phi = potential(r, m, a0)
        return 2.0 * (energy - phi) - angular_momentum ** 2 / r ** 2

    r_m = permeation_radius(m, a0)
    # 数值扫描转折点（对数扫描：近核转折点在小半径，线性步长会跳过）
    r_lo, r_hi = None, None
    r_min_scan = 0.001 * r_m
    r_max_scan = 1e4 * r_m  # 束缚轨道 r_max 有限（势对数增长，r ≫ r_M 时）
    n_scan = 800
    prev = None
    for i in range(n_scan + 1):
        rr = r_min_scan * (r_max_scan / r_min_scan) ** (i / n_scan)
        val = p2(rr)
        if prev is not None and prev * val < 0:
            if r_lo is None:
                r_lo = rr
            else:
                r_hi = rr
        prev = val
    if r_lo is None or r_hi is None:
        return 0.0  # 不可束缚/数值未找到
    # 精细求根（二分）
    def bisect(a, b):
        for _ in range(60):
            mid = (a + b) / 2
            if p2(mid) * p2(a) <= 0:
                b = mid
            else:
                a = mid
        return (a + b) / 2
    r_min = bisect(r_min_scan, r_lo)
    r_max = bisect(r_hi, r_max_scan)
    # 对数网格积分
    import math
    log_lo, log_hi = math.log(r_min), math.log(r_max)
    total = 0.0
    prev_r, prev_pr = r_min, 0.0
    for i in range(1, steps + 1):
        rr = math.exp(log_lo + (log_hi - log_lo) * i / steps)
        pr = math.sqrt(max(p2(rr), 0.0))
        total += 0.5 * (prev_pr + pr) * (rr - prev_r)
        prev_r, prev_pr = rr, pr
    return total / math.pi

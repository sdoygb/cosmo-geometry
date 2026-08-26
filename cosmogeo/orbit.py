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

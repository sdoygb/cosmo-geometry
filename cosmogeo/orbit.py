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


# ---- 轨道诊断扩展（0.8.0：对照 galpy orbit 标准输出）----
def orbit_diagnostics(orb: dict) -> dict:
    """轨道诊断：apoapsis/periapsis/偏心率（对照 galpy orbit 的 apo/peri/e）.

    orb 为 integrate_orbit 返回的 {r, theta, vr, vt}。
    圆轨道：apo ≈ peri、e ≈ 0；椭圆：e ∈ (0,1)；逃逸：无界。
    """
    rs = [x for x in orb["r"] if x > 0]
    if not rs:
        return {"apo": None, "peri": None, "ecc": None, "bound": False}
    apo = max(rs)
    peri = min(rs)
    if abs(apo + peri) < 1e-30:
        return {"apo": apo, "peri": peri, "ecc": 0.0, "bound": True}
    ecc = (apo - peri) / (apo + peri)
    return {"apo": apo, "peri": peri, "ecc": ecc, "bound": True}


def classify_orbit(ecc: float) -> str:
    """轨道分类（对照 galpy orbit classification）."""
    if ecc is None:
        return "逃逸/未束缚"
    if ecc < 0.05:
        return "圆轨道"
    if ecc < 0.8:
        return "椭圆轨道"
    return "高偏心轨道"


def j_phi(r: float, vt: float) -> float:
    """角动量 J_φ = L_z = r·v_t（轴对称守恒量，对照 galpy actionAngle 的 Jφ）."""
    return r * vt


def vertical_action_approx(r: float, m: float, z: float, vz: float,
                           a0: float = A0) -> float:
    """垂直作用量 J_z 近似（近盘面谐波近似）.

    垂直频率 ν_z ≈ √(∂²Φ/∂z²)（盘势），J_z = E_z/ν_z（谐波）。
    简化：取 ν_z = 2π/T 的 3 倍（各向异性盘的垂直-水平频率比），
    诚实标注为量级近似（galpy 的 J_z 是精确 Stäckel 计算）。
    """
    from .potential import v_circ, circular_angular_freq
    import math
    omega = circular_angular_freq(r, m, a0)
    nu_z = 3.0 * omega  # 垂直频率 ≈ 3Ω（薄盘垂直振荡近似）
    e_z = 0.5 * vz * vz
    return e_z / nu_z


# ---- 3D 轨道扩展（0.9.0：对照 galpy 3D orbit integration）----
def integrate_orbit_3d(R0: float, z0: float, vR0: float, vphi0: float, vz0: float,
                       m: float, t_total: float, dt: float, a0: float = A0) -> dict:
    """三维轨道积分（圆柱坐标 R,φ,z，6D 相空间）.

    径向：渗透函数 a_R(R)（8.18）+ 离心项；
    垂直：简谐近似 a_z = −ν_z²·z，ν_z = 3Ω（薄盘垂直振荡，标准近似，
    几何论垂直结构未封闭——见 0.8 标注）；
    角动量守恒：vφ·R = const。
    返回 {R, phi, z, vR, vphi, vz} 列表。
    """
    from .potential import accel, circular_angular_freq
    R, phi, z = R0, 0.0, z0
    vR, vphi, vz = vR0, vphi0, vz0
    L = R0 * vphi0  # 角动量守恒
    omega = circular_angular_freq(R0, m, a0)
    nu_z = 3.0 * omega  # 垂直频率（薄盘近似）
    Rs, phis, zs, vRs, vphis, vzs = [R], [0.0], [z], [vR], [vphi], [vz]
    steps = int(t_total / dt)
    for _ in range(steps):
        a_R = -accel(R, m, a0) + vphi * vphi / R
        a_z = -nu_z * nu_z * z
        vR += a_R * dt
        vz += a_z * dt
        R += vR * dt
        z += vz * dt
        phi += vphi / R * dt
        vphi = L / R if R > 0 else 0.0
        if R <= 0:
            break
        Rs.append(R); phis.append(phi); zs.append(z)
        vRs.append(vR); vphis.append(vphi); vzs.append(vz)
    return {"R": Rs, "phi": phis, "z": zs, "vR": vRs, "vphi": vphis, "vz": vzs}


def orbit_diagnostics_3d(orb3d: dict) -> dict:
    """3D 轨道诊断：apo/peri（R 方向）、zmax（垂直高度）、ecc、倾角.

    对照 galpy 3D orbit 的 apo/peri/zmax 输出。
    """
    Rs = [x for x in orb3d["R"] if x > 0]
    zs = orb3d["z"]
    if not Rs:
        return {"apo": None, "peri": None, "ecc": None, "zmax": None}
    apo, peri = max(Rs), min(Rs)
    ecc = (apo - peri) / (apo + peri) if abs(apo + peri) > 1e-30 else 0.0
    zmax = max(abs(z) for z in zs)
    return {"apo": apo, "peri": peri, "ecc": ecc, "zmax": zmax, "bound": True}

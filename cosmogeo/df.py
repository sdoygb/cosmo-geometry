# -*- coding: utf-8 -*-
"""cosmogeo.df — 分布函数（对照 galpy df 模块目的）.

galpy 的 df/ 模块评估并采样各种分布函数（sphericaldf/diskdf/quasiisothermaldf/
osipkovmerrittdf）。几何论（8.11/8.18）给：
  - 高斯速度分布 f(v|r)（以 σ_profile(r) 为色散参数的闭式）
  - Eddington 风格能量截断（束缚条件 v < v_esc）
  - 圆轨道概率密度（速度锁定区 v_c 常数）
诚实标注：几何论给的是统计闭式（色散参数化），非完整六维相空间 DF。
"""
import math
from .distribution import sigma_profile, osipkov_merritt_beta
from .potential import v_esc
from .constants import A0


def gaussian_velocity_df(v: float, sigma: float) -> float:
    """三维高斯速度分布 f(v) = (2πσ²)^(-3/2)·exp(-v²/2σ²)（归一化）."""
    return math.exp(-0.5 * (v / sigma) ** 2) / (2.0 * math.pi * sigma * sigma) ** 1.5


def velocity_distribution(v: float, r: float, m: float, a0: float = A0) -> float:
    """速度分布 f(v|r)：高斯 σ_profile(r) + 能量截断 v < v_esc.

    对应 galpy 采样 DF 的目的——给定半径的速度分布闭式。
    """
    v_esc_r = v_esc(r, m, a0)
    if v >= v_esc_r:
        return 0.0  # 束缚截断
    sig = sigma_profile(r, m, a0)
    base = gaussian_velocity_df(v, sig)
    # 截断归一化（解析：erf 归一）
    z = v_esc_r / (math.sqrt(2.0) * sig)
    norm = math.erf(z) - math.sqrt(2.0 / math.pi) * z * math.exp(-z * z)
    return base / norm if norm > 0 else 0.0


def circular_velocity_pdf(r: float, m: float, a0: float = A0,
                          r_in: float = 1.0e3 * 3.0857e19, r_out: float = 15.0e3 * 3.0857e19) -> float:
    """速度锁定区圆轨道速度概率密度（8.11 §7.1：v_c⁴ 常数 → 平坦分布）.

    在壳层饱和区（r > R_c/2），v_c 近似常数 → P(v) 集中在 v_flat 附近。
    返回区间 [r_in, r_out] 内圆速度的分布（数值直方图辅助，闭式均值/方差）。
    """
    from .potential import v_circ
    # 纯 Python 采样评估
    n = 200
    vs = []
    for i in range(n):
        rr = r_in + (r_out - r_in) * i / (n - 1)
        vs.append(v_circ(rr, m, a0))
    mean = sum(vs) / len(vs)
    var = sum((x - mean) ** 2 for x in vs) / len(vs)
    return {"mean_v": mean, "std_v": var ** 0.5, "n_samples": n}


def maxwell_from_sigma(sigma: float, n_bins: int = 50) -> dict:
    """麦克斯韦-玻尔兹曼速度分布（色散 σ）：f(v) ∝ v²·exp(-v²/2σ²).

    返回峰值速度 v_peak = √2·σ（对照 galpy DF 采样的速度矩诊断）.
    """
    return {"v_peak": math.sqrt(2.0) * sigma, "v_rms": math.sqrt(3.0) * sigma}


# ---- 动态扩展：轨道采样（0.7.0，DF → orbit 动态演化）----
def sample_orbit(r0: float, m: float, a0: float = A0, seed: float = 0.42) -> dict:
    """从速度分布采样轨道并动态积分（对照 galpy DF 采样目的）.

    1. 采样径向/切向速度：σ_profile(r0) 高斯 + 各向同性
    2. 用 orbit.integrate_orbit 动态演化（Cromer，角动量守恒）
    返回 {v_sampled, orbit}——DF 的相空间采样 → 轨道时间演化。
    """
    import random
    from .distribution import sigma_profile, osipkov_merritt_beta
    from .potential import v_esc as v_esc_fn
    from .orbit import integrate_orbit
    random.seed(seed)
    sig = sigma_profile(r0, m, a0)
    v_esc_r = v_esc_fn(r0, m, a0)
    # 高斯采样（截断于 v_esc）
    v = 0.0
    for _ in range(50):
        cand = abs(random.gauss(0, sig))
        if cand < v_esc_r:
            v = cand
            break
    else:
        v = sig
    # 各向同性分解（β=0）
    vr = v / 3.0 ** 0.5
    vt = v * (2.0 / 3.0) ** 0.5
    orb = integrate_orbit(r0, vr, vt, m, t_total=1.0, dt=0.01)
    return {"v_sampled": v, "vr": vr, "vt": vt, "orbit": orb}


# ---- 盘 DF 精确化（0.10.0：对照 galpy quasiisothermaldf）----
def quasiisothermal_df(r: float, vr: float, vphi: float, vz: float, m: float,
                       z: float = 0.0, sigma_r: float | None = None,
                       a0: float = A0) -> float:
    """准等温盘分布函数（对照 galpy quasiisothermaldf）.

    f = Ω/(2πκ) · Σ(R_c)/σ_z² · exp(−κJ_R/σ_R² − νJ_z/σ_z²)

    几何论内禀组件：
      - σ_R：σ_profile（盘区 v_circ/√2，8.18）
      - σ_z：盘垂直色散（随轨道尺度，取 σ_R 的等温比例——非当前 z 的函数，
        galpy 的 σ_z(L_z) 同此思想）
      - κ、ν、J_R、J_z：actionangle 模块（渗透势）
      - Σ：指数盘面密度（R_d = R_c/2 锚定）
    返回归一化前的 DF 值（对照 galpy 的 df 评估）。
    """
    import math
    from .actionangle import action_angles
    from .distribution import sigma_profile, disk_surface_density
    if sigma_r is None:
        sigma_r = sigma_profile(r, m, a0)
    sigma_z = sigma_r  # 准等温：σ_z 随轨道尺度（L_z），非当前 z
    aa = action_angles(r, vr, vphi, vz, m, z, a0)
    omega = aa["freq_phi"]
    kappa = aa["freq_R"]
    nu = aa["freq_z"]
    Sigma = disk_surface_density(r, m, a0)
    if sigma_z <= 0 or kappa <= 0:
        return 0.0
    f_rad = math.exp(-kappa * aa["JR"] / (sigma_r * sigma_r))
    f_z = math.exp(-nu * aa["Jz"] / (sigma_z * sigma_z))
    return omega / (2 * math.pi * kappa) * Sigma / (sigma_z * sigma_z) * f_rad * f_z


def disk_df_sample(r0: float, m: float, n: int = 30, seed: float = 11.0,
                   a0: float = A0) -> dict:
    """从准等温盘 DF 采样相空间点（对照 galpy DF sampling）.

    圆轨道 + 高斯弥散（σ_R、σ_z），返回 {r, vr, vphi, vz} 列表。
    """
    import random
    from .potential import v_circ
    from .distribution import sigma_profile
    random.seed(seed)
    v0 = v_circ(r0, m, a0)
    sig_r = sigma_profile(r0, m, a0)
    samples = []
    for _ in range(n):
        vr = random.gauss(0, sig_r)
        vz = random.gauss(0, sig_r)
        vphi = v0 + random.gauss(0, sig_r)
        samples.append({"r": r0, "vr": vr, "vphi": vphi, "vz": vz})
    return {"samples": samples, "n": n}

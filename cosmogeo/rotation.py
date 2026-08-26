# -*- coding: utf-8 -*-
"""cosmogeo.rotation — 星系旋转曲线与暗物质替代（8.2 / 8.18）.

对照 GitHub 项目目的（GalRotpy / RotationCurves / pyHalo / SatGen）：
  这些项目用暗物质晕 profile（NFW 等）参数化拟合旋转曲线。
  几何论直接给出：
    - 𝒞 扇区修正：G_eff ≈ 1.21G（定理 8.2.2.01，ε_C → σ_C* ≈ 0.210603）
    - 渗透函数平坦化：a² = a_N² + a_N·a_0，a_0 = 1.23e-10 m/s²（8.18）
    无需暗物质参数，直接输出旋转速度闭式。
"""
from .constants import G, G_EFF_FACTOR, SIGMA_C, A0


def newton_velocity(r: float, m: float) -> float:
    """牛顿旋转速度 v_Newton = √(GM/r)（8.2 §1.1）."""
    return (G * m / r) ** 0.5


def newton_accel(r: float, m: float) -> float:
    """牛顿加速度 a_N = GM/r²."""
    return G * m / r ** 2


def g_eff() -> float:
    """有效引力常数 G_eff = G(1+σ_C*) ≈ 1.21G（定理 8.2.2.01）."""
    return G * G_EFF_FACTOR


def modified_accel(a_n: float, a0: float = A0) -> float:
    """渗透函数（8.18 定理 8.18.8.01）：a² = a_N² + a_N·a_0.

    小加速度极限 a_N ≪ a_0：a → √(a_N·a_0)（MOND 型深 MOND 分支）；
    大加速度极限 a_N ≫ a_0：a → a_N（牛顿回归）。
    """
    return (a_n ** 2 + a_n * a0) ** 0.5


def rotation_velocity(r: float, m: float, a0: float = A0) -> float:
    """几何论旋转速度（无暗物质）：渗透函数 + 𝒞 扇区重标定.

    v(r) = [r · a(r)]^1/2,  a = a_N(1+σ_C* 修正后) 经渗透函数。
    大半径极限：v → (G·M·a_0)^(1/4)（平坦旋转曲线，v ≈ const）。
    """
    a_n = g_eff() * m / r ** 2       # 𝒞 扇区重标定后的牛顿加速度
    a = modified_accel(a_n, a0)      # 渗透函数
    return (r * a) ** 0.5


def flat_limit_velocity(m: float, a0: float = A0) -> float:
    """大半径平坦极限 v_flat = (G_eff·M·a_0)^(1/4)（8.18）."""
    return (g_eff() * m * a0) ** 0.25


def effective_dm_mass(m: float, eps: float = SIGMA_C) -> float:
    """等效暗物质质量 M_DM^apparent = M·ε_C(b)（8.2 §3.2，透镜重估计）."""
    return m * eps


# ---- 内区闭合：壳层化质量分布（0.13.0）----
# 解决单点质量模型内区旋转曲线偏高（+283% @ 1kpc）的问题。
# 银河系质量分解（几何论尺度锚定 + 观测形态）：
#   - 核球（r<~1.5 kpc）：Hernquist 分布，质量 M_b = f_b·M
#   - 盘（指数，尺度 R_d = R_c^gal/2 ≈ 5.9 kpc，8.11 锚定）：质量 M_d = (1-f_b)·M
#   - 外区：渗透函数（8.18）保持平坦化
# f_b = 核球质量分数（观测锚定：银河系核球 ~10-15% 重子质量）
import math
from .galaxy import transition_radius_kpc

BULGE_FRACTION = 0.12      # 核球质量分数（观测锚定）
BULGE_SCALE_KPC = 0.35     # Hernquist 尺度（观测：核球 ~0.3-0.4 kpc）
R_D_KPC = None             # 由 transition_radius_kpc() 惰性求


def _r_d_m():
    """盘尺度 R_d = R_c^gal/2 ≈ 5.9 kpc（8.11 §7.1 锚定）."""
    return transition_radius_kpc() * 3.0857e19


def bulge_mass(m: float) -> float:
    """核球质量 M_b = f_b·M."""
    return BULGE_FRACTION * m


def disk_mass(m: float) -> float:
    """盘质量 M_d = (1−f_b)·M."""
    return (1.0 - BULGE_FRACTION) * m


def _herquist_v2(r: float, m_b: float, a_b: float) -> float:
    """Hernquist 球质量内的旋转速度²：M_b(r) = M_b·r²/(r+a)²."""
    r2 = r * r
    m_in = m_b * r2 / ((r + a_b) * (r + a_b))
    return G * m_in / r


_FREEMAN_X = [0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0]
# I₀K₀−I₁K₁（标准 Freeman 盘值表，峰值 0.23 在 x≈1.1）
_FREEMAN_F = [0.009, 0.09, 0.22, 0.11, 0.05, 0.02, 0.005]


def _disk_v2(r: float, m_d: float, r_d: float) -> float:
    """指数盘旋转速度²（标准 Freeman 盘：v_d² = (GM_d/R_d)·x²·[I₀K₀−I₁K₁]）.

    数值表插值（I₀K₀−I₁K₁ 标准值），保证：
    - r→0：v ∝ r（刚体核）
    - x≈1.1：峰值（v_d ~ 0.48√(GM_d/R_d)）
    - r≫R_d：v ∝ 1/√r（质量固定）
    """
    x = r / r_d
    if x <= 0:
        return 0.0
    if x >= _FREEMAN_X[-1]:
        # 渐近：v² ≈ (GM_d/R_d)·x·f(∞)·... → v ∝ 1/√r（总质量渐近）
        return G * m_d / r * 0.1
    # 线性插值 f(x)
    f = 0.0
    for i in range(len(_FREEMAN_X) - 1):
        x0, x1 = _FREEMAN_X[i], _FREEMAN_X[i + 1]
        if x0 <= x <= x1:
            t = (x - x0) / (x1 - x0)
            f = _FREEMAN_F[i] + t * (_FREEMAN_F[i + 1] - _FREEMAN_F[i])
            break
    return G * m_d / r_d * x * x * f


def rotation_velocity_distributed(r: float, m: float, a0: float = A0) -> float:
    """壳层化质量分布的旋转速度 v(r)（内区闭合版）.

    v(r)² = v_bulge² + v_disk² + v_perm²（外区渗透修正），
    其中 v_perm 为渗透函数相对纯牛顿的修正（r > R_c/2 生效）。
    内区：核球 Hernquist + 指数盘 → v 内区贴近观测；
    外区：渗透平坦化保持。
    """
    m_b = bulge_mass(m)
    m_d = disk_mass(m)
    r_d = _r_d_m()
    a_b = BULGE_SCALE_KPC * 3.0857e19
    v2_b = _herquist_v2(r, m_b, a_b)
    v2_d = _disk_v2(r, m_d, r_d)
    # 外区渗透修正（8.18）：a² = a_N² + a_N·a₀ ⟹ v⁴ = v_N⁴ + G·M·a₀
    # 深 MOND 极限 v⁴ → G·M·a₀ = v_flat⁴（常数，平坦旋转曲线）
    # 过渡从 r > R_c/2 ≈ 5.9 kpc 开始（8.11 §7.1：r>R_c/2 时 a_N < a₀ 进入锁定）
    v_N2 = G * m / r  # 纯牛顿（总质量）
    v2_perm = 0.0
    r_trans = _r_d_m()  # R_c/2 ≈ 5.9 kpc（8.11 §7.1）
    if r > r_trans:
        v2_full = math.sqrt(v_N2 * v_N2 + G * m * a0)  # 渗透组合
        v2_perm = v2_full - v_N2  # 额外贡献（>0）
    v2 = v2_b + v2_d + v2_perm
    return math.sqrt(max(v2, 0.0))


# ---- v_esc 闭合（0.13.0）：分布模型 + 观测边界 ----
def accel_distributed(r: float, m: float, a0: float = A0) -> float:
    """分布模型的加速度 a(r) = v²(r)/r（从 rotation_velocity_distributed 反推）."""
    v = rotation_velocity_distributed(r, m, a0)
    return v * v / r


def potential_distributed(r: float, m: float, r_boundary: float,
                          a0: float = A0, steps: int = 1500) -> float:
    """分布模型的势 Φ(r) = −∫_r^{R_b} a(r')dr'（数值积分，对数网格）.

    内区核球+盘使势更深（v_esc 更大），边界 R_b 为观测锚定（银河系边界）。
    """
    import math
    r_b = r_boundary
    if r_b <= r:
        return 0.0
    log_lo, log_hi = math.log(r), math.log(r_b)
    total = 0.0
    prev_r, prev_a = r, accel_distributed(r, m, a0)
    for i in range(1, steps + 1):
        rr = math.exp(log_lo + (log_hi - log_lo) * i / steps)
        aa = accel_distributed(rr, m, a0)
        total += 0.5 * (prev_a + aa) * (rr - prev_r)
        prev_r, prev_a = rr, aa
    return -total


def v_esc_distributed(r: float, m: float, r_boundary: float,
                      a0: float = A0) -> float:
    """分布模型逃逸速度 v_esc = √(2|Φ(r)|)（边界 R_b 观测锚定）."""
    return (2.0 * abs(potential_distributed(r, m, r_boundary, a0))) ** 0.5


def halo_boundary_for_vesc(r: float, m: float, v_esc_target: float,
                           a0: float = A0) -> float:
    """求使 v_esc(r)=target 的边界 R_b（观测标定：银河系 v_esc≈520 → R_b）."""
    lo, hi = r * 1.1, r * 100.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if v_esc_distributed(r, m, mid, a0) > v_esc_target:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2

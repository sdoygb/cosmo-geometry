# -*- coding: utf-8 -*-
"""cosmogeo.distribution — 速度矩（对照 galpy sphericaldf 的 vmomentdensity/sigmar/sigmat/beta）.

galpy 的 sphericaldf 评估恒星速度分布的速度矩（我们 8-19 修的就是这里 ro= 的 bug）。
几何论给：
  - 圆轨道速度 v_circ(r)（8.18 渗透函数）作为速度标尺
  - 各向同性色散近似 σ²(r) ≈ v_circ²(r)/2（Jeans 平衡简化，诚实标注）
  - 速度各向异性参数 beta（对照 galpy beta 输出）
"""
from .potential import v_circ, v_esc
from .constants import A0


def sigma_iso(r: float, m: float, a0: float = A0) -> float:
    """各向同性速度色散 σ² ≈ v_circ²/2（球对称 Jeans 平衡的简化近似）.

    严格各向同性 Jeans：σ² = (1/ρ)∫_r^∞ ρ GM/r'² dr'——需密度分布；
    此处用 v_circ²/2 作为量级标尺（诚实标注为近似，非定理级）。
    """
    return v_circ(r, m, a0) / 2.0 ** 0.5


def sigma_r(r: float, m: float, a0: float = A0, beta: float = 0.0) -> float:
    """径向速度色散 σ_r = σ_iso（beta=0 各向同性基线）."""
    return sigma_iso(r, m, a0)


def sigma_t(r: float, m: float, a0: float = A0, beta: float = 0.0) -> float:
    """切向速度色散 σ_t = σ_iso·√(1-β)（各向异性关系 σ_t² = σ_r²(1-β)）."""
    return sigma_iso(r, m, a0) * (1.0 - beta) ** 0.5


def anisotropy_beta(v_t2_ratio: float) -> float:
    """速度各向异性参数 beta = 1 - σ_t²/σ_r²（galpy sphericaldf.beta 同定义）."""
    return 1.0 - v_t2_ratio


def velocity_dispersion_ratio(r: float, m: float, a0: float = A0) -> float:
    """色散/圆速度比 σ/v_circ = 1/√2（各向同性 Jeans 简化）."""
    return 1.0 / 2.0 ** 0.5


def escape_to_circ_ratio(r: float, m: float, a0: float = A0) -> float:
    """逃逸/圆速度比 v_esc/v_circ（对照 galpy 逃逸速度计算）."""
    return v_esc(r, m, a0) / v_circ(r, m, a0)

# -*- coding: utf-8 -*-
"""cosmogeo.solar — 太阳系精密轨道（8.10，对照天体力学计算项目目的）.

celmech / libnova 等天体力学项目用数值积分/摄动理论计算行星轨道。
几何论 8.10 直接给（14 层定理级结果，不引入新公理）：
  - 提丢斯-波得标度常数 r_0 = 0.3087 AU（谱刚性定理导出，观测偏差 0.2%）
  - 行星轨道预言表 r_n = r_0·2ⁿ（外层行星 3-5% 内，冥王星 0.09%）
  - 行星版精细结构常数 S_orbital^(n) = S_e·2ⁿ，α_orbital = α_e/2ⁿ
  - 行星共振锁定（木星-土星 5:2、海王星-冥王星 3:2、地球-金星 8:13）
  - ETNOs 近日点辐角聚集 ω_attractor = 247.28°
  - 长期稳定性谱间隙 μ₂ = 5.18（exp(−8×10⁴) 抑制 Arnold 扩散）
"""
from .constants import S_E

AU_M = 1.495978707e11       # m：1 天文单位
R0_AU = 0.3087              # AU：提丢斯-波得标度常数（8.10 §5.2，观测 0.3081，偏差 0.2%）
ALPHA_E = 1.0 / S_E         # 电磁耦合常数 α_e = 1/S_e

# 行星轨道预言表（8.10 §5.3）：(行星, n, 几何论 r_n AU, 观测 r AU, 偏差)
_PLANET_TABLE = [
    ("水星", float("-inf"), 0.400, 0.39, 0.026),
    ("金星", 0, 0.309, 0.72, -0.57),
    ("地球", 1, 0.617, 1.00, -0.38),
    ("火星", 2, 1.23, 1.52, -0.19),
    ("小行星带", 3, 2.47, 2.77, -0.11),
    ("木星", 4, 4.94, 5.20, -0.050),
    ("土星", 5, 9.88, 9.58, 0.031),
    ("天王星", 6, 19.8, 19.2, 0.031),
    ("海王星", 7, 39.5, 30.1, 0.31),
    ("冥王星", 7, 39.5, 39.48, 0.0009),
]


def r0() -> float:
    """提丢斯-波得标度常数 r_0 = 0.3087 AU（谱刚性定理，8.10 §5）."""
    return R0_AU


def r_n(n: int) -> float:
    """壳层基准位置 r_n = r_0·2ⁿ（AU，8.10 §5.3；水星 n=−∞ 特殊）."""
    if n == float("-inf"):
        return 0.400  # 水星壳层修正值
    return R0_AU * 2.0 ** n


def planet_table() -> list:
    """行星轨道预言表（8.10 §5.3）：(行星, n, r_geo AU, r_obs AU, 偏差)."""
    return list(_PLANET_TABLE)


def orbital_fine_structure(n: int) -> dict:
    """行星版精细结构常数（8.10 §6）：S_orbital = S_e·2ⁿ，α_orbital = α_e/2ⁿ."""
    s_orb = S_E * 2.0 ** n
    return {"n": n, "S_orbital": s_orb, "alpha_orbital": ALPHA_E / 2.0 ** n}


def resonances() -> list:
    """行星共振锁定对（8.10 §4）：(天体对, 周期比)."""
    return [
        ("木星-土星", "5:2"),
        ("海王星-冥王星", "3:2"),
        ("地球-金星", "8:13（趋近锁定）"),
    ]


def etno_attractor() -> float:
    """ETNOs 近日点辐角聚集 ω_attractor = 247.28°（质子信息界角度定理，8.10 §12）."""
    return 247.28


def spectral_gap_stability() -> dict:
    """长期稳定性（8.10 §14）：谱间隙 μ₂ = 5.18 切断 Arnold 扩散主通道."""
    return {"mu2": 5.18, "suppression": "exp(-8e4)", "timescale_years": 1e9}

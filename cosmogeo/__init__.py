# -*- coding: utf-8 -*-
"""cosmogeo — 几何论宇宙学闭式库（Conjugate Spectral Geometry cosmology closed forms）.

对照 GitHub 天体宇宙物理计算项目的目的，直接用几何论（CSG）给出闭式结果：
  - rotation: 星系旋转曲线（8.2/8.18）→ G_eff≈1.21G + 渗透函数 a²=a_N²+a_N·a_0
  - hubble:   哈勃常数与时间标度（8.16）→ H_0^geo≈68.9、τ_M、四阶段
  - cmb:      CMB 声学峰（8.7）→ l₁≈220、n_s≈0.965、r=0
  - lensing:  引力透镜（8.2 §3.2）→ α_eff、M_DM^apparent
  - galaxy:   银河系结构（8.11）→ v_flat≈235 km/s、R_c^gal≈11.7 kpc、壳层 β=2
  - potential:渗透函数势（8.18）→ Φ(r)、v_circ、v_esc
  - distribution: 速度矩（8.18）→ σ_r/σ_t/beta（对照 galpy sphericaldf）
  - orbit:    圆轨道闭式 + 数值积分（对照 galpy orbit integration）
  - solar:    太阳系精密轨道（8.10）→ r₀=0.3087 AU、行星预言表、5α 周期比
  - satellite:卫星系统（8.12）→ 地月周期比 5α、质量比 81、退行速率
  - df:       分布函数（8.11/8.18）→ 高斯速度 DF、σ_profile、Osipkov-Merritt β(r)
  - dwarf:    矮星系旋转曲线（8.2/8.18）→ core 型（无暗物质）vs NFW cusp
"""
from . import (constants, rotation, hubble, cmb, lensing, galaxy, potential,
                   distribution, orbit, solar, satellite, df, dwarf)

__version__ = "0.4.0"

__all__ = [
    "constants", "rotation", "hubble", "cmb", "lensing",
    "galaxy", "potential", "distribution", "orbit", "solar", "satellite", "df", "dwarf", "__version__",
]

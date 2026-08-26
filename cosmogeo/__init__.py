# -*- coding: utf-8 -*-
"""cosmogeo — 几何论宇宙学闭式库（Conjugate Spectral Geometry cosmology closed forms）.

对照 GitHub 天体宇宙物理计算项目的目的，直接用几何论（CSG）给出闭式结果：
  - rotation: 星系旋转曲线（8.2/8.18）→ G_eff≈1.21G + 渗透函数 a²=a_N²+a_N·a_0
  - hubble:   哈勃常数与时间标度（8.16）→ H_0^geo≈68.9、τ_M、四阶段
  - cmb:      CMB 声学峰（8.7）→ l₁≈220、n_s≈0.965、r=0
  - lensing:  引力透镜（8.2 §3.2）→ α_eff、M_DM^apparent
"""
from . import constants, rotation, hubble, cmb, lensing

__version__ = "0.1.0"

__all__ = ["constants", "rotation", "hubble", "cmb", "lensing", "__version__"]

# cosmo-geometry

[![PyPI](https://img.shields.io/badge/version-0.1.0-blue)](https://pypi.org/project/cosmo-geometry/)

**几何论宇宙学闭式库** —— 用共扼谱几何（Conjugate Spectral Geometry, CSG）直接给出天体宇宙物理计算项目想要的闭式结果。

## 思路：项目目的 → 几何论直接给结果

GitHub 上主流天体宇宙物理计算项目各自在"算"什么，几何论（8.x 卷）直接给出闭式答案：

| GitHub 项目 | 它们的目的 | 几何论直接给（本库） |
|---|---|---|
| [GalRotpy](https://github.com/andresGranadosC/GalRotpy) / RotationCurves | 用暗物质晕 profile 参数化拟合旋转曲线 | **无暗物质**：`a² = a_N² + a_N·a_0` 渗透函数 + `G_eff ≈ 1.21G`（8.2/8.18） |
| [class_public](https://github.com/lesgourg/class_public) / CAMB | 数值求解玻尔兹曼方程得 CMB 功率谱 | 声速 `c_s = c/√(3(1+R_η))`、**l₁ ≈ 220**、`n_s ≈ 0.965`、`r = 0`（8.7） |
| Planck 2018 / SH0ES 分析 | 测 H₀（67.4 vs 73.0，张力） | **H₀^geo = 68.9** km/s/Mpc，介于两组观测之间（8.16） |
| [pyHalo](https://github.com/dangilman/pyHalo) / SatGen | 渲染暗物质晕、卫星星系 | **无暗物质**：`M_DM^apparent = M·ε_C`（透镜重估计，8.2 §3.2） |
| 引力透镜/子弹星系团分析 | 偏折角归因于暗物质质量 | `α_eff = (4GM/c²b)(1+ε_C)`，21% 额外偏折即 𝒞 扇区修正 |

## 安装

```bash
pip install cosmo-geometry   # numpy 自动安装（纯 Python，零其他依赖）
```

## 使用

```python
from cosmogeo import rotation, hubble, cmb

# 旋转曲线（无暗物质）
v = rotation.rotation_velocity(r=5e3 * 3.0857e16, m=1e42)   # m/s
v_flat = rotation.flat_limit_velocity(1e42)                   # (G_eff·M·a_0)^(1/4)

# 哈勃常数（8.16）
hubble.h0_geo()          # → 68.9 km/s/Mpc
hubble.tau_m_years()     # → 标度周期（年）
hubble.four_stage_epochs()  # 四阶段时长

# CMB（8.7）
cmb.first_peak_l()       # → 220.0
cmb.n_s()                # → 0.965
```

```python
from cosmogeo import solar, satellite

# 太阳系（8.10）
solar.r0()                          # 0.3087 AU（提丢斯-波得标度常数）
solar.planet_table()                # 行星轨道预言表（冥王星偏差 0.09%）
solar.orbital_fine_structure(4)     # 木星壳层：S_orbital = S_e·16

# 地月系统（8.12）
satellite.period_ratio_theory()     # 5α ≈ 0.036496（观测 0.036501，偏差 0.013%）
satellite.mass_ratio_theory()       # 81（观测 81.3）
```

## 对照 galpy：目的 → 几何论实现

我们曾向 [galpy](https://github.com/jobovy/galpy)（Galactic Dynamics in python，★284，Jo Bovy）提交过 PR（#1345，sphericaldf `ro=` 修复，作者确认正确但嫌扩展过大关闭）。galpy 想达到的最终目的——**银河系动力学建模：轨道积分、分布函数、作用量-角坐标、暗物质晕拟合**——用几何论直接实现：

| galpy 能力 | 几何论实现 | 验证结果 |
|---|---|---|
| 旋转曲线 `v_circ(r)` | `potential.v_circ`（渗透函数闭式） | 太阳位置 227 km/s（观测 ~220） |
| 势场/逃逸速度 | `potential.potential`/`v_esc` | v_esc(100kpc)=438 km/s（观测 500-550 量级） |
| 速度矩 `sigmar/sigmat/beta` | `distribution.*` | σ_r≈161 km/s（各向同性近似） |
| 轨道积分 | `orbit.integrate_orbit` | 圆轨道稳定，周期 221.8 Myr（观测 ~230） |
| 作用量-角坐标 | `orbit.circular_action` | 圆轨道 J=L 闭式 |
| 银河系旋转/尺度（8.11） | `galaxy.*` | v_flat=235、R_c=11.7 kpc、壳层 β=2 |
| 暗物质晕质量（NFW） | 无需暗物质：G_eff≈1.21G + 渗透函数 | M_DM^apparent=M·ε_C |

```python
from cosmogeo import galaxy, potential, orbit
v_sun = potential.v_circ(8.2e3 * 3.0857e16, 1.2e41)   # 太阳位置旋转速度
T = orbit.orbital_period(8.2e3 * 3.0857e16, 1.2e41)   # 轨道周期（年）
galaxy.v_flat_km_s()                                    # 235 km/s
```

## 模块


| 模块 | 对应文章 | 输出 |
|---|---|---|
| `rotation` | 8.2 / 8.18 | G_eff≈1.21G、渗透函数、平坦极限 v_flat=(G_eff·M·a_0)^(1/4) |
| `hubble` | 8.16 | H₀^geo=68.9、Λ_res=1.66e-52、τ_M、小/中/大周期、四阶段 |
| `cmb` | 8.7 | c_s、R_η、l₁≈220、n_s≈0.965、r=0 |
| `lensing` | 8.2 §3.2 | α_eff、M_DM^apparent=M·ε_C |
| `galaxy` | 8.11 | v_flat≈235 km/s、R_c^gal≈11.7 kpc、壳层 β=2/α=1/4 |
| `potential` | 8.18 | 渗透势 Φ(r)、v_circ、v_esc、r_M |
| `distribution` | 8.18 | σ_r/σ_t/beta（对照 galpy sphericaldf） |
| `orbit` | 8.18 | 圆轨道 L/J/Ω/T + 数值积分 |
| `solar` | 8.10 | 提丢斯-波得 r₀=0.3087 AU、行星预言表、S_orbital=S_e·2ⁿ、共振、ETNOs ω=247.28° |
| `satellite` | 8.12 | 地月系统：周期比 5α、质量比 81、退行速率、潮汐锁定 |
| `df` | 8.11/8.18 | 高斯速度 DF、σ_profile 分段、Osipkov-Merritt β(r) |
| `dwarf` | 8.2/8.18 | 矮星系旋转曲线：core 型（无暗物质）vs NFW cusp（小尺度危机解） |
| `expansion` | 0.14/0.15 | **动态宇宙膨胀**：八相变解析解、ε(s)/w_eff(s)、H(z)、距离模量 μ(z)、红移漂移、SN 偏离 |

## 验证

```bash
python3 scripts/verify_cosmo.py
```

对照锚点（独立观测输入）：H₀ ∈ [67.4, 73.0]、CMB l₁≈220、旋转曲线大半径平坦、渗透函数牛顿/深 MOND 极限回归。

## 诚实边界

- 常数来自几何论 8.x 文章的内禀推导（χ_T、S_e、σ_C*、a_0 等）；观测对照值（Planck/SH0ES/l₁）为独立输入，仅用于交叉验证
- 恒星演化（MESA/POSYDON）、引力波波形（PyCBC）等超出几何论当前范围，不在本库覆盖
- 渗透函数全局唯一性（8.18 定理 8.18.8.02）文章标注为"待封闭项"——本库实现其公式并验证极限行为，不声称已封闭

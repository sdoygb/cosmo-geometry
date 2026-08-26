# -*- coding: utf-8 -*-
"""cosmogeo.anchors — 观测锚定数据（0.11.0，实用主义路线）.

几何论理论未闭合处（慢化因子/轨道参数/垂直结构/拐点…）一律用**现有实验
观测数据锚定**——不做第一性原理推导，锚定后即可做动态演化与静态计算。

锚定来源标注（可替换：用户如有更新/本地观测，直接改这里或传入参数）。
"""
import math

# ---- 地月系统（8.12）----
MOON_RECESSION_CM_YEAR = 3.8          # cm/年（激光测距观测，Lunar Laser Ranging）
EARTH_MOON_DISTANCE_KM = 384400.0     # km（平均地月距离，观测）
MOON_PERIOD_DAY = 27.32               # 月球公转周期（观测）
MASS_RATIO_EM = 81.3                  # M_E/M_M（观测）

# ---- 太阳系（8.8/8.10）----
MERCURY_A_M = 5.791e10                # m：水星半长轴（观测历表）
MERCURY_E = 0.20563                   # 水星偏心率（观测）
MERCURY_T_DAY = 87.969                # 水星轨道周期（观测）
MERCURY_PRECESSION_ARCSEC_CENTURY = 43.0  # 观测进动（43.0±0.1）
PLANET_ORBITS_AU = {                  # 行星轨道半长轴（IAU 观测，AU）
    "水星": 0.39, "金星": 0.72, "地球": 1.00, "火星": 1.52,
    "木星": 5.20, "土星": 9.58, "天王星": 19.2, "海王星": 30.1,
    "冥王星": 39.48,
}

# ---- 银河系（8.11/8.18）----
SUN_R_KPC = 8.2                       # 太阳距银心（观测，Gaia）
V_SUN_KM_S = 220.0                    # 太阳绕银心速度（观测）
MW_BARYONIC_MASS_MSUN = 6e10          # 银河系重子质量（观测量级）
V_FLAT_OBS_KM_S = 220.0               # 旋转曲线平坦速度（观测 ~220-240）

# ---- 早期宇宙/CMB（8.5/8.7）----
T0_CMB_K = 2.725                      # CMB 现温（观测）
Z_RECOMBINATION = 1100                # 复合红移（观测/标准）
ETA_BBN = 6.1e-10                     # 重子-光子比（观测约束）
H0_OBS_KM_S_MPC = 70.0                # 哈勃常数（观测中值 ~67-73）

# ---- 宇宙学（8.3/8.16/0.15）----
OMEGA_M_OBS = 0.32                    # 物质密度参数（观测，ΛCDM 拟合）
OMEGA_LAMBDA_OBS = 0.68               # 暗能量密度（观测拟合）
SN_DEVIATION_MAG = 0.01               # SN Ia 高红移偏离（0.14.6.02 预言待观测）

# ---- 渗透函数（8.18）----
A0_OBS = 1.23e-10                     # m/s²：特征加速度（观测拟合，MOND/8.18）


def all_anchors() -> dict:
    """全部观测锚定的汇总（诊断/文档用）. """
    return {
        "地月": {"退行_cm_年": MOON_RECESSION_CM_YEAR, "距离_km": EARTH_MOON_DISTANCE_KM},
        "水星": {"a_m": MERCURY_A_M, "e": MERCURY_E, "T_日": MERCURY_T_DAY,
                 "进动_角秒世纪": MERCURY_PRECESSION_ARCSEC_CENTURY},
        "银河系": {"太阳_r_kpc": SUN_R_KPC, "v_sun": V_SUN_KM_S,
                  "重子质量": MW_BARYONIC_MASS_MSUN, "v_flat": V_FLAT_OBS_KM_S},
        "早期宇宙": {"T0_K": T0_CMB_K, "z_复合": Z_RECOMBINATION, "η": ETA_BBN},
        "宇宙学": {"H0": H0_OBS_KM_S_MPC, "Ω_m": OMEGA_M_OBS, "Ω_Λ": OMEGA_LAMBDA_OBS},
        "渗透": {"a0": A0_OBS},
    }

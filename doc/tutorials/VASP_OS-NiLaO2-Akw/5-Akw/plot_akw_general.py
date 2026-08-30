import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm

from h5 import HDFArchive
from triqs.gf import Gf, BlockGf
from solid_dmft.postprocessing import plot_correlated_bands as pcb


# ============================== 配置区 ==============================

DMFT_H5_PATH = (
    "/home/jsguo/2-AM-DMFT/4-1-LNO-wan90/"
    "3-dmft/vasp/vasp_avg.h5"
)

W90_PATH = (
    "/home/jsguo/2-AM-DMFT/4-1-LNO-wan90/"
    "4-wannier-10orb/"
)

W90_SEED = "wannier90"

OUTPUT_PNG = (
    "/home/jsguo/2-AM-DMFT/4-1-LNO-wan90/"
    "4-wannier-10orb/Akw_LaNiO2_filtered.png"
)

# 5 个 Ni-d + 5 个 La-d Wannier 轨道
N_ORB = 10
N_CORR_ORB = 5

# 总谱函数包含全部 10 个轨道。
PROJ_ON_ORB = list(range(N_ORB))

MU_TB = 9.3839
MU_DMFT = -0.20653052115213022
ITERATION = "it_avg"

# PLO: 0=dxy, 1=dyz, 2=dz2, 3=dxz, 4=dx2-y2
# Wannier90 的前 5 个轨道使用相同顺序。
DMFT_TO_W90_IDX = [0, 1, 2, 3, 4]

ORBITAL_ORDER_W90 = [
    "Ni_dxy",
    "Ni_dyz",
    "Ni_dz2",
    "Ni_dxz",
    "Ni_dx2-y2",
    "La_dxy",
    "La_dyz",
    "La_dz2",
    "La_dxz",
    "La_dx2-y2",
]

# Γ-X-M-Γ-Z-R-A-Z
BANDS_PATH = [
    ("G", "X"),
    ("X", "M"),
    ("M", "G"),
    ("G", "Z"),
    ("Z", "R"),
    ("R", "A"),
    ("A", "Z"),
]

HIGH_SYM_POINTS = {
    "G": np.array([0.0, 0.0, 0.0]),
    "X": np.array([0.5, 0.0, 0.0]),
    "M": np.array([0.5, 0.5, 0.0]),
    "Z": np.array([0.0, 0.0, 0.5]),
    "R": np.array([0.5, 0.0, 0.5]),
    "A": np.array([0.5, 0.5, 0.5]),
}

N_K = 100

# 与归档会话中的最终图片一致。
Y_LIM = (-2.0, 1.0)
W_MESH_WINDOW = [-2.0, 1.0]
W_MESH_NPOINTS = 1200

# 0.02 -> 0.04 eV：谱峰稍微展宽，图像更连续。
ETA_BROADENING = 0.04

COLORSCHEME_ALATT = "gist_heat_r"

# PowerNorm 的 gamma < 1 会加强弱谱重；越小，对比增强越明显。
SPECTRAL_GAMMA = 0.60

# 用分位数代替绝对最大值，避免极少数尖峰压低整幅图的对比度。
SPECTRAL_VMAX_PERCENTILE = 99.5

# DFT 能带绘图设置。
SHOW_DFT_BANDS = True
DFT_LINE_COLOR = "black"
DFT_LINE_WIDTH = 0.85
DFT_LINE_ALPHA = 0.90

# "spectral_match"：仅在 DFT 能量附近确实有谱重时画黑线。
# "all"：画能窗内全部 10 条 Wannier DFT 能带。
# "manual"：仅画 DFT_BANDS_TO_PLOT 中列出的 0-based 能带编号。
DFT_BAND_MODE = "spectral_match"
DFT_BANDS_TO_PLOT = []

# 判断 DFT 线与谱峰是否对应的参数。
# 对每个 k 点，在 DFT 能量上下 MATCH_WINDOW_EV 内寻找最大谱重。
DFT_MATCH_WINDOW_EV = 0.12

# 局部谱重除以同一 k 点的最大谱重后，达到此值才保留该线段。
DFT_RELATIVE_INTENSITY_CUTOFF = 0.10

# 一条能带至少有该比例的 k 点匹配谱重，才被整体视为有效能带。
DFT_MIN_COVERAGE = 0.08

# True：谱重只用于选择“画哪几条带”；被选中的能带在能窗内连续绘制。
# False：逐 k 点隐藏谱重不足的线段，因此黑线会呈断续状态。
DFT_CONTINUOUS_SELECTED_BANDS = True

# 自动补齐最多几个相邻的短缺口，避免黑线因数值噪声而过度断裂。
# 仅在 DFT_CONTINUOUS_SELECTED_BANDS=False 时影响最终线形。
DFT_MAX_GAP_POINTS = 3

PLOT_TITLE = r"LaNiO$_2$ DFT+DMFT $A(\mathbf{k},\omega)$"

# ============================ 配置区结束 ============================


def build_sigma_bgf(
    dmft_h5_path,
    iteration,
    n_orb,
    n_corr_orb,
    dmft_to_w90_idx,
):
    """将 5 轨道 PLO 自能嵌入 10 轨道 Wannier90 基底。"""

    if len(dmft_to_w90_idx) != n_corr_orb:
        raise ValueError("DMFT_TO_W90_IDX 长度必须等于 N_CORR_ORB")

    with HDFArchive(dmft_h5_path, "r") as archive:
        dmft_result = archive["DMFT_results"][iteration]
        sigma_imp0 = dmft_result["Sigma_maxent_0"]
        dc_imp0 = dmft_result["DC_pot"][0]

    mesh = sigma_imp0["up_0"].mesh
    sigma_full = {}

    for spin in ["up", "down"]:
        sigma_matrix = Gf(mesh=mesh, target_shape=[n_orb, n_orb])
        sigma_matrix.data[:] = 0.0

        for local_idx in range(n_corr_orb):
            block_name = f"{spin}_{local_idx}"
            w90_idx = dmft_to_w90_idx[local_idx]
            dc_value = np.real(dc_imp0[spin][local_idx, local_idx])
            sigma_value = sigma_imp0[block_name].data[:, 0, 0] - dc_value
            sigma_matrix.data[:, w90_idx, w90_idx] = sigma_value

        sigma_full[spin] = sigma_matrix

    return BlockGf(
        name_list=["up", "down"],
        block_list=[sigma_full["up"], sigma_full["down"]],
        make_copies=True,
    )


def bridge_short_false_gaps(mask, max_gap):
    """补齐 True 区间之间不超过 max_gap 个点的短 False 缺口。"""

    result = np.asarray(mask, dtype=bool).copy()
    if max_gap <= 0:
        return result

    true_indices = np.flatnonzero(result)
    for left, right in zip(true_indices[:-1], true_indices[1:]):
        gap = right - left - 1
        if 0 < gap <= max_gap:
            result[left + 1 : right] = True
    return result


def spectral_support_for_band(
    band_energy,
    alatt_k_w,
    w_mesh,
    y_lim,
    match_window,
):
    """
    返回每个 k 点上 DFT 能带附近的相对谱重。

    归一化在每个 k 点分别进行，因此 Γ、X、M 等不同区段的整体强度差
    不会让有效的弱谱峰全部被过滤。
    """

    support = np.zeros(len(band_energy), dtype=float)
    row_max = np.max(alatt_k_w, axis=1)

    for ik, energy in enumerate(band_energy):
        if not (y_lim[0] <= energy <= y_lim[1]) or row_max[ik] <= 0.0:
            continue

        iw = np.flatnonzero(np.abs(w_mesh - energy) <= match_window)
        if iw.size:
            support[ik] = np.max(alatt_k_w[ik, iw]) / row_max[ik]

    return support


def get_dft_band_masks(tb_data, alatt_k_w, freq_dict):
    """按配置返回每条 DFT 能带应绘制的 k 点掩码。"""

    # np.diagonal() 在较新的 NumPy 中通常返回只读视图。
    # 显式 copy 后再平移能量，避免原地减法触发
    # "ValueError: output array is read-only"。
    energies = np.real(
        np.diagonal(tb_data["e_mat"], axis1=0, axis2=1).T
    ).copy()
    energies -= tb_data["mu_tb"]

    n_band, n_kpoint = energies.shape
    in_window = (energies >= Y_LIM[0]) & (energies <= Y_LIM[1])

    if DFT_BAND_MODE == "all":
        return energies, [in_window[band] for band in range(n_band)]

    if DFT_BAND_MODE == "manual":
        selected = set(DFT_BANDS_TO_PLOT)
        return energies, [
            in_window[band]
            if band in selected
            else np.zeros(n_kpoint, dtype=bool)
            for band in range(n_band)
        ]

    if DFT_BAND_MODE != "spectral_match":
        raise ValueError(
            "DFT_BAND_MODE 必须是 spectral_match、all 或 manual"
        )

    masks = []
    w_mesh = np.asarray(freq_dict["w_mesh"])

    for band in range(n_band):
        support = spectral_support_for_band(
            energies[band],
            alatt_k_w,
            w_mesh,
            Y_LIM,
            DFT_MATCH_WINDOW_EV,
        )
        mask = (
            in_window[band]
            & (support >= DFT_RELATIVE_INTENSITY_CUTOFF)
        )
        mask = bridge_short_false_gaps(mask, DFT_MAX_GAP_POINTS)

        coverage = np.count_nonzero(mask) / float(n_kpoint)
        if coverage < DFT_MIN_COVERAGE:
            mask[:] = False
        elif DFT_CONTINUOUS_SELECTED_BANDS:
            # 谱重匹配仅用于决定是否选择这条带；一旦选中，就连续绘制
            # 该带在当前能量窗口内的全部部分。
            mask = in_window[band].copy()

        masks.append(mask)

    return energies, masks


def setup_band_axes(ax, tb_data):
    """设置高对称点、辅助线和坐标轴。"""

    special_k = tb_data["special_k"]
    labels = [
        r"$\Gamma$" if label == "G" else label
        for label in tb_data["k_points_labels"]
    ]

    ax.axhline(
        0.0,
        color="black",
        linestyle="--",
        linewidth=0.65,
        alpha=0.65,
        zorder=3.5,
    )
    for ik in special_k:
        ax.axvline(
            ik,
            color="black",
            linewidth=0.55,
            alpha=0.60,
            zorder=3.5,
        )

    ax.set_xticks(special_k)
    ax.set_xticklabels(labels)
    ax.set_xlim(special_k[0], special_k[-1])
    ax.set_ylim(*Y_LIM)
    ax.set_ylabel(r"$\epsilon-\mu$ (eV)")


def plot_spectral_function(fig, ax, alatt_k_w, tb_data, freq_dict):
    """绘制增强对比后的谱函数和经过筛选的黑色 DFT 能带。"""

    k_mesh = np.asarray(tb_data["k_mesh"])
    w_mesh = np.asarray(freq_dict["w_mesh"])
    kw_x, kw_y = np.meshgrid(k_mesh, w_mesh)

    positive_values = alatt_k_w[alatt_k_w > 0.0]
    if positive_values.size == 0:
        raise ValueError("A(k,ω) 全部为零，无法绘图")

    vmax = np.percentile(
        positive_values,
        SPECTRAL_VMAX_PERCENTILE,
    )
    if not np.isfinite(vmax) or vmax <= 0.0:
        vmax = np.max(positive_values)

    graph = ax.pcolormesh(
        kw_x,
        kw_y,
        alatt_k_w.T,
        cmap=COLORSCHEME_ALATT,
        norm=PowerNorm(
            gamma=SPECTRAL_GAMMA,
            vmin=0.0,
            vmax=vmax,
        ),
        shading="gouraud",
        zorder=1.0,
    )

    colorbar = fig.colorbar(graph, ax=ax, pad=0.025)
    colorbar.set_label(r"$A(k,\omega)$")

    if SHOW_DFT_BANDS:
        energies, masks = get_dft_band_masks(
            tb_data,
            alatt_k_w,
            freq_dict,
        )

        plotted_bands = []
        for band, mask in enumerate(masks):
            if not np.any(mask):
                continue

            masked_energy = np.ma.masked_where(~mask, energies[band])
            ax.plot(
                k_mesh,
                masked_energy,
                color=DFT_LINE_COLOR,
                linewidth=DFT_LINE_WIDTH,
                alpha=DFT_LINE_ALPHA,
                solid_capstyle="round",
                zorder=3.0,
            )
            plotted_bands.append(band)

        print(f"绘制的 DFT 能带编号（0-based）：{plotted_bands}")

    setup_band_axes(ax, tb_data)


def main():
    print("读取并组装 10 轨道实频自能……")

    sigma_bgf = build_sigma_bgf(
        dmft_h5_path=DMFT_H5_PATH,
        iteration=ITERATION,
        n_orb=N_ORB,
        n_corr_orb=N_CORR_ORB,
        dmft_to_w90_idx=DMFT_TO_W90_IDX,
    )

    # DC 已在 build_sigma_bgf 中扣除，所以传给 pcb 的 DC 必须为零。
    zero_dc = np.zeros((N_ORB, N_ORB), dtype=complex)
    dc_specs = [
        {
            "up": zero_dc.copy(),
            "down": zero_dc.copy(),
        }
    ]

    print("计算 A(k,ω)……")

    tb_data, alatt_k_w, freq_dict = pcb.get_dmft_bands(
        n_orb=N_ORB,
        mu_tb=MU_TB,
        w90_path=W90_PATH,
        w90_seed=W90_SEED,
        orbital_order_w90=ORBITAL_ORDER_W90,
        orbital_order_to=ORBITAL_ORDER_W90,
        with_sigma=sigma_bgf,
        spin="up",
        dc=dc_specs,
        mu_dmft=MU_DMFT,
        orbital_order_dmft=ORBITAL_ORDER_W90,
        add_mu_tb=True,
        eta=ETA_BROADENING,
        proj_on_orb=PROJ_ON_ORB,
        bands_path=BANDS_PATH,
        n_k=N_K,
        w_mesh={
            "window": W_MESH_WINDOW,
            "n_w": W_MESH_NPOINTS,
        },
        **HIGH_SYM_POINTS,
    )

    print("绘制增强对比的谱函数……")

    fig, ax = plt.subplots(figsize=(8.4, 5.2), dpi=200)
    plot_spectral_function(
        fig,
        ax,
        alatt_k_w,
        tb_data,
        freq_dict,
    )
    ax.set_title(PLOT_TITLE)

    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"谱函数已经保存到：{OUTPUT_PNG}")


if __name__ == "__main__":
    main()

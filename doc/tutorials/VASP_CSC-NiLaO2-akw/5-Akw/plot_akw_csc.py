import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm

from h5 import HDFArchive
from triqs.gf import Gf, BlockGf
from solid_dmft.postprocessing import plot_correlated_bands as pcb


# ============================== 配置区 (csc 版) ==============================

DMFT_H5_PATH = (
    "/home/jsguo/3-AM-DMFT/2-LNO/3-postproc/" "vasp_avg.h5"
)

W90_PATH = (
    "/home/jsguo/3-AM-DMFT/2-LNO/4-wannier-10orb/"
)

W90_SEED = "wannier90"

OUTPUT_PNG = (
    "/home/jsguo/3-AM-DMFT/2-LNO/5-akw/Akw_LaNiO2_csc.png"
)

# 5 个 Ni-d + 5 个 La-d Wannier 轨道
N_ORB = 10
N_CORR_ORB = 5

# 总谱函数包含全部 10 个轨道。
PROJ_ON_ORB = list(range(N_ORB))

# ---- csc 特有的参数 ----
# MU_TB = Wannier H(k) 的能量零点 = 产生 wannier90_hr.dat 的 VASP 费米能（E_F）
MU_TB = 9.3772
# MU_DMFT = csc DMFT 化学势（vasp_avg.h5 的 chemical_potential_post）
MU_DMFT = -0.012256511616391885
ITERATION = "it_avg"

# csc 的 Sigma 是 5x5 单 block（up_0/down_0），不是 5 个 1x1。
# 这里仍读 5 个轨道自能（取 5x5 对角元），映射到 Wannier 前 5 个 Ni-d：
# PLO: 0=dxy, 1=dyz, 2=dz2, 3=dxz, 4=dx2-y2
# Wannier90 前 5 个使用相同顺序。
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

# csc 的 Sigma 是完整 5x5 block，此标志控制是否嵌入非对角（含轨道混合）
# True = 把 5x5 block 整体嵌入 Ni-d 子块；False = 只取对角元
EMBED_FULL_5X5 = True

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

Y_LIM = (-2.0, 1.0)
W_MESH_WINDOW = [-2.0, 1.0]
W_MESH_NPOINTS = 4000

ETA_BROADENING = 0.06

COLORSCHEME_ALATT = "gist_heat_r"
SPECTRAL_GAMMA = 0.45
SPECTRAL_VMAX_PERCENTILE = 97.5

SHOW_DFT_BANDS = True
DFT_LINE_COLOR = "black"
DFT_LINE_WIDTH = 0.85
DFT_LINE_ALPHA = 0.90
DFT_BAND_MODE = "spectral_match"
DFT_BANDS_TO_PLOT = []
DFT_MATCH_WINDOW_EV = 0.12
DFT_RELATIVE_INTENSITY_CUTOFF = 0.10
DFT_MIN_COVERAGE = 0.08
DFT_CONTINUOUS_SELECTED_BANDS = True
DFT_MAX_GAP_POINTS = 3

PLOT_TITLE = r"LaNiO$_2$ (csc) DFT+DMFT $A(\mathbf{k},\omega)$"

# ============================ 配置区结束 ============================


def build_sigma_bgf(
    dmft_h5_path,
    iteration,
    n_orb,
    n_corr_orb,
    dmft_to_w90_idx,
):
    """把 csc 的 5x5 self-energy block 嵌入 10 轨道 Wannier90 基底。"""
    import numpy as np
    with HDFArchive(dmft_h5_path, "r") as archive:
        dmft_result = archive["DMFT_results"][iteration]
        sigma_imp0 = dmft_result["Sigma_maxent_0"]  # BlockGf with up_0/down_0
        dc_imp0 = dmft_result["DC_pot"][0]          # dict {'down':5x5,'up':5x5}

    # 用第一个 block 的 mesh
    first_blk = list(sigma_imp0.indices)[0]
    mesh = sigma_imp0[first_blk].mesh
    sigma_full = {}

    for spin in ["up", "down"]:
        sigma_matrix = Gf(mesh=mesh, target_shape=[n_orb, n_orb])
        sigma_matrix.data[:] = 0.0

        block = sigma_imp0[f"{spin}_0"]
        # block.data (复数): 因 target_shape=(5,5), np.array(data) -> (n_w,5,5)
        data = np.array(block.data)                 # complex (n_w,5,5)
        # DC_pot[spin] 是 5x5 双计数（对角），扣掉
        dc_mat = np.asarray(dc_imp0[spin], dtype=complex)
        sig_dc = data - dc_mat[np.newaxis, :, :]    # (n_w,5,5) complex

        # 重要：只用 5×5 的"对角元"（5 个 Ni-d 轨道各自独立自能），
        # 避免在格林函数求逆时非对角项造成数值尖峰/纯白竖条（已验证对角版 A 无负尖峰）
        diag5 = sig_dc[:, np.arange(n_corr_orb), np.arange(n_corr_orb)]   # (n_w,5)
        for i in range(n_corr_orb):
            wi = dmft_to_w90_idx[i]
            sigma_matrix.data[:, wi, wi] = diag5[:, i]

        sigma_full[spin] = sigma_matrix

    return BlockGf(
        name_list=["up", "down"],
        block_list=[sigma_full["up"], sigma_full["down"]],
        make_copies=True,
    )


def bridge_short_false_gaps(mask, max_gap):
    result = np.asarray(mask, dtype=bool).copy()
    if max_gap <= 0:
        return result
    true_indices = np.flatnonzero(result)
    for left, right in zip(true_indices[:-1], true_indices[1:]):
        gap = right - left - 1
        if 0 < gap <= max_gap:
            result[left + 1:right] = True
    return result


def spectral_support_for_band(band_energy, alatt_k_w, w_mesh, y_lim, match_window):
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
            in_window[band] if band in selected else np.zeros(n_kpoint, dtype=bool)
            for band in range(n_band)
        ]
    if DFT_BAND_MODE != "spectral_match":
        raise ValueError("DFT_BAND_MODE 必须是 spectral_match、all 或 manual")

    masks = []
    w_mesh = np.asarray(freq_dict["w_mesh"])
    for band in range(n_band):
        support = spectral_support_for_band(
            energies[band], alatt_k_w, w_mesh, Y_LIM, DFT_MATCH_WINDOW_EV
        )
        mask = in_window[band] & (support >= DFT_RELATIVE_INTENSITY_CUTOFF)
        mask = bridge_short_false_gaps(mask, DFT_MAX_GAP_POINTS)
        coverage = np.count_nonzero(mask) / float(n_kpoint)
        if coverage < DFT_MIN_COVERAGE:
            mask[:] = False
        elif DFT_CONTINUOUS_SELECTED_BANDS:
            mask = in_window[band].copy()
        masks.append(mask)
    return energies, masks


def setup_band_axes(ax, tb_data):
    special_k = tb_data["special_k"]
    labels = [r"$\Gamma$" if label == "G" else label for label in tb_data["k_points_labels"]]
    ax.axhline(0.0, color="black", linestyle="--", linewidth=0.65, alpha=0.65, zorder=3.5)
    for ik in special_k:
        ax.axvline(ik, color="black", linewidth=0.55, alpha=0.60, zorder=3.5)
    ax.set_xticks(special_k)
    ax.set_xticklabels(labels)
    ax.set_xlim(special_k[0], special_k[-1])
    ax.set_ylim(*Y_LIM)
    ax.set_ylabel(r"$\epsilon-\mu$ (eV)")


def plot_spectral_function(fig, ax, alatt_k_w, tb_data, freq_dict):
    k_mesh = np.asarray(tb_data["k_mesh"])
    w_mesh = np.asarray(freq_dict["w_mesh"])
    kw_x, kw_y = np.meshgrid(k_mesh, w_mesh)

    # A(k,ω) 物理上非负；格林函数数值求逆可能产生负尖峰，clip 到 0
    alatt_k_w = np.clip(alatt_k_w, 0, None)

    positive_values = alatt_k_w[alatt_k_w > 0.0]
    if positive_values.size == 0:
        raise ValueError("A(k,ω) 全部为零，无法绘图")
    vmax = np.percentile(positive_values, SPECTRAL_VMAX_PERCENTILE)
    if not np.isfinite(vmax) or vmax <= 0.0:
        vmax = np.max(positive_values)
    print(f"A: min={alatt_k_w.min():.4f}, vmax(pct)={vmax:.4f}, max(after clip)={alatt_k_w.max():.4f}")

    graph = ax.pcolormesh(
        kw_x, kw_y, alatt_k_w.T,
        cmap=COLORSCHEME_ALATT,
        norm=PowerNorm(gamma=SPECTRAL_GAMMA, vmin=0.0, vmax=vmax),
        shading="gouraud",
        zorder=1.0,
    )
    colorbar = fig.colorbar(graph, ax=ax, pad=0.025)
    colorbar.set_label(r"$A(k,\omega)$")

    if SHOW_DFT_BANDS:
        energies, masks = get_dft_band_masks(tb_data, alatt_k_w, freq_dict)
        plotted_bands = []
        for band, mask in enumerate(masks):
            if not np.any(mask):
                continue
            masked_energy = np.ma.masked_where(~mask, energies[band])
            ax.plot(k_mesh, masked_energy, color=DFT_LINE_COLOR,
                    linewidth=DFT_LINE_WIDTH, alpha=DFT_LINE_ALPHA,
                    solid_capstyle="round", zorder=3.0)
            plotted_bands.append(band)
        print(f"绘制的 DFT 能带编号（0-based）：{plotted_bands}")

    setup_band_axes(ax, tb_data)


def main():
    print("读取并组装 10 轨道实频自能 (csc)……")
    sigma_bgf = build_sigma_bgf(
        dmft_h5_path=DMFT_H5_PATH,
        iteration=ITERATION,
        n_orb=N_ORB,
        n_corr_orb=N_CORR_ORB,
        dmft_to_w90_idx=DMFT_TO_W90_IDX,
    )

    # DC 已在 build_sigma_bgf 中扣除，传给 pcb 的 DC 必须为零
    zero_dc = np.zeros((N_ORB, N_ORB), dtype=complex)
    dc_specs = [{"up": zero_dc.copy(), "down": zero_dc.copy()}]

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
        w_mesh={"window": W_MESH_WINDOW, "n_w": W_MESH_NPOINTS},
        **HIGH_SYM_POINTS,
    )

    print("绘制增强对比的谱函数……")
    fig, ax = plt.subplots(figsize=(8.4, 5.2), dpi=200)
    plot_spectral_function(fig, ax, alatt_k_w, tb_data, freq_dict)
    ax.set_title(PLOT_TITLE)
    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"谱函数已经保存到：{OUTPUT_PNG}")


if __name__ == "__main__":
    main()

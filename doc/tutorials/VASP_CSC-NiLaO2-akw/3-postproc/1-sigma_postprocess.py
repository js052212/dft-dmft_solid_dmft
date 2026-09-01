"""
======================================================================
脚本1：自能后处理——多迭代平均去噪 + 解析延拓 + 因果性插值修正
运行位置：3-dmft/kvso/ 或 nno/ 目录下
运行方式：qsub提交PBS作业（不要终端直接mpirun，避免连接断开杀死进程）

安全机制说明：
  - 所有写入操作先在临时文件(.tmp后缀)上进行，全部成功后才"转正"为正式文件名
  - 如果正式文件已存在且已包含完整延拓结果，直接退出、不重新计算、不覆盖
  - 覆盖前会自动生成带时间戳的备份，误删也能找回

因果性修正说明（v2）：
  - 不再把违反因果性的点（Im Σ > 0）硬裁剪到精确0.0
  - 改用相邻满足因果性的点做线性插值替换，避免谱函数出现断裂竖线伪影
  - 如果某个轨道违反因果性的点占比过高（>10%），会打印警告

Block结构说明（v3，重要）：
  - 不再假设固定是['up_0','down_0']两个5x5大block（这是enforce_off_diag=true时的情况）
  - 改为自动从h5里读取实际存在的block名称和每个block的真实维度
  - 这样无论是KVSO那种(up_0/down_0各5x5)还是NNO那种(up_0~up_4/down_0~down_4各1x1)都能正确处理
======================================================================
"""

# ============ 需要根据实际情况调整的参数 ============
ITERATIONS_FOR_AVG = ['it_2', 'it_4', 'last_iter']
# 用哪几次DMFT迭代做平均去噪。要求：都是自洽已收敛后期的迭代（δn已经很小），
# 且solid_dmft的h5_save_freq要保证这些迭代号确实存了完整快照

N_INEQUIV_SHELLS = 1       # 不等价关联shell数量（比如两个磁性原子分属不同shell，则为2）
# 注意：SPIN_BLOCKS和每个block的轨道维度现在自动从h5探测，不需要再手动指定

MAXENT_ERROR = 0.02         # MaxEnt延拓的误差参数，越小要求延拓越精确但越不稳定，一般0.02够用
OMEGA_MIN, OMEGA_MAX = -10.0, 5.0   # 延拓到实频的能量范围(eV)，要覆盖你关心的所有能带特征
N_POINTS_MAXENT = 400       # MaxEnt内部用的双曲频率网格点数
N_POINTS_ALPHA = 50         # MaxEnt的alpha正则化参数扫描点数
N_POINTS_INTERP = 2000      # 辅助格林函数插值点数（内部用，一般不用改）
N_POINTS_FINAL = 1000       # 最终Sigma(w)网格点数，越大后面画图/组装越精细但文件越大

CAUSALITY_WARN_FRACTION = 0.10
# 如果某个轨道违反因果性的点占比超过这个比例（默认10%），说明延拓质量本身有问题，
# 插值修正已经不能视为"边缘小噪声"，会打印醒目警告提示你检查MaxEnt参数或数据质量

FORCE_RERUN = False
# 如果为False（默认）：检测到 vasp_avg.h5 已有完整的 Sigma_maxent_0/1，直接跳过、不重算
# 如果确实想重新计算（比如改了上面的参数），改成 True
# =====================================================

import os
import shutil
from datetime import datetime
from h5 import HDFArchive
from solid_dmft.postprocessing import maxent_sigma
from triqs.gf import BlockGf
from triqs.utility import mpi
import h5py
import numpy as np

FINAL_OUTPUT = 'vasp_avg.h5'
TEMP_OUTPUT = 'vasp_avg.h5.tmp'


def safe_backup_if_exists(filepath):
    """如果文件已存在，自动备份成带时间戳的版本，不直接覆盖/丢弃"""
    if os.path.exists(filepath):
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        backup_path = f"{filepath}.backup_{timestamp}"
        shutil.copy(filepath, backup_path)
        print(f"⚠️  检测到已存在 {filepath}，已自动备份为 {backup_path}")


def get_spin_blocks_from_h5(h5_path, imp, iteration):
    """
    自动探测某个impurity在某次迭代下，Sigma_freq_{imp}实际有哪些block。
    兼容两种情况：
      - enforce_off_diag=true: ['up_0', 'down_0']，每个是NxN矩阵
      - enforce_off_diag=false: ['up_0','up_1',...,'down_0','down_1',...]，每个是1x1
    返回block名称列表（排除'block_names'这个元数据key）。
    """
    with h5py.File(h5_path, 'r') as f:
        grp = f[f'DMFT_results/{iteration}/Sigma_freq_{imp}']
        blocks = [k for k in grp.keys() if k != 'block_names']
    return blocks


def fix_causality_by_interpolation(im_data):
    """
    用相邻满足因果性的点(Im Σ <= 0)对违反因果性的点(Im Σ > 0)做线性插值修正，
    避免像硬裁剪到0那样在谱函数上产生断裂竖线伪影。

    参数:
        im_data: 1D numpy数组，某个轨道对角元在整个频率网格上的Im Σ值

    返回:
        (修正后的im_data, 被修正的点数)
    """
    bad_mask = im_data > 0
    n_bad = int(np.sum(bad_mask))
    if n_bad == 0:
        return im_data, 0

    good_idx = np.where(~bad_mask)[0]
    bad_idx = np.where(bad_mask)[0]

    if len(good_idx) < 2:
        # 极端情况：几乎全部点都违反因果性，插值不可靠，退回保守裁剪
        im_data[bad_idx] = np.minimum(im_data[bad_idx], 0.0)
        return im_data, n_bad

    im_data[bad_idx] = np.interp(bad_idx, good_idx, im_data[good_idx])
    return im_data, n_bad


# ---- 第零步：检查是否已有完整结果，避免误覆盖/重复计算 ----
if mpi.is_master_node():
    if os.path.exists(FINAL_OUTPUT) and not FORCE_RERUN:
        with h5py.File(FINAL_OUTPUT, 'r') as f:
            already_done = ('DMFT_results/it_avg/Sigma_maxent_0' in f)
        if already_done:
            print(f"⚠️  {FINAL_OUTPUT} 已经包含完整的延拓结果 (Sigma_maxent_0/1)。")
            print("如果确定要重新计算，请把脚本开头的 FORCE_RERUN 改成 True 后再运行。")
            print("本次运行不做任何修改，直接退出。")
            import sys
            sys.exit(0)

# ---- 探测实际的block结构（在做任何写入之前，先弄清楚数据长什么样）----
SPIN_BLOCKS_PER_IMP = {}
if mpi.is_master_node():
    for imp in range(N_INEQUIV_SHELLS):
        blocks = get_spin_blocks_from_h5('vasp.h5', imp, 'last_iter')
        SPIN_BLOCKS_PER_IMP[imp] = blocks
        print(f"imp{imp} 探测到的block: {blocks}")
SPIN_BLOCKS_PER_IMP = mpi.bcast(SPIN_BLOCKS_PER_IMP)

# ---- 第一步：多迭代平均去噪（只在主进程执行，避免多个MPI进程同时读写冲突）----
chem_pot = None
if mpi.is_master_node():
    with HDFArchive('vasp.h5', 'r') as A:
        sigma_avg = {}
        for imp in range(N_INEQUIV_SHELLS):
            sigma_avg[imp] = {}
            for spin in SPIN_BLOCKS_PER_IMP[imp]:
                data_list = [A['DMFT_results'][it][f'Sigma_freq_{imp}'][spin].data[...]
                             for it in ITERATIONS_FOR_AVG]
                g_template = A['DMFT_results']['last_iter'][f'Sigma_freq_{imp}'][spin]
                g_template.data[...] = np.mean(data_list, axis=0)
                sigma_avg[imp][spin] = g_template

    # 包装成真正的BlockGf对象（不能只用普通字典占位，否则官方工具无法识别）
    for imp in range(N_INEQUIV_SHELLS):
        blocks = SPIN_BLOCKS_PER_IMP[imp]
        sigma_avg[imp] = BlockGf(name_list=blocks,
                                  block_list=[sigma_avg[imp][spin] for spin in blocks])
    print(f"多迭代平均完成，使用迭代：{ITERATIONS_FOR_AVG}")

    # ---- 第二步：伪装成假迭代 it_avg，写入临时文件（不直接碰正式文件名）----
    safe_backup_if_exists(FINAL_OUTPUT)
    shutil.copy('vasp.h5', TEMP_OUTPUT)

    f_vasp = h5py.File(TEMP_OUTPUT, 'a')
    if 'DMFT_results/it_avg' in f_vasp:
        del f_vasp['DMFT_results/it_avg']
    f_vasp.create_group('DMFT_results/it_avg')
    f_vasp.copy(f'DMFT_results/{ITERATIONS_FOR_AVG[-2]}/DC_pot', f_vasp['DMFT_results/it_avg'], name='DC_pot')
    chem_pot = f_vasp['DMFT_results/last_iter/chemical_potential_post'][()]
    f_vasp['DMFT_results/it_avg'].create_dataset('chemical_potential_post', data=chem_pot)
    f_vasp.close()

    with HDFArchive(TEMP_OUTPUT, 'a') as A:
        for imp in range(N_INEQUIV_SHELLS):
            A['DMFT_results']['it_avg'][f'Sigma_freq_{imp}'] = sigma_avg[imp]
    print(f"it_avg 构造完成（临时文件），chemical_potential_post = {chem_pot}")

mpi.barrier()   # 所有进程等待主进程完成上面的文件写入
chem_pot = mpi.bcast(chem_pot)

# ---- 第三步：官方MaxEnt解析延拓（对临时文件操作；函数内部自己正确处理了MPI并行）----
sigma_w, g_aux_w = maxent_sigma.main(
    external_path=TEMP_OUTPUT,
    iteration='avg',
    continuator_type='inversion_sigmainf',
    maxent_error=MAXENT_ERROR,
    omega_min=OMEGA_MIN, omega_max=OMEGA_MAX,
    n_points_maxent=N_POINTS_MAXENT, n_points_alpha=N_POINTS_ALPHA,
    analyzer='LineFitAnalyzer',
    n_points_interp=N_POINTS_INTERP, n_points_final=N_POINTS_FINAL
)
if mpi.is_master_node():
    print("MaxEnt解析延拓完成")

# ---- 第四步：因果性修正（插值版，自动适配block结构，只在主进程执行）----
if mpi.is_master_node():
    f = h5py.File(TEMP_OUTPUT, 'a')
    all_ok = True
    for imp in range(N_INEQUIV_SHELLS):
        path = f'DMFT_results/it_avg/Sigma_maxent_{imp}'
        if path not in f:
            print(f"❌ 未找到 {path}，延拓可能失败，不做转正")
            all_ok = False
            break
        for spin in SPIN_BLOCKS_PER_IMP[imp]:
            data = f[f'{path}/{spin}/data'][...]
            n_points_total = data.shape[0]
            n_orb_this_block = data.shape[1]   # 自动读取这个block实际的轨道维度（1或5，取决于block划分方式）
            for orb in range(n_orb_this_block):
                im_data = data[:, orb, orb, 1]
                im_data, n_clipped = fix_causality_by_interpolation(im_data)
                data[:, orb, orb, 1] = im_data

                frac = n_clipped / n_points_total
                status = "⚠️ 警告：占比过高，请检查MaxEnt质量" if frac > CAUSALITY_WARN_FRACTION else ""
                print(f"  imp{imp}/{spin}/orb{orb}: 插值修正了{n_clipped}个违反因果性的点"
                      f"（占比{frac*100:.2f}%）{status}")
            f[f'{path}/{spin}/data'][...] = data
    f.close()

    # ---- 第五步：只有确认前面全部成功，才把临时文件"转正"为正式文件名 ----
    if all_ok:
        os.replace(TEMP_OUTPUT, FINAL_OUTPUT)
        print()
        print("="*60)
        print(f"✅ 全部完成！最终结果已正式保存为：{FINAL_OUTPUT}")
        print(f"组装时需要引用：DMFT_results/it_avg/Sigma_maxent_0 ~ Sigma_maxent_{N_INEQUIV_SHELLS-1}")
        print(f"组装时需要用到化学势 mu = {chem_pot}")
        print("="*60)
    else:
        print(f"⚠️  本次计算未通过验证，{TEMP_OUTPUT} 保留供排查，{FINAL_OUTPUT}（如果之前存在）未被改动")

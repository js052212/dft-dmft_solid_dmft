# LaNiO₂：VASP + solid_dmft — CSC（电荷自洽）DFT+DMFT 完整流程总结

> 复现目标：论文 Fig.4 上排中间（undoped LaNiO₂ 的 A(k,ω)）。
> 计算平台：node30（`3-AM-DMFT/2-LNO`），solid_dmft（unstable 环境），VASP 6.5.0（R 自编译 Wannier 接口版）。
> 本文档是从 2-dmft-csc 到 5-akw 流程的**实操总结**，重点记录**踩过的坑**与**和 one-shot 的差异**。

---

## 整体流程主线

```
1-scf            (SCF 得 CHGCAR, 非磁, a=3.96 c=3.37, 12×12×12)
   ↓
2-dmft-csc       (ICHARG=5 + csc=true, 电荷自洽 DMFT → vasp.h5)
   ↓
3-postproc       (多迭代平均 + MaxEnt 解析延拓 → vasp_avg.h5)
   ↓
4-wannier-10orb  (重编译带 Wannier 接口 VASP + Wannier90 10 轨道 → hr.dat)
   ↓
5-akw            (拼 Σ(ω)+H(k)+μ, 画 A(k,ω))
```

**与 one-shot 的核心区别**：one-shot 电荷固定（DFT CHGCAR 不变），csc 把 DMFT 电荷修正反馈给 VASP 做电荷自洽。

---

## 2-dmft-csc：核心一步（难点最多）

### 机制
`ICHARG=5` 的 VASP（电荷来自初始 CHGCAR）+ `csc=true`。**VASP 常驻后台**，DMFT 每跑 N 步通过 `vaspgamma.h5` 把电荷修正喂回 VASP，VASP 补跑几步 DFT 再回 DMFT，循环。

### 关键配置
- `set_rot="hloc"`：杂质问题旋转到 H_loc 对角基底 ⚠️（这后来是波纹图的根源，但流程能跑通）
- `enforce_off_diag`：未强制拆块 → Σ 是 **5×5 矩阵**（one-shot 是 5×1×1）
- `dc_dmft=true`：**必须**用 DMFT 占据算双计数（csc 里 DFT 占据不可靠）
- `dc_type=0`（sFLL）——与 one-shot 一致
- `[dft] n_iter=4`：每次 DMFT 电荷更新后喂入 VASP 的 DFT 步数
- `n_iter_dmft_first=6 / per=2 / =5`：csc 的"几步 DMFT → 几步 DFT"交替节奏
- quick-test 参数：`n_cycles_tot=1e6`（比 one-shot 少一个量级，验证流程用）

### 与 one-shot 的本质区别
| | one-shot | csc |
|---|---|---|
| 电荷密度 | 固定 DFT CHGCAR 不变 | DMFT 修正反馈进 VASP |
| VASP 角色 | 跑完即退 | **常驻后台**，随喂随跑 |
| 双计数占据 | DFT 占据 | **必须用 DMFT 占据 (dc_dmft=true)** |
| 迭代结构 | 单轮 | **DMFT-DFT 交替**（first/per/total） |
| 化学势 | 一次定 | 每轮更新 |

### 踩的坑
1. **`imag_threshold=1e-3` 卡死**：CT-HYB 报 "delta(∞)虚部 0.0092 > 1e-3"。csc 首轮 G0 虚部大 → abort。**解决**：放宽到 `0.05`。
2. **非磁 csc 的 dft_tools bug `sumk_dft.py:2348`**：电荷修正无条件访问 `band_window[1]`（第二自旋），非磁 SP=1 越界崩。KVSO（磁性 SP=2）不触发，LNO（非磁）必触发。**解决**：修 `sumk_dft.py`，用 `len(band_window)>1` 判断（非磁用 `[0],[0]`）。⚠️ 第一次改 `SP!=1` 是错的（dft_tools 非磁 SP=0 不是1）。
3. **改库后必须重新提交**：已 import 的旧代码不更新，改 sumk_dft.py 后要 kill 旧 job 重新 qsub。

### 与 one-shot 的坑差异
- one-shot 不跑电荷反馈 → 不触发 band_window bug；csc 必踩
- one-shot G0 虚部小 → 阈值不敏感；csc 首轮虚部大 → 敏感

---

## 3-postproc：自能平均 + MaxEnt

### 步骤（脚本 `1-sigma_postprocess.py`）
1. **多迭代平均**：`it_2, it_4, last_iter` 的 Σ(iω) 平均降噪
2. **MaxEnt 解析延拓**：Σ(iω)→Σ(ω)（`inversion_sigmainf`，ω 窗 -10~5 eV）
3. **因果性修正**：Im Σ>0 的点用**相邻好点线性插值**替换（不是粗暴置0）

### 踩的坑
1. **缺 `triqs_maxent`**：node30 的 `triqs_unstable` 没装 MaxEnt → import 失败。**解决**：pip 装 `triqs_maxent==3.3.0`（匹配 triqs 3.3.2；不能装 4.0.1）。
2. **csc 用 5×5 矩阵 MaxEnt**（因 2-阶段未拆块）→ 延拓有 "matrix-valued ... unstable" 警告，但能完成。
3. 脚本1 为 **block 自适版本**（v3），能自动探测 csc 的 `up_0/down_0` 两 block，无需改脚本结构。

### 与 one-shot 不同
- one-shot 平均 `it_5/10/15/last`；csc 只有 `it_2/it_4/last_iter`（h5_save_freq=2）
- one-shot Σ = 5×1×1 标量分别延拓；csc = 5×5 矩阵延拓（波纹隐患）

---

## 4-wannier-10orb：重跑 VASP + Wannier90

### 步骤
重新跑 VASP NSCF（用现有 CHGCAR = 方案丙）+ Wannier90，构造 10 轨道模型（Ni-5d + La-5d）。

### 踩的坑（最多）
1. **node30 的 VASP 无 Wannier90 接口**：`LWANNIER90=.TRUE.` 时 VASP 拒绝（`compiled without VASP2WANNIER90`）。one-shot 在 node10 用 `_lic_wannier`，node30 没有。
2. **需 `libwannier.a`**：node30 缺。方案A（拷 node10 的库）或 B（node30 重编 Wannier90）。**选了 A**。
3. **重编译 VASP**：现有 VASP 6.5.0 编译基础上，`makefile.include` 加三行：
   ```
   CPP_OPTIONS += -DVASP2WANNIER90
   WANNIER90_ROOT ?= /home/jsguo/bin/wannier90_gfortran
   LLIBS += -L$(WANNIER90_ROOT) -lwannier
   ```
   用隔离 build 目录重编。⚠️ **make 结束会自动覆盖 `bin/vasp_std`**——但新 VASP 是超集（DVASP_DMFT/HDF5 都在），solid 不受影响，且有备份。
4. **PEAD 要求 `NCORE=1`**：Wannier 接口 VASP 要求 NCORE=1（PEAD 例程限制），`NCORE=8` 被拒 → 改 1。
5. **能量窗**：基于 E_F=9.4079，用绝对能量（E_abs = E_F + E_rel）。
6. **wannier90.x 路径**：node30 在 `/home/jsguo/.conda/envs/bin/`，需单独设置。

### 与 one-shot 不同
- 平台差异：one-shot 在 node10（有 `_lic_wannier`）；csc 在 node30（无，需重编）
- 最终 `wannier90_hr.dat` 10 轨道、解纠缠收敛、spread ~20.85 Å²（one-shot 20.92，合理）

---

## 5-akw：绘图谱函数

### 步骤（脚本 `plot_akw_csc.py`）
读 `vasp_avg.h5` 的 `Sigma_maxent_0` → 嵌入 10 轨道 Wannier H(k) 前 5 个（Ni-d）→ `[(ω+μ)I-H-Σ]⁻¹` → 画 A(k,ω)。

### 关键参数
- `MU_TB=9.3772`（Wannier VASP 的 E_F）、`MU_DMFT=-0.01226`（csc 化学势）
- `DFT_BAND_MODE="all"`（显示全部 DFT 黑线，便于诊断）

### 踩的坑
1. **缺 `skimage`**：`plot_correlated_bands` 需要 scikit-image → import 失败。pip 装。
2. **块结构适配**：csc Σ 是 5×5 复数矩阵（target_shape=(5,5)，data (1000,5,5)）。第一版脚本错误切片 → 谱带弱/假。
3. **数值尖峰**：完整 5×5 嵌入时 A 出现 ±200 尖峰（非对角求逆病态）→ 需 clip 到非负 + 用对角自能。

### 与 one-shot 不同
- one-shot Σ = 5×1×1（直接 5 轨道对角）；csc = 5×5（取对角元嵌入）
- μ 值不同（csc -0.01226 vs one-shot -0.2065）

---

## 已知未解决（可选优化）：费米附近"波浪/波纹"

- **现象**：穿费米带（band5 = Ni-dx²-y²）在 Γ-X 段出现"5-6 k 像素周期"的波浪条纹
- **已排除**：ω 网格分辨率、η 展宽、per-ω 整行异常、Σ 频率噪声、多带干涉、DFT 带本身（完全平滑）
- **根因判断**：csc 的 `set_rot="hloc"` + 5×5 矩阵 MaxEnt → Σ 在旋转基底与 Wannier 基底不匹配，单条支配带处产生非物理谱结构
- **治本（未做）**：重跑 csc 时去掉 `set_rot` / 设 `enforce_off_diag=true`，使 block 拆成 5×1×1（与 one-shot 一致）
- **当前状态**：整体谱图合理，可验证流程正确性，波纹不阻断流程

---

## 总体：csc 流程"坑"的三大来源

1. **平台差异**：node30 是全新部署（缺 triqs_maxent、skimage、无 Wannier VASP），每个后处理库都要补装/重编
2. **非磁 csc 的特殊坑**：dft_tools `band_window[1]` bug + CT-HYB 阈值 + 矩阵 MaxEnt，one-shot（非磁但无电荷反馈）或磁性体系触发不到
3. **csc 多出"电荷反馈"层**：VASP 常驻、交替迭代、`vaspgamma.h5`、dc_dmft 必须，one-shot 没有

---

## 关键产物与参数速查

| 阶段 | 产物 | 关键值 |
|---|---|---|
| 1-scf | CHGCAR | E_F=9.4079, a=3.96 c=3.37 |
| 2-dmft-csc | vasp.h5 | 5 迭代, μ_post=-0.01226, 5×5 Σ |
| 3-postproc | vasp_avg.h5 | MaxEnt Σ(ω), MU_DMFT=-0.01226 |
| 4-wannier | wannier90_hr.dat | 10 轨道, dis_win基于E_F=9.3772 |
| 5-akw | Akw_LaNiO2_csc.png | MU_TB=9.3772, Γ-X-M-Γ-Z-R-A-Z |

**对比 one-shot 关键参数**：one-shot MU_TB=9.3839, MU_DMFT=-0.2065, Σ=5×1×1；csc MU_TB=9.3772, MU_DMFT=-0.01226, Σ=5×5(取对角)。

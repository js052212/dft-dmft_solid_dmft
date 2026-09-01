# CSC（电荷自洽）DFT+DMFT 标准操作流程模板（VASP + solid_dmft）

> 用途：作为计算其它材料时从头到尾的标准模板。
> 本模板在 LaNiO₂（node30, `3-AM-DMFT/2-LNO`）上实测跑通，可复用。

## 全局前置条件
- 服务器：node30
- conda 环境：
  - `vasp_build`：跑/编译 VASP
  - `triqs_unstable`：跑 solid_dmft / 后处理
- VASP：6.5.0，**带 Wannier90 接口**（原生无接口时需重编，见 4-wannier）
- solid 环境需补装：`triqs_maxent`（匹配 triqs 版本）、`scikit-image`
- 缺的库按需：`pip install triqs_maxent==3.3.0 scikit-image`

## 目录结构（每个材料一套）
```
<材料>/  （如 2-LNO/）
├── 1-scf/            初始自洽
├── 2-dmft-csc/       CSC DMFT
├── 3-postproc/       多迭代平均 + MaxEnt
├── 4-wannier-10orb/  Wannier 10 轨道
└── 5-akw/            A(k,ω) 绘图
```

## 作业脚本对应表（统一风格：每步一个作业）
| 步骤 | 作业脚本 | 程序 | 风格 |
|---|---|---|---|
| 1-scf | `sub.pbs` | vasp_std | 一个作业 |
| 2-dmft-csc | `run.pbs` | solid_dmft | 一个作业 |
| 3-postproc | `run.pbs` | python3 1-sigma_postprocess.py | 一个作业 |
| 4a-wannier导出 | `run_export.pbs` | vasp_std_wannier | 同 1-scf 风格 |
| 4b-wannier90 | `run_w90.pbs` | wannier90.x | 同 3-postproc 风格 |
| 5-akw | `run.pbs` | python3 plot_akw | 一个作业 |

> **4a 和 4b 是人工衔接**：先提交 4a（等跑出 .mmn/.amn/.eig），再提交 4b（wannier90.x 解纠缠）。

---

## 1-scf：初始自洽

**文件**：POSCAR、POTCAR、INCAR、KPOINTS、sub.pbs

**POSCAR**：放结构优化(relax)好的结构。**注意：`Direct` 前不要留空行**（旧版 triqs plo 解析器会因空行 IndexError）。
```
<材料>
1.0
a 0 0
0 b 0
0 0 c
La Ni O        # 元素
1 1 2          # 原子数
Direct
0.5 0.5 0.5
0.0 0.0 0.0
0.5 0.0 0.0
0.0 0.5 0.0
```

**INCAR（非磁金属）**：
```
ISTART = 0
ICHARG = 2
ENCUT = 600
PREC = Accurate
EDIFF = 1E-6
NELM = 120
ISPIN = 1        # 非磁
ISMEAR = 1       # 金属
SIGMA = 0.05
LREAL = F
LASPH = T
ISYM = -1
ALGO = Normal
LORBIT = 14
LMAXMIX = 6
LCHARG = T
LWAVE  = T
NCORE = 8
```

**KPOINTS**：均匀 Gamma 网格（如 `12 12 12`）

**sub.pbs**：`vasp_build` 环境，`mpirun -np 32 vasp_std`

**得到**：CHGCAR(给csc)、OUTCAR（记 E-fermi，供 Wannier 用）

**验收**：SCF 收敛、记录 E_F

---

## 2-dmft-csc：CSC DMFT（核心，坑最多）

**文件**：INCAR、plo.cfg、dmft_config.toml、KPOINTS、POSCAR/POTCAR(link)、run.pbs

**链接**：`ln -sf ../1-scf/POSCAR ../1-scf/POTCAR`；提交前 `cp ../1-scf/CHGCAR CHGCAR`

**INCAR（csc 版核心）**：
```
ICHARG = 5        # 电荷来自外部(GAMMA)，DMFT 反馈
ISYM = -1         # VASP-PLO CSC 不能用对称
NELM = 80
NELMIN = 10
NELMDL = -1
LSYNCH5 = .TRUE.
LWRTPRJ = .TRUE.
IMIX = 1
AMIX = 0.02
BMIX = 0.001
ALGO = Normal
ISMEAR = -5
SIGMA = 0.05
LMAXMIX = 6
LORBIT = 14
NEDOS = 3001
EMIN = -8
EMAX = 8
LOCPROJ = 2 : d   # 第N个离子的d轨道（改原子编号）
NBANDS = 30
NCORE = 1
```

**plo.cfg**：
```
[General]
DOSMESH = -10 10 3001
efermi = <E_F>     # 手动给费米能（关键字是 efermi 不是 fermi！）

[Group 1]
SHELLS = 1
NORMALIZE = True
EWINDOW = -8 4     # 投影窗(相对E_F)

[Shell 1]
LSHELL = 2
IONS = <原子编号>
```

**dmft_config.toml**：
```
[general]
seedname = "vasp"
set_rot = "hloc"        # ⚠️ 会导致谱波浪；干净优先则去掉或 enforce_off_diag=true
csc = true
h_int_type = "density_density"
U = 5.0
J = 0.8
beta = 40
prec_mu = 0.001
mu_initial_guess = <E_F>
g0_mix = 0.3
n_iter_dmft_first = 6
n_iter_dmft_per = 2
n_iter_dmft = 5          # 验证5；收敛需加大
h5_save_freq = 2
dc = true
dc_type = 0              # sFLL
dc_dmft = true           # 必须

[solver]
type = "cthyb"
n_l = 40
legendre_fit = true
length_cycle = 200
n_warmup_cycles = 1e4
n_cycles_tot = 1e6       # 验证1e6；收敛1e7+
imag_threshold = 0.05    # ⚠️ 默认1e-3会崩，要放宽

[dft]
n_iter = 4
n_cores = 32
dft_code = "vasp"
dft_exec = "/home/jsguo/bin/vasp.6.5.0/bin/vasp_std"
mpi_env = "default"
plo_cfg = "plo.cfg"
projector_type = "plo"
```

**run.pbs**：`triqs_unstable` 环境，`mpirun -np 32 solid_dmft dmft_config.toml`

**机制**：VASP 常驻后台，DMFT 每 N 步通过 vaspgamma.h5 喂电荷给 VASP，VASP 补跑再回 DMFT，交替循环。

**得到**：vasp.h5（DMFT 结果）

**必踩两坑**：
1. `imag_threshold` 放宽到 0.05（csc 首轮 G0 虚部大）
2. 非磁 csc 的 dft_tools bug：`sumk_dft.py` 访问 `band_window[1]` 越界 → 改库为 `len(band_window)>1` 判断

---

## 3-postproc：多迭代平均 + MaxEnt

**文件**：`1-sigma_postprocess.py`、run.pbs

**脚本1 参数**：
```
ITERATIONS_FOR_AVG = ['it_2','it_4','last_iter']   # 按 h5_save_freq
N_INEQUIV_SHELLS = 1
MAXENT_ERROR = 0.02
OMEGA_MIN, OMEGA_MAX = -10, 5
...
```
> block 自适版本可自动处理各种 block 结构。

**步骤**：① 平均 Σ(iω) ② MaxEnt: Σ(iω)→Σ(ω) ③ 因果性修正（插值非删除）

**run.pbs**：`triqs_unstable` 环境，`mpirun -np 32 python3 1-sigma_postprocess.py`，自动 `cp ../2-dmft-csc/vasp.h5 .`

**得到**：vasp_avg.h5（Sigma_maxent_0、DC_pot、chemical_potential_post）

**记录**：chemical_potential_post = MU_DMFT（供绘图）

**坑**：缺 triqs_maxent → pip 装匹配版本

---

## 4-wannier-10orb：Wannier 10 轨道

**文件**：INCAR、KPOINTS、POSCAR/POTCAR/CHGCAR、run_export.pbs、run_w90.pbs

**KPOINTS**：`12 12 12` Gamma

**INCAR（Wannier 接口版）**：
```
ISTART = 1
ICHARG = 11
ENCUT = 600
ISPIN = 1
ISMEAR = 0
SIGMA = 0.05
ISYM = -1
LMAXMIX = 6
NBANDS = 40
NUM_WANN = 10
LWANNIER90 = .TRUE.
LWANNIER90_RUN = .FALSE.
LWRITE_MMN_AMN = .TRUE.
NCORE = 1          # Wannier接口要求NCORE=1

WANNIER90_WIN = "
num_wann = 10
dis_win_min = <E_F - 8>
dis_win_max = <E_F + 10>
dis_froz_min = <E_F - 1>
dis_froz_max = <E_F + 1>
dis_num_iter = 2000
dis_conv_tol = 1.0d-8
dis_conv_window = 10
num_iter = 500
conv_tol = 1.0d-10
begin projections
f=<Ni坐标>:dxy,dyz,dz2,dxz,dx2-y2
f=<La坐标>:dxy,dyz,dz2,dxz,dx2-y2
end projections
guiding_centres = true
write_hr = true
"
```
> 能量窗用绝对能量 = E_F + 相对值；投影坐标来自 POSCAR。

**run_export.pbs（4a，VASP导出，风格同 1-scf）**：
```
#PBS -N <X>-w90-export
#PBS -l nodes=1:ppn=32
source .../activate vasp_build
mpirun -np 32 /home/jsguo/bin/vasp.6.5.0/bin/vasp_std_wannier
```

**run_w90.pbs（4b，解纠缠，风格同 3-postproc）**：
```
#PBS -N <X>-w90-min
#PBS -l nodes=1:ppn=1
source .../activate triqs_unstable
/home/jsguo/.conda/envs/bin/wannier90.x wannier90
```

**得到**：wannier90_hr.dat（核心）、.mmn/.amn/.eig/.chk/.wout

**验收**：hr.dat 第2行=轨道数(10)；解纠缠 `criteria satisfied`；spread 合理(~20.85)

**坑**：node30 原生 VASP 无 Wannier 接口 → 需重编（makefile.include 加 `-DVASP2WANNIER90` + `-lwannier`，libwannier.a 从 node10 拷或重编 Wannier90）

---

## 5-akw：绘图谱函数

**文件**：plot_akw.py、run.pbs

**plot 参数**：
```
MU_TB = <wannier VASP的E_F>
MU_DMFT = <chemical_potential_post>
N_ORB = 10
N_CORR_ORB = 5
DMFT_TO_W90_IDX = [0,1,2,3,4]
DFT_BAND_MODE = "all"
ETA_BROADENING = 0.06
```

**Sigma 嵌入（关键）**：从 Sigma_maxent_0 取 **5×5 的对角元** → 嵌入 Wannier 前5个(Ni-d)；其它轨道自能为0；DC 已扣。
> ⚠️ 用对角元（非对角求逆数值病态→A 出现±200尖峰→纯白条）

**格林函数**：`G = [(ω+μ)I - H(k) - Σ]⁻¹`，`A = -1/π · Im Tr G`，`with_sigma=..., add_mu_tb=True`

**run.pbs**：`triqs_unstable` 环境，`python3 plot_akw.py`

**得到**：Akw_<X>.png

**坑**：缺 scikit-image → pip 装

---

## 关键参数对照（LaNiO₂ 示例）
| 参数 | one-shot | csc |
|---|---|---|
| MU_TB | 9.3839 | 9.3772 |
| MU_DMFT | -0.2065 | -0.01226 |
| Σ 结构 | 5×1×1 | 5×5（取对角嵌入） |
| VASP 角色 | 跑完即退 | 常驻+电荷反馈 |
| dc_dmft | 可 false | 必须 true |

## 换材料时的改动点
1. POSCAR/POTCAR：目标结构
2. LOCPROJ / plo.cfg：投影原子编号、轨道
3. U/J、相关子空间：按材料
4. E_F：每阶段后更新（plo.cfg efermi + Wannier 窗 + MU_TB）
5. mu_initial_guess：用当前 E_F

## 本流程已知问题（可选优化）
- 用 `set_rot="hloc"` + 5×5 矩阵 MaxEnt 时，穿费米带附近可能出现"波浪/波纹"伪影；根治：重跑 csc 去掉 set_rot 或 `enforce_off_diag=true`，使 block 拆成 5×1×1。

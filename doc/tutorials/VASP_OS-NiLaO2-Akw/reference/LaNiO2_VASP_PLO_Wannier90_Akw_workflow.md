# LaNiO₂：VASP/PLO one-shot DMFT + 独立 Wannier90 \(H(\mathbf{k})\) 绘制 \(A(\mathbf{k},\omega)\)

## 精确工作流、Wannier90 操作细节与验收标准

本文档总结已经实际跑通的路线：

```text
VASP SCF
  ├─→ VASP/LOCPROJ → PLO → solid_dmft one-shot → Σ(iω)
  │                                      └─→ MaxEnt → Σ(ω)
  │
  └─→ VASP fixed-charge NSCF → Wannier90
                         └─→ 10 轨道 wannier90_hr.dat → H(k)

Σ(ω) + H(k) + μ + DC
        ↓
solid_dmft.postprocessing.plot_correlated_bands
        ↓
A(k,ω), Γ-X-M-Γ-Z-R-A-Z
```

这是一条“双轨后处理”路线：

- DMFT 的物理自能来自 VASP/PLO；
- 动量依赖哈密顿量来自另一套独立构造的 Wannier90 模型；
- 两者在实频后处理阶段拼接。

本文档以已经使用的环境和文件格式为准：

- VASP 6.5.0（另行编译了 Wannier90 接口版本）；
- Wannier90 独立可执行文件；
- TRIQS/DFTTools；
- solid_dmft 3.3.5；
- 非磁性 one-shot DMFT；
- LaNiO₂；
- 目标图：论文 arXiv:2606.00223 Figure 4 上排中间的 undoped \(A(\mathbf{k},\omega)\)。

> **当前接口说明**：最新 TRIQS/DFTTools 文档已经描述了
> `KPOINTS_OPT + LOCPROJ_OPT + vaspout.h5` 的原生 VASP 路径投影流程。若所用
> VASP/DFTTools 组合能够完整写出并转换 `LOCPROJ_OPT`，优先采用该原生路线可以避免
> PLO 与独立 MLWF 之间的基底拼接问题。本文仍详细记录独立 Wannier90 路线，因为这是
> 本项目已经实际跑通、并且需要总结的工作流。

---

## 一页版执行顺序

1. 用固定结构、POTCAR、ENCUT 和 PBE 完成非磁性 VASP SCF，保存 `CHGCAR`。
2. 从同一 `CHGCAR` 做 PLO fixed-charge NSCF，生成 `LOCPROJ`。
3. 用正确的 `plo.cfg` 转换；根据目标选择 Ni-\(d\) 五轨道或论文 Ni-\(e_g\) 两轨道，
   两者不能混用。
4. 运行 solid_dmft one-shot，确认占据、化学势和低频自能收敛。
5. 只平均已收敛的后期迭代，用 `maxent_sigma` 得到 `Sigma_maxent_0`。
6. 新建干净 `4-wannier-10orb` 目录，使用均匀 12×12×12 k 网格。
7. 用 Wannier 接口版 VASP 导出 40 条带、10 个投影的 `.mmn/.amn/.eig`。
8. 以当前 `E_F` 将相对能窗换算为 `.eig` 使用的绝对能量；检查每个 k 点的外窗和冻结窗
   带数。
9. 独立运行 `wannier90.x wannier90`，要求解纠缠明确收敛。
10. 检查 `wannier90_hr.dat` 为 10 轨道、WF 中心和 spread 合理，并让 Wannier 插值带
    在冻结窗内与 VASP 带重合。
11. 检查 HDF5 自能 block 的真实维度，建立 PLO→W90 映射；若基底存在旋转，使用
    \(R\Sigma R^\dagger\)。
12. 扣除 DC 一次、统一 `MU_TB/MU_DMFT` 能量参考，再调用
    `pcb.get_dmft_bands()` 绘制 Γ-X-M-Γ-Z-R-A-Z 的 \(A(\mathbf{k},\omega)\)。

## 0. 必须先分清的三套“轨道数”

这是整个工作流中最容易混淆、也最重要的地方。

### 0.1 论文的相关子空间

论文使用：

\[
\mathrm{Ni}\text{-}e_g=
\left(d_{z^2},d_{x^2-y^2}\right),
\]

即两个相关轨道。论文还给出了：

- PBE；
- \(8\times8\times10\) k 网格；
- 投影窗 \(E_F\pm10\) eV；
- \(U=5\) eV；
- \(J=0.8\) eV；
- Kanamori 相互作用；
- Held 双计数；
- \(\beta=40\ {\rm eV}^{-1}\)；
- 40 阶 Legendre 表示。

### 0.2 已经实际跑通的绘图脚本

最终成功脚本使用：

- 10 轨道 Wannier90 哈密顿量：

\[
5\times\mathrm{Ni}\text{-}d+
5\times\mathrm{La}\text{-}d;
\]

- 从 HDF5 中读取 5 个 Ni-\(d\) 自能分量；
- 将这 5 个自能放入 10×10 自能矩阵的前 5 个对角元；
- La-\(d\) 的 5×5 自能块置零。

因此，“已经跑通的版本”在实现上是：

\[
\Sigma_{10\times10}(\omega)=
\begin{pmatrix}
\Sigma_{\mathrm{Ni},5d}(\omega)&0\\
0&0_{\mathrm{La},5d}
\end{pmatrix}.
\]

### 0.3 若严格按论文的两轨道 \(e_g\) 自能

仍可保留 10 轨道 Wannier 哈密顿量，以保留 La 电子口袋；但只应嵌入两个自能：

```text
Wannier 顺序：
0 Ni_dxy
1 Ni_dyz
2 Ni_dz2
3 Ni_dxz
4 Ni_dx2-y2
5 La_dxy
6 La_dyz
7 La_dz2
8 La_dxz
9 La_dx2-y2
```

论文 \(e_g\) 自能应映射到：

```text
Ni_dz2     → Wannier index 2
Ni_dx2-y2  → Wannier index 4
```

其他 Ni-\(d\) 和 La-\(d\) 轨道的自能均置零。

不能把“两轨道 PLO 配置”和“五轨道绘图循环”混用。运行前必须以 HDF5 中的实际 block
名称与维度为准。

---

## 1. 数学上实际计算的对象

Wannier90 提供实空间跃迁矩阵 \(H(\mathbf{R})\)，傅里叶变换得到：

\[
H(\mathbf{k})=\sum_{\mathbf{R}}
H(\mathbf{R})e^{i\mathbf{k}\cdot\mathbf{R}}.
\]

实频格林函数为：

\[
G(\mathbf{k},\omega)=
\left[
(\omega+\mu+i\eta)I
-H(\mathbf{k})
-\Sigma(\omega)
\right]^{-1}.
\]

总谱函数为：

\[
A(\mathbf{k},\omega)=
-\frac{1}{\pi}
\operatorname{Im}\operatorname{Tr}G(\mathbf{k},\omega).
\]

其中：

- \(H(\mathbf{k})\) 来自 `wannier90_hr.dat`；
- \(\Sigma(\omega)\) 来自 `Sigma_maxent_0`；
- \(\eta\) 是额外绘图展宽；
- \(\mu\) 必须在与 \(H(\mathbf{k})\) 相同的能量零点下使用。

本次实际参数为：

```python
MU_TB = 9.3839
MU_DMFT = -0.20653052115213022
add_mu_tb = True
```

因此 `get_dmft_bands()` 中使用的总化学势参考为：

\[
\mu_{\mathrm{used}}=\mu_{\mathrm{TB}}+\mu_{\mathrm{DMFT}}.
\]

这些数值只能用于产生它们的那一套 VASP/Wannier/DMFT 数据。换材料、换结构、换赝势、
换电荷数或重跑 SCF 后都必须重新读取。

---

## 2. 形式严格性：PLO 基底与 Wannier 基底

这条双轨路线能够工作，但必须明确其基底假设。

PLO 自能严格属于 PLO 局域基底，Wannier 哈密顿量属于 MLWF 基底。形式上应存在一个局域
变换矩阵 \(R\)：

\[
\Sigma_{\mathrm{W90}}(\omega)
=R\,\Sigma_{\mathrm{PLO}}(\omega)R^\dagger.
\]

当前成功脚本使用的是：

\[
R\approx I,
\]

即假设：

1. PLO 和 MLWF 使用相同的实立方谐函数；
2. 两边轨道顺序完全一致；
3. 局域轴方向一致；
4. Ni-\(d\) 子空间内部没有显著的额外旋转；
5. 非对角自能可以忽略或已经在所选基底中对角化。

这是一种可用且已得到合理谱图的近似，但如果要把工作流称为“形式上严格”，还应计算
PLO 与 MLWF 之间的重叠/旋转矩阵，或至少完成以下验证：

- 比较两套局域哈密顿量；
- 比较轨道占据；
- 比较投影 DOS；
- 检查每个 MLWF 的中心、轨道形状与局域轴；
- 检查自能变换后是否出现不能忽略的非对角元。

因此本文档区分：

- **操作上已经跑通**：按相同轨道名称和顺序直接嵌入；
- **形式上完全严格**：额外求出 \(R\) 并执行 \(R\Sigma R^\dagger\)。

---

## 3. 推荐目录结构

```text
LaNiO2_Akw/
├── 1-scf/
│   ├── INCAR
│   ├── POSCAR
│   ├── POTCAR
│   ├── KPOINTS
│   ├── CHGCAR
│   ├── OUTCAR
│   └── vasprun.xml
│
├── 2-plo-nscf/
│   ├── INCAR
│   ├── CHGCAR
│   ├── LOCPROJ
│   ├── OUTCAR
│   └── vasprun.xml
│
├── 3-dmft/
│   ├── plo.cfg
│   ├── dmft_config.toml
│   └── vasp/
│       ├── vasp.h5
│       └── vasp_avg.h5
│
├── 4-wannier-10orb/
│   ├── INCAR
│   ├── POSCAR
│   ├── POTCAR
│   ├── KPOINTS
│   ├── CHGCAR
│   ├── wannier90.win
│   ├── wannier90.mmn
│   ├── wannier90.amn
│   ├── wannier90.eig
│   ├── wannier90.wout
│   ├── wannier90_hr.dat
│   └── wannier90_band.*
│
└── 5-akw/
    ├── plot_akw_general_recovered.py
    └── Akw_LaNiO2_filtered.png
```

不同阶段必须放在不同目录中。尤其不要把旧的 5 轨道和新的 10 轨道
`wannier90.*` 文件混在一起。

---

## 4. 阶段 A：VASP SCF

### 4.1 结构

论文的未掺杂 LaNiO₂ 晶格常数：

```text
a = 3.97 Å
c = 3.37 Å
```

必须固定一套原点和原子坐标约定，并在 SCF、PLO 和 Wannier 三条路径中使用完全相同的
POSCAR。

注意：归档记录中最终成功 Wannier 计算采用的坐标约定是：

```text
Ni 约在 (0,0,0)
La 约在 (0.5,0.5,0.5)
```

本地早期模板的 POSCAR 使用了另一种原点：

```text
La = (0,0,0)
Ni = (0.5,0.5,0)
```

两者不能混用。Wannier `projections` 中的坐标必须来自实际运行的 POSCAR。

### 4.2 SCF k 网格

论文使用 \(8\times8\times10\)。用于复现论文 DMFT 参数时可采用：

```text
LaNiO2 8x8x10 Gamma centered
0
Gamma
8 8 10
0 0 0
```

### 4.3 SCF 基本要求

- PBE；
- 非磁性：`ISPIN = 1`；
- `LREAL = .FALSE.`；
- `LASPH = .TRUE.`；
- 足够严格的 `EDIFF`，建议 \(10^{-8}\)；
- 保存 `CHGCAR`；
- `ENCUT`、POTCAR、结构必须记录；
- 检查 SCF 真正收敛，而不是只检查作业是否正常退出。

建议检查：

```bash
grep "reached required accuracy" OUTCAR
grep "E-fermi" OUTCAR | tail -n 1
tail -n 10 OSZICAR
ls -lh CHGCAR vasprun.xml OUTCAR
```

---

## 5. 阶段 B：VASP/PLO 和 one-shot DMFT

### 5.1 PLO fixed-charge NSCF

从 SCF 复制：

```bash
cp ../1-scf/POSCAR .
cp ../1-scf/POTCAR .
cp ../1-scf/KPOINTS .
cp ../1-scf/CHGCAR .
```

若可靠复用同一套 `WAVECAR`，可用 `ISTART=1`；否则使用 `ISTART=0` 重新对固定电荷密度
对角化。核心设置：

```ini
ISTART = 0
ICHARG = 11

PREC   = Accurate
ENCUT  = <与 SCF 一致>
EDIFF  = 1E-8
LREAL  = .FALSE.
LASPH  = .TRUE.

ISPIN  = 1
ISYM   = -1

LORBIT  = 14
LOCPROJ = <实际 Ni 原子编号> : d : Ni 1

EMIN = -10.0
EMAX =  10.0
```

对文本型旧 VASP PLO 接口，`ISYM=-1` 是必要的安全设置。VASP 6.5.0 HDF5 接口已改善
对称性支持，但同一项目中不要混用两种转换方式。

检查：

```bash
ls -lh LOCPROJ OUTCAR vasprun.xml
grep "E-fermi" OUTCAR | tail -n 1
```

### 5.2 `plo.cfg`：五轨道版与两轨道版

#### 已跑通脚本对应的五轨道顺序

VASP/TRIQS 常用 d 轨道顺序为：

```text
0 dxy
1 dyz
2 dz2
3 dxz
4 dx2-y2
```

若 DMFT HDF5 中确实存在 5 个标量自能 block，应使用完整 5×5 恒等变换：

```ini
TRANSFORM = 1 0 0 0 0
            0 1 0 0 0
            0 0 1 0 0
            0 0 0 1 0
            0 0 0 0 1
```

#### 严格论文 \(e_g\) 两轨道版

```ini
TRANSFORM = 0 0 1 0 0
            0 0 0 0 1
```

对应：

```text
row 0 → dz2
row 1 → dx2-y2
```

投影归一化窗按论文为：

```ini
EWINDOW = -10 10
```

PLO 的 `EWINDOW` 相对于 \(E_F\)，而本次 `wannier90.eig` 中使用的是 VASP 绝对本征值。
这两个能量约定不要混淆。

### 5.3 DMFT 参数

严格论文参数：

```toml
beta = 40
h_int_type = "kanamori"
U = 5.0
J = 0.8
dc_type = 1
magnetic = false
```

Held 双计数公式中的轨道数 \(d\) 必须与实际相关子空间一致：

- 两轨道 \(e_g\)：\(d=2\)；
- 五轨道 \(d\)：\(d=5\)。

五轨道和两轨道计算的双计数、占据和自能并不等价。

CT-HYB 的 `n_cycles_tot`、warmup、tail fit 和 DMFT 迭代数属于数值收敛参数；论文未全部
公开，不能把任意一组 Monte Carlo 参数称为论文原值。必须通过误差条和收敛测试确定。

### 5.4 DMFT 验收

至少检查：

- 化学势是否稳定；
- 总占据与各轨道占据是否稳定；
- 最近数次迭代的 \(\Sigma(i\omega_n)\) 是否重合；
- 低频自能噪声是否可接受；
- `Im Σ(iω_n)` 的低频行为是否物理；
- CT-HYB 平均符号、测量统计和误差；
- 结果是否对增加 Monte Carlo 步数稳定。

只有已收敛的后期迭代可以用于平均。不能为了减小噪声而平均仍在漂移的迭代。

---

## 6. 阶段 C：自能解析延拓

使用：

```python
from solid_dmft.postprocessing import maxent_sigma
```

目标是在 HDF5 中获得：

```text
DMFT_results/it_avg/Sigma_maxent_0
DMFT_results/it_avg/DC_pot
DMFT_results/it_avg/chemical_potential_post
```

### 6.1 正确做法

1. 选择已收敛的后期迭代；
2. 平均 \(\Sigma(i\omega_n)\) 以降低随机噪声；
3. 对自能而不是格林函数调用 `maxent_sigma`；
4. 延拓能窗必须覆盖最终绘图能窗，并留出足够边界；
5. 保存原始 HDF5，不直接覆盖；
6. 检查延拓后的因果性：

\[
\operatorname{Im}\Sigma(\omega)\le 0.
\]

### 6.2 因果性处理

少量由数值插值造成的孤立正值点可以标记并检查；若大范围
\(\operatorname{Im}\Sigma(\omega)>0\)，应重新调整 MaxEnt、输入噪声和延拓参数，
而不应靠大面积裁剪或插值把问题隐藏起来。

后处理“修正”只能用于可视化的小型数值瑕疵，不能替代物理可靠的解析延拓。

---

## 7. 阶段 D：Wannier90 详细工作流

## 7.1 为什么最终选择 10 轨道

早期只用 5 个 Ni-\(d\) Wannier 轨道会丢失 La-\(d\) 电子口袋，尤其是论文中 Γ 和 A
附近的重要低能结构。

最终选择：

\[
\mathrm{Ni}\text{-}3d(5)+\mathrm{La}\text{-}5d(5)=10
\]

个 Wannier 函数。

这样：

- Ni-\(d\) 子空间可接收 DMFT 自能；
- La-\(d\) 子空间保留 DFT 色散；
- Γ、A 附近的电子口袋仍存在；
- `wannier90_hr.dat` 可以在任意 k 路径上插值。

## 7.2 一次性编译

### 7.2.1 编译 Wannier90

Wannier90 库与 VASP 应使用兼容的 Fortran 编译器、MPI 和数学库。

```bash
mpif90 -show
gfortran --version
```

示意：

```bash
cd /path/to/wannier90
cp config/make.inc.gfort make.inc

# 检查并修改 F90、COMMS、LIBS 等
make clean
make lib
make wannier
```

产物：

```text
libwannier.a
wannier90.x
```

如果编译的是 MPI 版 Wannier90，可用 MPI 运行；本次实际使用的是串行
`wannier90.x`，因此 Wannier 最小化阶段申请单核即可。

### 7.2.2 编译 VASP Wannier 接口版

VASP 必须启用：

```make
CPP_OPTIONS += -DVASP2WANNIER90
```

或与所用 VASP/Wannier90 版本匹配的：

```make
CPP_OPTIONS += -DVASP2WANNIER90v2
```

并链接 Wannier90 库，例如：

```make
WANNIER90_ROOT ?= /path/to/wannier90
LLIBS += -L$(WANNIER90_ROOT) -lwannier
```

重新编译：

```bash
make clean
make std
```

检查：

```bash
/path/to/vasp_wannier/bin/vasp_std
```

若 `LWANNIER90=.TRUE.` 时提示接口不存在，说明编译宏或链接没有成功。

## 7.3 建立完全干净的 10 轨道目录

```bash
cd /home/jsguo/2-AM-DMFT/4-1-LNO-wan90
mkdir 4-wannier-10orb
cd 4-wannier-10orb

cp ../1-scf/CHGCAR .
cp ../1-scf/POSCAR .
cp ../1-scf/POTCAR .
```

如果采用 `ISTART=0`，不需要复制 `WAVECAR`。

第一次计算时不要复制任何旧文件：

```text
wannier90.win
wannier90.mmn
wannier90.amn
wannier90.eig
wannier90.chk
wannier90_hr.dat
wannier90.wout
```

最安全的规则是：每次改变轨道数、投影原子、NBANDS 或 k 网格，都新建目录。

## 7.4 Wannier k 网格

实际成功使用：

```text
Automatic mesh
0
Gamma
12 12 12
0 0 0
```

这里必须是均匀三维网格，不是高对称线。高对称路径由 Wannier 插值阶段另外指定。

Wannier k 网格不必等于 DMFT 的 \(8\times8\times10\)，但必须对 MLWF 和插值收敛。
应至少测试：

```text
10×10×10
12×12×12
14×14×14
```

并比较低能插值能带。

## 7.5 能量参考：本次最关键的诊断

本次：

```text
E_F = 9.3839 eV
```

`wannier90.eig` 使用 VASP 的绝对本征值，而不是自动减去费米能。

因此相对 \(E_F\) 的能窗必须转换为绝对能量：

\[
E_{\rm abs}=E_F+E_{\rm rel}.
\]

例如：

```text
相对外窗 [-8,+4] eV
→ 绝对外窗 [1.3839,13.3839] eV

相对冻结窗 [-1,+1] eV
→ 绝对冻结窗 [8.3839,10.3839] eV
```

最终成功的解纠缠计算把外窗上限扩展到了：

```text
dis_win_max = 19.3839  # E_F + 10 eV
```

最终成功参数为：

```text
dis_win_min  = 1.3839   # E_F - 8 eV
dis_win_max  = 19.3839  # E_F + 10 eV
dis_froz_min = 8.3839   # E_F - 1 eV
dis_froz_max = 10.3839  # E_F + 1 eV
```

这些数值是本次数据专用值。换计算后应先运行：

```bash
grep "E-fermi" OUTCAR | tail -n 1
head -n 10 wannier90.eig
```

并重新换算。

## 7.6 外窗、冻结窗和 NBANDS 的硬约束

在每一个 k 点：

1. 外窗内必须至少有 `num_wann` 条带；
2. 冻结窗内的带数不能超过 `num_wann`；
3. `NUM_WANN <= NBANDS`；
4. `NBANDS` 必须覆盖外窗中希望纳入的最高状态。

本次使用：

```text
NBANDS = 40
NUM_WANN = 10
```

初始外窗上限为 \(E_F+4\) eV 时，40 条带有充分余量。后来只修改 Wannier 解纠缠外窗为
\(E_F+10\) eV，并复用原有 40 条带数据。对于新的材料，不能假设 40 条带一定覆盖
\(E_F+10\) eV，必须统计 `.eig`。

建议写一个检查程序，逐 k 点统计：

```text
n_outer(k) = 外窗内能带数
n_frozen(k) = 冻结窗内能带数
```

验收条件：

```text
min_k n_outer(k) >= num_wann
max_k n_frozen(k) <= num_wann
```

## 7.7 VASP Wannier NSCF：本次成功模板

下列模板体现本次最终路线。`ENCUT` 必须与实际 SCF 保持一致；归档指导中使用过 600 eV，
但如果 SCF 是 500 eV，就应统一为 500 eV，或先以 600 eV 重跑 SCF。

```ini
SYSTEM = LaNiO2 Wannier Ni-d + La-d

# 固定 CHGCAR 的非自洽计算
ISTART = 0
ICHARG = 11

ENCUT = <与 SCF 相同>
PREC  = Accurate
EDIFF = 1E-8
NELM  = 80
ALGO  = Normal

ISPIN = 1
ISMEAR = 0
SIGMA  = 0.05

ISYM = -1
LMAXMIX = 6

NBANDS = 40

LCHARG = .FALSE.
LWAVE  = .FALSE.

KPAR  = 1
NCORE = 8

# Wannier90 接口
NUM_WANN = 10
LWANNIER90 = .TRUE.
LWANNIER90_RUN = .FALSE.
LWRITE_MMN_AMN = .TRUE.

WANNIER90_WIN = "
# 本次 .eig 使用 VASP 绝对能量

dis_win_min = 1.3839
dis_win_max = 19.3839

dis_froz_min = 8.3839
dis_froz_max = 10.3839

dis_num_iter = 2000
dis_conv_tol = 1.0d-8
dis_conv_window = 10

num_iter = 500
conv_tol = 1.0d-10
conv_window = 10

begin projections
f=0.0,0.0,0.0:dxy,dyz,dz2,dxz,dx2-y2
f=0.5,0.5,0.5:dxy,dyz,dz2,dxz,dx2-y2
end projections

guiding_centres = true
write_hr = true

fermi_energy = 9.3839
bands_plot = true
bands_num_points = 100

begin kpoint_path
G 0.0 0.0 0.0 X 0.5 0.0 0.0
X 0.5 0.0 0.0 M 0.5 0.5 0.0
M 0.5 0.5 0.0 G 0.0 0.0 0.0
G 0.0 0.0 0.0 Z 0.0 0.0 0.5
Z 0.0 0.0 0.5 R 0.5 0.0 0.5
R 0.5 0.0 0.5 A 0.5 0.5 0.5
A 0.5 0.5 0.5 Z 0.0 0.0 0.5
end kpoint_path
"
```

### 投影坐标的严格规则

上面两行坐标只适用于实际成功计算所用的原点约定。

如果 POSCAR 是：

```text
La = (0,0,0)
Ni = (0.5,0.5,0)
```

则必须改为：

```ini
begin projections
f=0.5,0.5,0.0:dxy,dyz,dz2,dxz,dx2-y2  # Ni
f=0.0,0.0,0.0:dxy,dyz,dz2,dxz,dx2-y2  # La
end projections
```

不要凭元素名称或旧文件猜坐标。先读实际 POSCAR。

### 投影语法

同一位置的多个轨道用逗号：

```text
dxy,dyz,dz2,dxz,dx2-y2
```

不要用分号。投影总数必须与 `NUM_WANN=10` 一致。

## 7.8 为什么不要复用旧 `wannier90.win`

VASP 在 `LWANNIER90=.TRUE.` 时会：

- 生成或补充 `unit_cell_cart`；
- 生成或补充 `atoms_cart`；
- 写入 `mp_grid`；
- 写入完整 `kpoints`；
- 设置 `num_bands`；
- 按 `NUM_WANN` 设置 `num_wann`；
- 生成 `.mmn/.amn/.eig`。

若旧 `.win` 已存在，VASP 只补缺失块，不会验证已有 `kpoints`、晶格或原子是否与当前 VASP
计算一致。因此最安全的方式是在确认目录后先备份旧文件：

```bash
mkdir -p old_w90_input
mv wannier90.win old_w90_input/
```

然后让当前 VASP 运行重新生成。移动前必须确认处于新建的
`4-wannier-10orb` 目录并已备份旧结果。

不要在 `WANNIER90_WIN` 中写：

```text
dis_win_placeholder
frozen_window_placeholder
```

这些不是注释，而是非法关键字。本次曾因此在 `wannier_setup` 中退出，只生成
`wannier90.win/wout`，没有正确完成 `.mmn/.amn/.eig`。

## 7.9 VASP 作业脚本

MPI 进程数必须与调度器申请一致。示例：

```bash
#!/bin/bash
#PBS -N LNO-w90-export
#PBS -q normal
#PBS -o jsguo.out
#PBS -e jsguo.err
#PBS -l nodes=1:ppn=32
#PBS -l walltime=48:00:00

cd "$PBS_O_WORKDIR" || exit 1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

mpirun -np 32 \
  /home/jsguo/bin/2-VASP/vasp.6.5.0_lic_wannier/bin/vasp_std
```

若 `NCORE=8`，MPI 总进程数应能被 8 整除。

## 7.10 VASP 导出阶段验收

VASP 结束后先不要马上运行 Wannier90：

```bash
grep -iE \
"num_bands|num_wann|dis_win|dis_froz|fermi_energy" \
wannier90.win

grep "E-fermi" OUTCAR | tail -n 1

ls -lh \
wannier90.win \
wannier90.eig \
wannier90.amn \
wannier90.mmn
```

本次曾得到：

```text
wannier90.amn   35 MB
wannier90.eig   3.1 MB
wannier90.mmn   586 MB
wannier90.win   104 KB

num_bands = 40
num_wann = 10
```

文件大小仅作为完整性线索，不是物理验收标准。

还应检查：

```bash
grep -iE "error|unrecognised|aborting" wannier90.wout OUTCAR
```

如果修改以下内容，必须重跑 VASP：

- POSCAR/POTCAR/ENCUT；
- KPOINTS 或 k 网格；
- NBANDS；
- NUM_WANN；
- projections；
- 自旋或 SOC 设置。

如果只修改以下内容，通常可直接复用 `.mmn/.amn/.eig`，只重跑 Wannier90：

- `dis_win_*`（但必须仍在已有本征值覆盖范围内）；
- `dis_froz_*`；
- `dis_num_iter`；
- `dis_conv_tol/window`；
- `num_iter`、`conv_tol/window`；
- `bands_plot` 和高对称路径。

## 7.11 运行独立 Wannier90

因为 VASP 已调用 `wannier_setup` 并生成 `.mmn/.amn/.eig`，本次直接运行：

```bash
~/bin/wannier90_gfortran/wannier90.x wannier90
```

通常不需要再次 `wannier90.x -pp wannier90`。`-pp` 主要用于传统流程中先生成
`.nnkp` 再交给外部 DFT 接口。

本次使用串行 Wannier90，PBS 示例：

```bash
#!/bin/bash
#PBS -N LNO-w90-min
#PBS -q normal
#PBS -o w90.out
#PBS -e w90.err
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00

cd "$PBS_O_WORKDIR" || exit 1
~/bin/wannier90_gfortran/wannier90.x wannier90
```

如果确认编译的是 MPI 版 Wannier90，才使用 `mpirun`。

## 7.12 本次 Wannier 收敛过程

### 第一次不合格结果

曾出现：

```text
Omega D  ≈ 29.59
Omega OD ≈ 60.38
Omega Total ≈ 108.57 Å²
```

La 相关 WF 展宽约 14–24 Å²。这不是合格结果。

### 第二次：展宽改善，但解纠缠未收敛

曾出现：

```text
Omega D = 0
Omega OD ≈ 0.424
Omega Total ≈ 20.920 Å²
```

但日志仍有：

```text
Maximum number of disentanglement iterations reached
Disentanglement convergence criteria not satisfied
```

因此仅看 spread 小还不能验收。

### 最终合格设置

```ini
dis_num_iter = 2000
dis_conv_tol = 1.0d-8
dis_conv_window = 10
```

最终日志出现：

```text
Disentanglement convergence criteria satisfied
```

最终 spread：

```text
WF 1   0.88252185 Å²
WF 2   1.02450153 Å²
WF 3   1.02449661 Å²
WF 4   2.26937509 Å²
WF 5   1.11834551 Å²
WF 6   3.93540191 Å²
WF 7   2.21624293 Å²
WF 8   2.21624369 Å²
WF 9   2.29905687 Å²
WF 10  3.93419116 Å²

Omega I     = 20.496619882 Å²
Omega D     = 0.000000000 Å²
Omega OD    = 0.423757261 Å²
Omega Total = 20.920377143 Å²
```

WF 中心：

```text
WF 1–5  接近 (0,0,0)
WF 6–10 接近 (-1.98,-1.98,-1.685) Å
```

后者与 \((0.5,0.5,0.5)\) 相差晶格平移，属于同一晶格等价位置。

## 7.13 Wannier 最终检查命令

```bash
grep -iE \
"warning|error|not converged|convergence criteria" \
wannier90.wout | tail -n 30

grep -A 16 "Final State" wannier90.wout | tail -n 17

head -n 4 wannier90_hr.dat

ls -lh wannier90_hr.dat wannier90_band.*
```

本次 `wannier90_hr.dat` 头部为：

```text
written on ...
10
2197
...
```

第二行 `10` 是轨道数；若不是 10，不能交给 10 轨道绘图脚本。

## 7.14 最重要的物理验收：插值能带

spread 和收敛日志都通过后，还必须比较：

- VASP 原始能带；
- Wannier90 插值能带；
- 特别是冻结窗 \([E_F-1,E_F+1]\) eV 内的曲线；
- Γ 与 A 附近 La-\(d\) 电子口袋；
- Ni-\(d_{x^2-y^2}\) 穿越 \(E_F\) 的能带；
- 高对称点顺序和坐标。

`bands_plot=true` 会产生：

```text
wannier90_band.dat
wannier90_band.gnu
wannier90_band.kpt
wannier90_band.labelinfo.dat
```

若集群 gnuplot 没有 `pngcairo`，可输出 SVG：

```bash
gnuplot -e \
"set terminal svg size 1600,1000; set output 'wannier90_band.svg'" \
wannier90_band.gnu
```

验收标准不是“曲线大致像”，而是冻结窗内应几乎重合；外窗边缘允许误差逐渐增大。

---

## 8. 阶段 E：轨道映射、自能嵌入和双计数

## 8.1 先检查 HDF5 实际 block

不要根据文件名猜自能维度。运行类似：

```python
from h5 import HDFArchive

with HDFArchive("vasp_avg.h5", "r") as ar:
    sigma = ar["DMFT_results"]["it_avg"]["Sigma_maxent_0"]
    print(sigma.indices)
    for name, gf in sigma:
        print(name, gf.target_shape)
```

可能出现：

```text
up_0 ... up_4，每个 1×1
```

也可能出现：

```text
up_0，一个 2×2 或 5×5 block
```

两者必须使用不同嵌入代码。

## 8.2 已跑通五轨道版

```python
DMFT_TO_W90_IDX = [0, 1, 2, 3, 4]
```

对应：

```text
PLO dxy     → W90 Ni_dxy
PLO dyz     → W90 Ni_dyz
PLO dz2     → W90 Ni_dz2
PLO dxz     → W90 Ni_dxz
PLO dx2-y2  → W90 Ni_dx2-y2
```

## 8.3 论文两轨道版

若 PLO block 顺序为：

```text
0 dz2
1 dx2-y2
```

则：

```python
DMFT_TO_W90_IDX = [2, 4]
N_CORR_ORB = 2
```

绝不能继续循环 5 次。

## 8.4 双计数只能扣一次

本次脚本使用：

```python
sigma_value = sigma_imp0[block_name].data[:, 0, 0] - dc_value
```

已经把 `DC_pot` 从自能中扣除。因此传给 `get_dmft_bands()` 的 `dc` 必须是零矩阵：

```python
zero_dc = np.zeros((N_ORB, N_ORB), dtype=complex)
```

如果传入非零 DC，就会重复扣除。

`DC_pot` 是按自旋保存的矩阵，例如：

```python
dc_imp0["up"][i, i]
```

不能把它误当作 `up_0`、`up_1` 等独立 block。

## 8.5 非对角自能

若 HDF5 中的自能 block 是 2×2 或 5×5，必须把完整矩阵嵌入：

```text
Σ_W90[w90_indices, w90_indices] = R Σ_PLO R†
```

不能只取对角元，否则会丢失轨道混合。

---

## 9. 阶段 F：计算与绘制 \(A(\mathbf{k},\omega)\)

最终恢复版脚本：

```text
plot_akw_general_recovered.py
```

关键配置：

```python
DMFT_H5_PATH = ".../3-dmft/vasp/vasp_avg.h5"
W90_PATH = ".../4-wannier-10orb/"
W90_SEED = "wannier90"

N_ORB = 10
N_CORR_ORB = 5  # 若严格两轨道 e_g，必须改成 2

MU_TB = 9.3839
MU_DMFT = -0.20653052115213022
ITERATION = "it_avg"

BANDS_PATH = [
    ("G", "X"),
    ("X", "M"),
    ("M", "G"),
    ("G", "Z"),
    ("Z", "R"),
    ("R", "A"),
    ("A", "Z"),
]
```

最终绘图设置：

```python
Y_LIM = (-2.0, 1.0)
W_MESH_WINDOW = [-2.0, 1.0]
W_MESH_NPOINTS = 1200
ETA_BROADENING = 0.04

COLORSCHEME_ALATT = "gist_heat_r"
SPECTRAL_GAMMA = 0.60
SPECTRAL_VMAX_PERCENTILE = 99.5

DFT_LINE_COLOR = "black"
DFT_BAND_MODE = "spectral_match"
```

### 绘图参数与物理参数的区别

- `ETA_BROADENING` 会改变谱峰宽度，是数值/可视化展宽；
- `PowerNorm`、gamma、percentile 只改变颜色映射；
- DFT 黑线筛选只改变叠图显示；
- 这些参数不能用来掩盖错误的自能、化学势、轨道映射或 Wannier 模型。

运行：

```bash
python -m py_compile plot_akw_general_recovered.py
python plot_akw_general_recovered.py
```

---

## 10. 最终验收清单

只有以下项目全部通过，才认为流程完整：

### VASP/DMFT

- [ ] SCF 收敛；
- [ ] POSCAR/POTCAR/ENCUT 在各阶段一致；
- [ ] PLO 投影窗与目标模型一致；
- [ ] HDF5 中的相关轨道数已确认；
- [ ] DMFT 化学势、占据和自能已收敛；
- [ ] CT-HYB 统计误差可接受；
- [ ] MaxEnt 自能基本满足因果性。

### Wannier90

- [ ] 使用均匀 12×12×12 网格；
- [ ] `NBANDS=40`、`NUM_WANN=10`；
- [ ] 投影坐标来自实际 POSCAR；
- [ ] `wannier90.mmn/.amn/.eig` 完整；
- [ ] 解纠缠明确显示 `criteria satisfied`；
- [ ] `wannier90_hr.dat` 第二行是 10；
- [ ] WF 中心落在 Ni 和 La 等价位置；
- [ ] spread 与最终成功值量级一致；
- [ ] VASP 与 Wannier 插值能带在冻结窗内重合；
- [ ] Γ/A 的 La 电子口袋存在。

### 拼接

- [ ] PLO 与 W90 轨道顺序一致；
- [ ] 若存在基底旋转，已经处理 \(R\Sigma R^\dagger\)；
- [ ] 自能维度与 `N_CORR_ORB` 一致；
- [ ] DC 只扣一次；
- [ ] `MU_TB` 与当前 `OUTCAR` 一致；
- [ ] `MU_DMFT` 与当前 `vasp_avg.h5` 一致；
- [ ] La 自能置零是明确的物理近似，而不是数组初始化事故。

### 图像

- [ ] 路径是 Γ-X-M-Γ-Z-R-A-Z；
- [ ] 能量零点是最终 DMFT 化学势；
- [ ] Ni-\(d_{x^2-y^2}\) 带显示明显重整化；
- [ ] Γ/A 电子口袋基本跟随 DFT；
- [ ] 改变 \(\eta\) 或色标后主要色散位置不发生非物理跳动。

---

## 11. 常见错误与处理

| 症状 | 原因 | 处理 |
|---|---|---|
| `Unrecognised keyword dis_win_placeholder` | `.win` 中写了占位符 | 删除非法行，重新跑 VASP 导出 |
| 只有 `win/wout`，没有 `mmn/amn/eig` | `wannier_setup` 中途失败 | 查 `wout/OUTCAR`，修复后重跑 VASP |
| `num_wann` 不是 10 | `NUM_WANN` 错或复用了旧文件 | 新建干净目录，设 `NUM_WANN=10` |
| 投影数与 `num_wann` 不符 | projections 少了一行或语法错误 | Ni、La 各列 5 个 d 轨道 |
| 冻结窗没有任何带 | 把绝对本征值误当成相对 \(E_F\) | 用 \(E_{\rm abs}=E_F+E_{\rm rel}\) |
| spread 很小但日志未收敛 | 只看了最终 spread | 增加 `dis_num_iter` 并要求 criteria satisfied |
| La WF spread 14–24 Å² | 解纠缠/投影/能窗不佳 | 检查外窗、冻结窗、投影和收敛 |
| 谱图缺少 Γ/A 口袋 | 只做了 Ni 5 轨道 H(k) | 使用 Ni-d + La-d 10 轨道模型 |
| 谱线整体平移 | `MU_TB`、`MU_DMFT` 或 `add_mu_tb` 错 | 统一能量零点 |
| 谱线异常偏移 | DC 重复扣除 | 自能中扣过 DC 后传零 DC |
| 只显示两个自能却循环五次 | e_g 与 5d 配置混用 | 根据 HDF5 block 改 `N_CORR_ORB` 和映射 |
| 图像好看但形式不严格 | PLO/W90 基底未做旋转 | 求 \(R\)，或完整验证基底一致性 |

---

## 12. 哪些结论可以称为“完全正确”

可以确定：

1. 10 轨道 Wannier 模型的构建流程已实际跑通；
2. 最终解纠缠明确收敛；
3. `wannier90_hr.dat` 为 10 轨道；
4. Ni/La WF 中心和 spread 合理；
5. 该 \(H(\mathbf{k})\) 已成功与实频自能一起生成合理 \(A(\mathbf{k},\omega)\)；
6. 化学势、DC 和轨道顺序在成功脚本中有明确处理。

仍需谨慎表述：

1. 论文用 WIEN2k，而这里用 VASP；
2. 论文相关子空间是两轨道 \(e_g\)，成功脚本读取了五个 Ni-\(d\) 自能分量；
3. PLO 与 MLWF 是独立局域基底；
4. 若没有显式 \(R\) 矩阵，直接按轨道名称嵌入是经过物理验证的近似，不是严格恒等变换；
5. 因此这条路线可以称为“VASP/PLO + 独立 W90 的成功工作流”，但不能在未完成基底旋转验证前宣称与论文逐点数值完全相同。

---

## 13. 官方资料

- [VASP Wiki：LWANNIER90](https://vasp.at/wiki/LWANNIER90)
- [VASP Wiki：WANNIER90_WIN](https://vasp.at/wiki/WANNIER90_WIN)
- [VASP Wiki：NUM_WANN](https://vasp.at/wiki/NUM_WANN)
- [VASP Wiki：Constructing Wannier orbitals](https://vasp.at/wiki/Constructing_Wannier_orbitals)
- [Wannier90 User Guide](https://wannier.org/user-guide/)
- [Wannier90 参数源码文档](https://www.wannier.org/ford/sourcefile/parameters.f90.html)
- [TRIQS/DFTTools：VASP interface](https://triqs.github.io/dft_tools/latest/guide/conv_vasp.html)
- [TRIQS/DFTTools：PLOVasp](https://triqs.github.io/dft_tools/latest/guide/plovasp.html)
- [solid_dmft 文档](https://triqs.github.io/solid_dmft/)
- [论文 arXiv:2606.00223](https://arxiv.org/abs/2606.00223)

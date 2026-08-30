# 4. Wannier90 Hamiltonian

0.本目录用于从固定电荷密度的 VASP 计算中构造 Wannier90 哈密顿量，为后续动量分辨谱函数计算提供 (H(\mathbf{k}))。

1.主要文件：

- `INCAR`
- `KPOINTS`
- `sub10.pbs`
2.具体流程
这一步跟前边dmft 并行， 前边的3-dmft 已经得到自能，这一步是为了计算谱函数的并行流程使用wannier进行计算。 故这一步要求vasp 是编译了wannier的功能的。
从1-scf 中复制过来 CHGCAR WAVECAR  POSCAR POTCAR  然后运行vasp 和wannier，之后得到wannier的核心文件。一组wannier90 开头的文件。
这一步计算结束，这一步之后得到了包含路径的H（k）接下来就是组装 把 4- 与3- 计算的结果汇总到一起得到谱函数A(k，w)

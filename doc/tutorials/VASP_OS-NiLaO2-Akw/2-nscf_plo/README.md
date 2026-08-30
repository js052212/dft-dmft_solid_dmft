# 2. NSCF and PLOVASP

0.本目录用于在固定电荷密度下进行 NSCF 计算，并生成 PLOVASP 所需的投影数据。

1.主要文件：

- `POSCAR`
- `INCAR`
- `KPOINTS`
2.具体操作：
将1-scf 中的CHG、CHGCAR、WAVECAR 复制到2-nscf 中，然后修改输入文件，进行nscf计算。 最终是要得到 核心文件 LOCPROJ 

# 3. DFT+DMFT

0.本目录用于进行 one-shot DFT+DMFT 计算，并对 DMFT 自能进行后处理。

1.主要文件：

- `dmft_config.toml`
- `1-sigma_postprocess.py`
- `plo.cfg`
2.具体操作
将2-nscf 计算后的文件夹进行复制命名为3-dmft， 创建上边dmft_comfig.toml和plo.cfg 这两个文件。先进行PLO 投影,执行脚本dmft_convert plo.cfg 生成 .h5文件。 
运行solid 进行计算， 结束之后，在输出的文件夹中有另一个.h5 文件，那是dmft 的输出文件，
运行1-sigma_postprocess.py 这个脚本，对计算结果进行一个初步处理， 然后得到一个平均之后的. h5 文件。 这一步计算结束。

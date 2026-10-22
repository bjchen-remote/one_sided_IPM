# 服务器直接部署

把本目录整个 `+ipm` 和 `main.m` 放到：

```text
/data/user/hd58131/ipm/one_sided_ipm_wall_integral_20261007/
```

这是完整包，可直接复制文件，不用zip。目录名沿用原发布地址，当前内容为source26。可选的tests和数据读取脚本一起复制即可；不需要本机runs目录或Mac MEX。

用服务器已有MATLAB执行：

```bash
matlab -batch "run('/data/user/hd58131/ipm/one_sided_ipm_wall_integral_20261007/main.m')"
```

也可把main.m复制到CHESS作业目录执行；它在同目录找+ipm，找不到时使用上面的固定发布地址。输出默认保存到入口脚本所在目录下的新时间戳目录。可用`IPM_OUTPUT_DIR`指定一个新的绝对目录。

`IPM_PREFLIGHT_ONLY=1`只检查配置；`IPM_SMOKE_ONLY=1`只跑一个新初值真实步；`IPM_PORTABLE_ACCEPTANCE=1`运行本包Gaussian两步、重放和面积触发验收。默认完整运行保持1e-5目标及真实0.5τ持稳条件。

续算通过`IPM_RESTART_CHECKPOINT`指定绝对检查点路径。它保留原初值及参考量，不把旧运行改成Gaussian；不兼容的数值后端会明确拒绝。要采用Gaussian须新开计算。

本机R2026a已执行所列测试，服务器R2022b尚未执行；服务器对接由用户处理。长期目标仍在测试，不把短测通过当作收敛。

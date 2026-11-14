# 运行统一 B 样条求解器

```matlab
c = main('config');
c.run.tauEnd = 5;
c.run.outputDir = 'result/medium_0_5';
[summary,s,ctx] = main(c);
```

默认非均匀中网格、窄重叠阶梯张量核心、极坐标环带、紧支撑 Gaussian 源、CFL 0.4 和 RK4。默认目标为 200τ，单次预算 600 秒。 `c.method.geometry.preset` 取 `coarse`、`medium` 或 `fine`，角向分别为 91、181 或 361 点。主入口构建本机采样内核；处理同一个算例时可保持原输出目录冷恢复。`main(c,@(x,y) ...)` 将自定义密度物化为原生样本及带 SHA 的初值清单。

## 形态触发 remesh

`c.method.remesh` 控制 `enabled`、`level`、`threshold`、`minCheckEveryTau`、`checkEveryTau` 和 `stopTau`。默认 25% 与 80% 主峰连墙轮廓相对上次实际换网格锚点变化达到 25% 后更新。观测间隔在 0.005τ 至 0.5τ 内由形态变化率调整。粗网格默认 `stopTau=50`，传入空数组允许后续换网格。参见 [架构](ARCHITECTURE.md#时间与重网格)。

## 冻结运行与恢复

```bash
python3 tools/run_explicit_experiment.py start --case examples/bspline4.json --output result_bspline4 --target 76 --hours 8
```

监督器以本机采样 MEX 和源文件制作冻结运行时，记录实际网格、运行时段及退出原因。运行中在监督目录放置 `STOP` 可在已完成的步端点暂停。移除 `STOP` 后用 `resume --output <原目录>` 恢复，保留原截止时间。更换方法、初值或实际网格使用新目录。显式 seed 通过 `case.run.seed.file/sha256/tau` 声明，`tau` 必须是源检查点的真实时间。

`latest.json` 指向当前接受 MAT 和 SHA。`initial.mat/json` 保存实际起态。`mesh/` 保存换网格事件。`history.csv` 及 `data/` 保存残差历史和原生数据。 `F=W_\tau`，`Q=q_\tau`。原生数据包含规范尺度和四字物理时钟。

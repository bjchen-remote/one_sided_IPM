# ipm_cr · 单侧 IPM 混合网格求解器

统一五次 B 样条离散用于 Poisson、速度恢复和迎风输运。默认非均匀中网格采用窄重叠的阶梯张量核心与极坐标环带。原生点上的完整源进入 Poisson 弱式，外圆采用齐次 DtN。时间用经典 RK4，默认 CFL 0.4。

```matlab
c = main('config');
c.run.tauEnd = 5;
c.run.outputDir = 'result/medium_0_5';
[summary,s,ctx] = main(c);
```

按 25% 和 80% 主峰连墙轮廓相对最后换网格锚点的变化更新网格。纵向尺度由穿主峰的剖面给出，任一层任一方向变化达 25% 即触发。每次更新保留原规范和完整时钟，重建空间算子并求一次完整 RHS。

[使用](docs/USAGE.md) · [39 文件架构](docs/ARCHITECTURE.md) · [独立验证](docs/VALIDATION.md) · [主干契约](docs/MAINLINE_REQUIREMENTS.md) · [分支与复现](REPOSITORY.md)

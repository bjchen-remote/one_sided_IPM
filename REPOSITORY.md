# 候选主干、参考分支与复现边界

当前公开源码以 `main.m` 为唯一运行入口，默认非均匀中网格；核心主节点与 Poisson 节点逐位共用。接口微分修复、形态触发和高阶 remesh 的资格分别登记，不从源码发布推导完整轨道或全域精度已经通过。唯一交付契约是 [MAINLINE_REQUIREMENTS](docs/MAINLINE_REQUIREMENTS.md)，当前接口见 [README](README.md)/[架构](docs/ARCHITECTURE.md)，实际测试范围见 [VALIDATION](docs/VALIDATION.md)。

本次 GitHub 发布由用户明确授权，源自隔离的含时求解器分支；并发稳态研究、实验状态和本机文件不纳入。NAS 同步检查未成功，公开发布不表示 NAS `main` 已同步或完成集成。

默认 FEM 与对流的全部原生节点、物理坐标和数组 shape 已逐位共用，局部 Gauss 点仅作积分；独立 GLL 只留作显式组件参考。旧 GLL 的失败记录保留，不能混作当前原生节点结果。

## 维护面与数据身份

| 路径 | 职责 |
|---|---|
| `main.m`、`src/+ipm_cr` | 唯一路径包装、公开 API、核心数值模块及同一 `run` 演化入口 |
| `data/initial.json`、`data/reference/gaussian_t0.mat` | 原 Gaussian authority、物理参数和固定原墙 I/P/C |
| `data/mesh.json` | 真实非均匀中网格物理轴/度量；粗网格精确逐轴隔点取得 |
| `tools` | 可选本机原生构建、唯一既有冻结监督器、统一验收登记 |
| `tests / examples` | 验证、显式参考比较与只读展示；旧 case/图形仅代表原配置 |
| `docs / test/ALGORITHM_REPORT.tex` | 当前接口、需求与验收；历史推导和记录保留其来源 |

默认值只由 `config` 解析。FEM、谱、FD 明确记录后端；采样 MEX 可选，原生二进制不入库。`src` 与所需工具进入 MATLAB 路径，不使用 `genpath` 混入另一 worktree、旧包或 frozen runtime。`main` 核对所选包属于本 checkout。

源码身份包括相对路径/内容 SHA、`data/initial.json` 和 `data/mesh.json`；选中初始 MAT 与 manifest 另核 SHA。原 authority MAT 的 SHA256 为 `fa80e7d2b4fedb6b28186dcaee844975569fac41da06d26b35a2519d7521db22`。完整时钟、实际网格和原规范缺失不能合法恢复，几何计划在同身份下冷重建。

## 协作与集成

NAS 是协作集成点；每个候选独立分支/worktree，保留主目录并发稳态工作的索引、未提交内容和 pi 工作树。工作会话及推送前执行 `git sync-check`，失败不等于同步；只在干净可快进条件下拉 NAS，不 reset、自动 stash、force-push 或清理另一任务。

只有满足需求并附实际源码/配置/输入和验收证据的候选，才能作为已交付主干集成。GitHub `origin` 是 push-only 发布目的地，禁止 fetch/pull；公开发布须另行授权，私有研究分支及结果不自动上传。

`run.seed` 的跨方法/网格导入是新运行，记录原文件/SHA/实际 τ 和转移差；不是旧源码的同方法恢复。改变方法或初始输入建立新目录。生成数据、锁、MEX、日志、运行状态与大文件不进入公开源码 Git；科学包、接受态和 remesh 谱系按 [USAGE](docs/USAGE.md) 保存。

## 历史与参考算法

完整旧代码保存在 `legacy/solvers`，历史分离点 `2010a4e2e0714eb092ee0a6207e834b52844ced8`。其 Cartesian、纯 polar、SDIRK2、SSPRK3、旧 remesh 和管理脚本仅按原冻结身份复现；不以当前包重标旧轨道。旧 Q4 FEM、谱和 FD 的证据同样仅适用于其当时配置。

本地 `time_evolution/`、`solver/`、`experiments/`、`test/prototypes/`、历史 inputs/results/runtime 和共享 pole/jet helper 保留原位，不能因“统一主干”而删除。独立稳态不动点以原研究文档为入口，不替代本含时主干契约；旧 native steady-Q 保留为历史，不重新成为默认。

旧均匀核心的 `0→10τ / 490→500τ / 粗500`、旧 FD 小粗网格的 `0→400τ` 都不能充当本轮缩短的非均匀中网格 `0→5τ / 390→395τ` 与精确半密度粗网格 `0→200τ` 的验收。过去文件的详述可从原提交、[清理记录](docs/CLEANUP.md)和 [VALIDATION](docs/VALIDATION.md)追溯；缩短当前接口文档不改变历史资格。

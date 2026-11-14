# 维护与验证

- 生产主链是 `main → run → step → rhs`。网格更新由 `remesh` 组织，空间算子由 `context` 绑定到实际坐标。文件职责见 [架构](ARCHITECTURE.md)。
- 修改空间离散时复核原生坐标、域外模型、单侧墙源、人工迹和活动 donor。空间、时间和转移误差由独立测试量测，正式演化记录被动诊断。
- 修改存储时复核 owner、读回、SHA、原子指针、移动结果目录及同版本冷恢复。以 [主干契约](MAINLINE_REQUIREMENTS.md) 为部署验收依据。
- `tools/run_validation.py --suite full` 调用现有 API 的独立测试。报告目录保存每项状态、源清单和日志。早期参考算法及其专用检查随对应 Git 版本保存。
- 提交前运行 `git sync-check`。改动在候选工作树完成并独立验证后评审，保留并发工作和冻结运行时。

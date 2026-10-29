# 求解器维护入口

维护对象是 `src/+ipm_cr` 的含时求解器。数值输入、历史结果、冻结运行时和独立稳态研究不随代码整理重写。

## 代码调用链

```text
main -> config / problem -> init -> context
run -> step -> time.step -> time.trial -> rhs
rhs -> mesh.fill / mesh.diff -> poisson.apply -> psolver
    -> gauge.wall / gauge.rates -> op.prepare -> 完整 F / Q
run -> io.checkpoint / diag.overlap -> pack
```

| 位置 | 维护职责 |
| --- | --- |
| `+mesh` | 实际节点、坐标度量、活动区域、ghost 交换和局部重构 |
| `+poisson`、`psolver` | 谱源映射、分片径向系统、core Sylvester 解和 overlap Schur |
| `+gauge` | 固定原墙上的 I/P/C、规范闭合和真实 RHS 审计 |
| `+model`、`rhs` | 解析背景、残差变量、完整物理量和标量速率 |
| `+time` | RK 阶段、完整 CFL、峰事件、试步事务和补偿时钟 |
| `+io`、`+diag` | 来源身份、接受检查点、独立观测和数据交付 |

公开接口及数据含义见 [ARCHITECTURE](ARCHITECTURE.md)，实际运行见 [USAGE](USAGE.md)。默认测试由 `tools/validation.json` 唯一定义；默认检查不依赖本机私有历史结果。

## 维护顺序

先从已发布版本建立独立维护分支，执行默认完整测试，再提交整理基线。数值修复使用后续提交，每项包含能区分旧行为的回归；历史运行保留原源码身份。

2026-10-29 审查确定后续维护范围：端点补偿时钟与 RHS 一致性、规范根条件数、背景相消与 Q/F 离散相容性、谱源与绝对误差审计、外边界和固定墙规范的误差界限，以及热路径重复计算。数学近似的剩余误差必须明确报告，不能以小线性残差或测试完成替代精度认证。

效率修改需要先测再保留。固定差分计划缓存曾在旧版本微基准中变慢，因此不能仅凭减少分配次数认定更快。独立端点审计、完整阶段 CFL、Poisson 回退、写者锁和双槽检查点继续保留。

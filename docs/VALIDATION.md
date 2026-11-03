# 当前公开默认与接口物理微分（2026-11-03）

本版默认非均匀中网格（core 129×241、polar 549×193）、原生节点 P5 FEM、linear5 五阶迎风、平滑速度系数 .1、RK4 和累积形态触发 remesh。central4 保留为显式参考。main/config/case/公开示例解析同一默认；原生内核的构建依赖与冻结管理器已同步。

人工接口附近 q 和 Q 使用同 patch 活动节点上的七点物理导数；内部中心模板、物理墙法向、外圆及五阶迎风 ghost 模板不变。不投影活动重叠，未加生产精度门。q/Poisson 源和 Q 同时改变，geometry key 标识新闭合；旧轨道不作为新闭合的同方法恢复。

在原有实际 quarter/coarse/medium 轴及相同两个解析 scalar fixture 上，coarse→medium 的接口 Bx 观测阶从3.81561/3.25808变为4.79969/4.08208，medium 接口绝对误差从6.32951e−6/2.94596e−5降为2.06350e−6/1.07520e−5。全活动 Bx 阶4.97477/4.44768。B 模板未变；对照和误差分解表明原来的主要项是已离散 B 的截断误差经跨片同步后再求导，并非连续解析 B 的值插值误差。条带随网格移动，这些有限对观测不构成通用渐近阶证明，也不是完整 IPM Q 的连续误差。

新默认的12项相关注册检查分别在 validation_03（前6项）及 validation_04（其余6项）满足原门，生产数值源清单完全相同；入口测试改为检查解析后的默认方法，不直接读取尚未解析的 case 可选字段。validation_04 再做原完整 RK4 三步长时间试验：W=3.962435阶、q=3.955132阶、Q=3.810904阶；F低于原分辨屏，不赋阶。时间屏通过，所有原阈值保留。该次只请求 temporal，registry 的整体 mainline_accuracy 仍明确为测量完成但空间资格未请求/未认证。另7项Python冻结管理器检查通过。

同一个 τ143.0074 实际粗网格态的固定几何0.1τ短诊断，两种闭合各13步完成；新闭合 Q∞从8.12720e−5降至7.00191e−5，最大阶段 CFL=.285（配置.3），原始参考及全部时钟保留。旧闭合相同诊断从5.42975e−5到5.41187e−5；同态不同离散残差不作真误差比较。此短诊断不证明长期稳定。

本地原始证据位于 results/interface_accuracy_20261103_01；validation_03 execution SHA256 `d9636ed7b8b1631e35e7e924a8a909f96b2d1950eda5a3180881d3fc3f2fc107`，validation_04 execution SHA256 `6fa679d32c37430bcd37273b495098322bc1ff2b713a331e4844b79ceba6bf9a`。完整新闭合轨道、统一空间/远场精化、速度导数连续精度和整套主干交付仍未认证。原冻结bc61655的粗τ200轨道与此新闭合分别保存，不能混写来源。

---

# 当前形态触发修复（2026-11-03）

每2τ强制重分布是实现自行加入的错误规则，不是用户要求。旧报告使用25%正源主峰连墙区域，与最后实际 remesh 后尺度比较，20%累积变化才触发。新分支恢复此规则；首次只校准，检查不重置参考，未到形态阈值不转移、不重建RHS。以下已有组件与长轨道属于其冻结祖先，不构成本次修复资格。新检查完成后单独登记。

# 维护验收记录

## 形态触发修复的实际回归

- 首轮7项相关检查通过，62.6613秒：region几何、累计触发与同τ转移、三网格转移阶数、run/seed/冷恢复、main、config、arch；源码与测试输入身份未变。
- 随后恢复旧的自适应观测节奏（.005至.5τ），仅改变控制器检查调度/记忆；相关remesh、mainline_runtime、config三项复测通过，50.3975秒，源码与测试输入身份未变。
- 固定形状运行至任意τ无插值；10%/15%累计不足不触发，25%累计触发；真实context的非force阈值事件、一次fresh RHS、参考仅在实际换网格后更新、时钟及冷恢复检查通过。
- 证据分别为 `results/remesh_trigger_20261103_01/validation_01/execution.json` 与 `validation_02/execution.json`。该组件结果不证明长期残差复现；错误周期迎风实验已在τ79.714停止，数据保留。

## 当前 remesh 修复候选（2026-11-03）

候选尚未集成或部署到 NAS `main`。本次修复已有同 τ 转移、完整 RHS 三层阶数和实际单事件复算证据；32项 `repair_release_01` 统一检查已通过（282.715秒），源码与测试输入身份均未改变。**修复版长期独立残差是否恢复至旧实验的 `1e-9` 仍未验证，三段修复版复现与整套交付尚不能宣称通过。**

本轮统一检查证据为 `results/mainline_adaptive_20261102_01/repair_release_01/`，execution SHA256 `fd5ea36a15003f247be0b89a650b9ddd2c6e02db6519efa8adf3af36ef3080ac`。32项请求检查通过；默认空间与完整RK4时间屏通过，W时间观测阶仍为3.9606。旧radial与非默认P4的资格限制保持原样，未改门限。

修复前 `15efb09` 的周期 remesh 把略越出旧活动边界的新节点切换到另一独立 patch 的 donor。独立常数 `core=.125/polar=.375` 在新暴露条带上直接换值，旧测试只检查旧活动覆盖内的点，漏掉该路径。修复后，节点数及完整轴端点相同的再分布在各片自己的活动轴上作最近七点张量重构，边界条带继续本片端部多项式；不读 ghost、不覆盖活动重叠、不混投影。通用映射默认仍禁止外推，改变 patch 框或节点数的导入仍按覆盖路由处理。没有增加运行时误差门、拒绝或回滚。

| 独立检查 | 已执行结果与范围 |
|---|---|
| 真实 coarse `τ=155.329` η 的同 τ 往返 | core 最大值误差 `6.6531e-3→2.8246e-5`，D7 导数差 `1.2482e-2→2.9947e-5`；这是转移一致性诊断，不是解析真误差 |
| 完整同 τ remesh，L=12/24/48、非零解析 η、同新网格/context 的直接解析态参考 | 4.285秒；W/q/F/Q 的原生活动点和可分辨公共查询误差对通过近四阶屏；见下表 |
| `test_remesh` | 5.75秒，通过全活动独立常数、六次多项式、平滑非多项式、原参考/四字时钟、周期及冷恢复检查 |
| 实际 coarse `99.408875847→100.35` 单事件复算 | 179步；在同一 `τ=100.148927761`，旧 frozen 版本完整活动点 `Q∞: .0029180875→.036561329`，修复版 `.0029180875→.0029185270`，增量 `4.395e-7`；此为纯同 τ 前后比较 |

三层解析 η 为 `.02*sin(.8*X+.5*Y).*exp(-.3*(X.^2+Y.^2))`，不要求旧离散态残差为零；各层先 transfer 到新网格，再与同网格直接解析构造态的完整 FEM/输运/I/P/C RHS 比较，无时间推进。原生活动点完整最大误差如下，阶数只表示这些有限层、该平滑源与该小幅移格试验：

| 场 | L12 / L24 / L48 绝对最大误差 | 相邻观测阶 |
|---|---|---|
| W | `4.2417e-10 / 7.9943e-12 / 5.4345e-14` | `5.730 / 7.201` |
| q | `2.2581e-9 / 7.4113e-11 / 1.1144e-12` | `4.929 / 6.055` |
| F | `1.1793e-8 / 2.3160e-11 / 4.1661e-13` | `8.992 / 5.797` |
| Q | `1.1657e-7 / 3.7379e-10 / 9.0501e-12` | `8.285 / 5.368` |

固定公共查询只由活动误差重构；common W 的细层进入声明舍入屏，不赋阶。普通 `mesh.sample` 在活动边缘可读同步 ghost，其公共读数差不能直接叫纯 transfer 误差。实际网格非多项式对照中，wall-layer 值误差由 `2.22e-4` 降至 `1.06e-6`；Gaussian/振荡例的小截断项有所增加，窄层仍未分辨、值误差约 `.10`，均保留原报告，不据七点次数宣称普遍高阶、无过冲或每项误差都更低。

证据：`results/mainline_adaptive_20261102_01/repair_remesh_01/{interpolation_accuracy.json,same_patch_interpolation_accuracy.json}`、`repair_order_01/remesh_order.json`、`repair_remesh_unit_01/remesh.json`、`repair_event_frozen_baseline_01/` 和 `repair_event_replay_01/`。修复前 HW3 在 `τ=155.329` 由 owner STOP 保存，未达到200；旧结果及其扰动记录不删除，也不改标为修复版成功。

## 修复前直接求解版本（2026-11-02）

`compute_only_release_05` 的31项请求检查通过（281.328秒），当时源码与测试身份保持不变；默认空间屏和完整 RK4 时间屏均通过，W 时间观测阶3.9606。该冻结版本随后暴露上述周期 remesh 条带问题，组件通过不能代替轨道复现。

证据：`results/mainline_adaptive_20261102_01/compute_only_release_05/`；execution SHA256 `792ebfa1fbfbf9563f155254582420f9302b0d9318888e22edeabaf68e2cbeee`。旧 radial 参考完整资格仍为 false、非默认 P4 两层筛查仍失败，不能冒用默认 P5 的资格。F 时间差细层低于 `1e-12` 屏，不造阶；Q 为3.8398阶、q为3.9537阶。

生产主链已移除数值精度门、诊断拒绝、备用算法、重试/回滚、最小步长和在线 steady 停止；CFL 只参与当前速率的步长公式，阶段速率被动记录。Poisson 每次 patch 线性解固定作一次原弱残差修正。周期 remesh 使用六项配置，首个完成端点即适配，此后按 τ 间隔，根据当前解与固定参考核心有界重分布；活动七点转移、重建 context、一次 fresh RHS 后直接继续。

`compute_only_components_02` 的 11 项检查已通过，76.071 秒。随后修正 arch 的重复 Newton 末位差测试语义（实测 F 差 `4.996e-16`、W/q/速度与隔离 context 逐位一致）、入口测试混用初始/首步 remesh 几何，以及参考 FD 测试期待旧缺状态保护的错误用例；独立物理精度门未放宽。

`compute_only_release_04` 的31项请求检查已执行完成，266.390秒，源码/测试身份未变；**时间资格仍为 false**。默认 P5 实际网格空间屏通过。未分段 RK4 的原 `τ=.4`、`.01/.005/.0025` 对照中 W 差为 `7.482e-12 / 1.840e-12`，有限对阶2.024；logC_l阶1.801。旧/新端点差在 `τ≤.11` 约 `1e-15`，主要增量出现在原样条节点跨越 `.115346` 和 `.337327`。不能把它改标为未分辨、通用二阶或完整通过。

当前仅恢复原 IPC 分片光滑问题所需的数学分段：普通连续单峰跨相邻 knot 时固定四次弦截构造 RK4 段；不恢复旧残差门、阶段折半、括根备用法、重试、拒绝或 minDt。节点 ULP 仅作坐标身份。新原全窗缩步已由独立 temporal-only 及 suite05 复测通过；不借光滑小窗或旧3.9566阶继承资格。
下面 `b890305` 的统一、空间/时间及独立 Green 证据保留原来源。新版本需要独立完整验证；正式实验不设置精度门。用户已恢复并减半窗口：中网格 `0→5τ`、真实原态 `390→395τ`、精确半密度粗网格 `0→200τ`，主要与已保存的同 τ 旧数据独立比较，不拟合平移、振幅或时间。已完成 Green 积分不重复计算。

## b890305 原生网格版本的已执行验收（2026-11-02）

`b890305` 冻结版本的统一 `mainline_release` 已完成：31 项请求检查满足，单线程 MATLAB 含启动耗时 237.630 秒，生产源码、输入与测试 SHA 前后不变。旧谱径向测量的完整空间资格仍为 false；参考 P4 的实际网格有限对照也未通过其空间阶屏。它们与默认 P5 资格分别记录。

该冻结版本采用非均匀原生节点 P5/P6 Galerkin、Gauss8、缓存均衡稀疏 Cholesky、central4＋Q6 耗散、RK4 和旧多门自适应核心 remesh。主体节点、边界节点、物理坐标及 raw 系数均与输运网格逐位一致；用户允许的低维边界节点例外本实现未使用。源无关的几何条件优化只重排 P6 余单元；最外径向整单元为 P6，所有速度仍由同一 FE 流函数求导。

| 测量 | b890305 冻结源的结果与范围 |
|---|---|
| 默认 P5，实际精确粗／中网格，固定查询云速度最大误差 | 原生源 `7.13824e-5 → 9.15766e-7`，有限两层阶 6.284；高于预登记屏，但不是通用六阶证明 |
| 非零外墙角源，64／128，解析／原生源及 Dirichlet／DtN 四组合 | 完整速度最大误差约 `8.45e-8 → 4.30e-9`；完整阶 4.293–4.296、角点阶 4.024–4.130，原门限通过 |
| 完整 RK4 固定粗网格，τ=.4，最大步长 .01/.005/.0025 | W 3.95663、q 3.95347、Q 3.93621；F 与部分标量细层低于舍入屏，保留未分辨标记 |
| R=1e10 常数＋解析源等价拆分 | 场／系数差 `1.42e-11`，常数采样精确；原门限通过 |
| 暖原生源 Poisson，粗／中 | 中位 .01384/.07821 秒；节点约3.968倍、时间约5.651倍，经验指数1.257；稀疏因子存储指数1.129 |
| 旧换网格与存储 | 该版本的真实触发接受、拒绝回滚、控制器重试、完整时钟、冷恢复、种子换位置及事件 SHA 检查通过；不认证新简化 remesh |

上述成本不含构建、Green 全局缓存、完整 RHS、输出和 remesh；没有 `O(N)` 证明。C0 ψ 不保证速度 C0/C1；实际单元面速度跳跃和 native D7 速度导数均有独立原始测量，不从 ψ 次数推导导数阶。有限盘零外源模型不认证未知物理尾部。

证据根为本地 `results/mainline_adaptive_20261102_01/native_release_01/`，`execution.json` SHA256 为 `88708dde2f3286ffed77a0bbce668696a7114787c8f42f523c6d87fad695b49e`，`components.json` 为 `3ad06a277c36eae327855d643de1382686c1276108d99f8440fdd53eb0a33dbb`。生成数据不进入源码 Git。

旧多门 remesh 的轨道在接受态边界正常停止并保留：

| 旧轨道 | 实际停止状态 | 解释 |
|---|---|---|
| T1 中网格 0→10 | 正常 STOP，接受至约 `τ=1.2878` | 未达到 10；不是数值完成或简化版通过 |
| T3 精确半密度粗网格 0→400 | 正常 STOP，接受至约 `τ=23.297`、3680 步 | 保留旧往返转移误差门拒绝 9 次；未达到 400，不移标为新版本通过 |

缩短后的 T1 中网格 0–5、T2 真实 390–395、T3 精确半密度粗网格 0–200 仍须在同一新源码下完成管理器试跑、恢复与执行。真实 τ390 输入已从 NAS 获取并逐位核对原中网格，SHA 为 `f2688332431fb7193830a73cc03c0805ba0a9fddbbdf44ad1a122c45f0229f78`。目前不能将组件通过写成主干全部交付。

失败证据保持不变：旧 GLL/P5 64／128 完整速度仅约2.6阶；第一版 native P5 约3.116阶。70位精确节点插值复现了最外单元导数的粗层 Taylor 抵消，支持整外单元 P6 策略；它不是独立 Poisson 求解资格。第一版拉伸转折单元的病态放大源拆分舍入，修复后按原门复测通过。详见 `gll_stage_snapshot_01/`、`p5_corner_diagnostic/run_01/`、`native_wall_interpolation_diagnostic_01/`、`native_split_diagnostic_01/` 和前序 `native_components_* / native_precision_01/`。

## 历史 GLL/Q4 FEM qualification (2026-11-02)

The selectable Q4/Gauss8 FEM backend preserves the original weak-row, trace and essential-boundary gates. A shared analytic polynomial Green contribution removes the finite-source-disk outer-wall corner terms; radial elements are grouped inward from the outer circle.

All **24 registered full-suite checks passed** against the release source, including the original solver checks and new `fem`, `fem_green`, and `fem_wall` tests. Independent image-Green quadrature fixtures validate all 32 supported polynomial monomials; normalized streamfunction/velocity differences are below `2e-16` in those fixtures.

For `q=1+0.1*(1-r^2/R^2)^4` inside `R=8`, zero outside, both outer-wall source values equal 1. Full-array velocity errors at levels 64/128/256 are `9.4954e-7`, `5.7377e-8`, `4.0674e-9` (orders **4.049, 3.818**). The last wall-only error is `1.4136e-9` (order 3.826). Native and analytic source paths and Green Dirichlet/DtN boundaries agree on the 64/128 order measurement. Off-mesh and corner-neighborhood errors, the failed unenriched controls, coarse pre-asymptotic orders, and the core roundoff floor remain recorded. This is a tested smooth-source qualification, not a universal continuum bound.

Evidence: owner-local `results/fem_release_20261102_01/regression_01/` and `results/fem_wall_order_20261102_01/`. The three requested trajectories use the original upwind transport and unchanged RK4/normalization. All bounded pilots and cold recovery checks passed.

| Requested interval | Actual result | Accepted steps | Retries | Maximum stage CFL |
|---|---|---:|---:|---:|
| Medium 0 to 10 | Reached 10 | 1542 | 0 | 0.286661 |
| Medium 490 to 500 | Reached 500 from the actual original 490 checkpoint | 1140 | 0 | 0.285001 |
| Coarse 0 to 500 | Failed at 138.10811799 | 9258 | 0 | 0.288467 |

The coarse run hit the original `PeakEventUnresolved` guard. A separately identified warm half-CFL probe, with every peak/gauge/CFL guard retained, accepted one more step and failed again at 138.11055034. Two distant local maxima have closing height gaps; the observed behavior supports an impending global-maximum branch exchange, not a continuously moving peak crossing the requested interpolation knot. Neither saved state is an accepted exact tie. This is not evidence that the continuous PDE is unstable. The control is not spliced into a claimed 0-to-500 success. The requested coarse target remains unmet, so deployment to NAS main and GitHub is held.

The completed medium intervals are finite-time integration results, not steady-state or continuum-accuracy certificates. The late endpoint has near-field RF `7.82e-5` and RQ `1.329e-3`; bounded completion does not imply a residual convergence gate passed. See [FEM algorithm and reproduction commands](FEM_POISSON.md).

## 独立同节点 FD 后端（2026-11-01）

`same_mesh_fd` 是显式选用的实验后端；默认仍为 `spectral`。本轮在独立分支运行最终 `full` 的 21 项请求检查，全部通过，含 MATLAB 启动共 80.053 秒。原 radial 项仍保留 `qualificationSatisfied=false`，未将舍入受限测量改成完整空间资格。新增 `test_same_mesh_poisson` 覆盖均匀/sinh core、两种非零墙源制造解、活动点/ghost/重叠误差、解析源分离、同一流函数速度、后端身份及四词时钟。

公开 FD API 与实验版本在两个实际初始状态上的 W、q、完整 F/Q、两分量速度及规范速率逐项相同，最大差为零。另以原始 Gaussian 初值在粗 sinh 网格完成 τ=.003 的三步运行；第一步后保存并冷恢复，与不中断运行的尺度、W、q、完整 F/Q 逐项相同，导出后端身份及数据包 SHA 检查通过。这些是接口、短程恢复与代码提取等价性证据。

最终自足测试还覆盖 FD+sinh+非零时钟尾字的真实 checkpoint/pack 保存读回、context 重建和两个微小步对照，以及完整尾字的采样和 summary 导出。这是合成持久化回归，不是另一条物理轨道实验。混合精确节点和非节点的两点采样暴露了原张量插值器索引块行/列方向不一致的问题；修复只统一稀疏三元组形状，保留权重和组装次序，最终全套检查包含该实际触发路径。

制造解和同精度成本的范围见 [同节点 Poisson 后端](SAME_MESH_POISSON.md)。有限离散方程残差通过不等于连续精度、全局四阶或长期稳态资格；实验后端未继承谱后端的连续残差验收含义。

固定粗网格 τ400 已完成，六项区域 F/Q 残差和最后三个 W/q 变化窗口通过，新建因子完整 RHS 与保存末态最大差 1.42e-13；原生 F 变化窗口仍未通过，故预先要求的双判据未通过。此结论与组件回归全部通过分别记录，未以回归通过替代稳态确认，也未延长运行或执行推送。

## 六项完整性维护（2026-10-29）

先核对既有 GitHub 维护对话和已发布 `00403aa`，在独立工作树整理调用链；16 项默认维护检查通过后发布文档基线 `4a05e86`。本轮不纳入并行稳态研究源码，不重写既有检查点或冻结运行。

修复覆盖补偿阶段/端点时钟、规范根条件数及非唯一峰、活动点完整 F 的背景时间项相消、Q 表示缺陷与原生源 Poisson 离网格诊断、缺失外圆入流拒绝，以及峰扫描/固定谱插值/稀疏径向装配的重复成本。原 I/P/C 定义、固定参考、完整阶段 CFL 和原 Poisson 门不放宽。新 native defect 门和绝对/去平移速度门仅可附加；同一源重构的缺陷不冒充独立连续误差。

最终 `full` 的 20 项维护检查全部满足请求条件；再显式加入只跑时间的三层短缩步，共 21 项请求测量完成，含 MATLAB 启动总计 88.109 秒。仅时间测试耗时 39.044 秒，固定 h=.06、nr=193、ntheta=64、谱阶 48/64/16，τ=1，名义步数 64/128/256、实际 68/132/260；W 观测阶 3.943378，耦合标量阶 3.941980，全部五个标量阶 3.927–3.951，均通过原舍入屏及时间阶窗口；最大阶段 CFL=.247595，配置为 .3。未运行空间大实验，radial 的完整空间资格和 order 的全局资格仍为 false。

新增回归确认：2^53 时钟高部仍保留小增量；RK4/SSPRK3 各一短步保存的端点与重新求值 W/q 逐位相同，F/Q 差为零；Q 与独立四阶中心差分的 q 时间导数差约 8.0e-13。测试粗网格的 `Q-D_X F` 缺陷仍可达到 core 2.41e-7、polar 9.51e-4，已明确报告而非隐藏。固定墙 I/C 的连续求积偏差仍呈二阶；其定义与参考没有静默更换。

性能以同一 MEX、两个独立单线程进程、固定 t0 状态、20 次暖机和 5×10 次完整 RHS 测量。原版两轮内部阶段中位成本 25.41–26.23 ms，最终版 26.45 ms；完整诊断原版 19.31–19.33 ms，最终版 20.61 ms。张量收缩回收了新增审计的部分开销，但完整端点诊断仍约贵 6.6–6.7%，未宣称净加速或轨迹加速。差分缓存候选只做 ABBA 微基准，未采用为生产优化。

本地证据在 `results/maintenance_final_20261029_01/`，性能记录在 `results/maintenance_performance_20261029/`，结果与二进制不入公开 Git。最终 execution SHA256 为 `e5bb46a181b5d16e44731d68d1250936cab74c805ab85de7faeca7da34782062`，components SHA256 为 `e0dd52225e5e8f74fc4512a76ea4ddce21d8a227c196469b4ed15ae4373b7f7a`。这些是本轮维护、短时间阶和成本证据，不认证最终残差机器精度、长期稳态或全域空间精度。

## 根目录 MATLAB 入口与 Windows 路径（2026-10-25）

本轮入口为编辑器 Run、`[c,rho0]=main('config')`、`main(c)` 和 `main(c,rho0)`。顶部直接手写 `rho0=@(x,y)(1+tanh(x)).*exp(-y.^2)`，要求向量化、同形状、有限实数。`set_initial` 物化真实 core/polar 和原固定参考轴数值，保存独立 SHA 输入并重算 I/P/C；不保存句柄，不做 Gaussian 拟合。表达式路径为 background none、无 eta，换输运网格须用原表达式重建。默认目标τ=200、预算600秒、项目根result，中层网格及谱次数64/80/20保持入口设置。

新增入口验收需单独保存以下范围及来源，不能继承下面的既有17项结果：

| 范围 | 应核对的行为 | 当前证据状态 |
|---|---|---|
| 只读配置与定位 | 从项目目录、其他 cwd 和含空格项目路径调用；配置正确，`pwd` 不变，不构建/求解 | 待本轮入口实测记录 |
| 本平台内核 | MEX 已有与缺失分支；从源码构建，未选 C++ 编译器给出指引，不自动切换算法 | 待本轮入口实测记录 |
| 函数表达式初值 | 实际输运值、固定参考轴、SHA、I/P/C；background none、无 eta；不拟合、不保存句柄 | Mac 短程与恢复已完成，见下文 |
| 短程与恢复 | 小预算接受态、同目录恢复、含空格输出和目录迁移；错误初值/方法身份拒绝 | 待本轮入口实测记录 |
| Windows 文件与路径 | 发布文件无保留名、大小写冲突及非法组件；盘符/UNC 路径识别 | 静态审查与实际 Windows 执行分开报告 |

Mac 已用表达式初值推进至 τ=.002、3个接受步，并成功恢复，记录为 `result/windows_main_20261025/quick_initial.json`。这是入口与短程恢复证据，不是精度资格；其他入口检查的证据位置和 SHA 待补充。Windows 实机的 MEX 编译、短程运行与恢复尚未验证，macOS 路径或文件名静态检查不能替代。纯 MATLAB `main` 不要求 Python；Python 监督器和验证启动器仍限 macOS/Linux。此次修改不新增空间阶、长期轨迹或稳态资格。

## 原始参考与网格收拢（2026-10-24）

本轮仅修改参考初始化链与无消费者的基础网格派生量。源码从 66 个文件、5709 行减至 65 个文件、5496 行。

- 新 `reference(cfg)` 的原墙、尺度、I/P/C 参考和初猜与旧三载体输出逐字节一致；完整初始化数值态及全部非身份/计时 RHS 诊断一致。
- 两种输运网格的完整几何（仅排除 buildSeconds）、基础网格的全部保留字段及网格导出数据逐字节一致。
- 两种网格各 τ=.006 的完整短轨道及公共查询逐字节一致；上下文、原始 authority 正反例、冷恢复和目录迁移通过。
- 径向重分配结构检查通过；与冻结原版复算的 18,413 个数值及 174 个逻辑值完全一致。所有可分辨观测阶均≥4；唯一未分辨项细层误差 2.56e−13，低于声明分辨阈值 1e−11，保留未分辨资格标记。

单线程、每次新建上下文的完整初始化，3 次暖机后测 7 次，中位耗时由 0.2861 s 变为 0.2358 s。原版样本存在时间漂移，此记录只反映本次初始化测量，不推导长轨道加速比。

本地原始记录位于 `results/solver_condense_20261024/`：`before/` 为修改前快照，`authority_before/` 为参考/初态/几何基线及 checker，`final/` 为最终验证，`redistribution_before/comparison_with_candidate.json` 为静态逐数字复算对照。原历史与冻结输入均保留。该维护验证没有新增 τ200 或全局空间四阶资格。

## 求解器收拢（2026-10-23）

初始化与背景分发直接收进公开 `init/rhs`，推进、诊断和输出全程使用 `W/eta/rhs`。移除 5 个源文件及不可达的旧初值分支，`.m/.cpp` 源文件由 **71 减至 66**、总行数由 **5863 减至 5709**。Gaussian 计算式、运算顺序、全部阶段 CFL 和原 Poisson 验收保持不变。

| 检查 | 最终源码结果 |
|---|---|
| 原网格 / 远区重分配，各 τ=.006、6 步 | 初末完整数值快照、1587 个公共查询及诊断与修改前冻结版逐字节一致 |
| 无背景 / 自定义背景宽度 / SSPRK3，各 τ=.003、3 步 | 初态、第一次暖 RHS、完整重建密度、末态及各步全部非计时信息逐字节一致 |
| 状态与上下文 | 统一字段；上下文隔离、错误身份拒绝、同身份 refresh、独立时间控制通过 |
| v2 运行事务 | 冷恢复与目录迁移差为零；STOP、SHA 损坏、过期游标拒绝及 v1 原始读取/只读规范化通过 |
| 椭圆与背景组件 | 独立制造解、原门限及无效源拒绝、两种外边界的线性叠加、解析背景、主尾及原规范通过 |
| 完整 RK4 时间缩步 | 固定空间 τ=.4，20/40/80 个名义时间步；W 观测阶 **3.9715**，耦合尺度约 **3.96–3.97**；X_shift 细层进入舍入区，未声明全部标量达阶 |

时间检查仍通过真实峰事件拆步，最大阶段 CFL 为 .244323（配置 .3）。无背景初始化后第一次暖求根的 F 差约 4.42e−17，在修改前源码中同样存在；新旧两次调用各自逐字节相同。验收比较同样调用序列，不以重复求解必须恒等替代原求根契约。

单线程 MATLAB R2026a、相同 MEX，在 h=.03、基准 nr=385、ntheta=129、谱阶 64/80/20 的同一 τ=.003 接受态，50 次暖机后测 40 次完整诊断 RHS：

| 成本 | 修改前 | 精简后 |
|---|---:|---:|
| 暖 RHS 中位数 | 40.078 ms | 37.662 ms |

计时前后完整状态与返回数值逐字节一致。仅测固定状态的 RHS，排除启动、建图和接受步成本；此测量未见性能回退，不能外推为整段轨道加速比。最终相关检查在一个独立 MATLAB 进程中完成，内部合计约 34.77 秒。

本地证据位于 `results/solver_condense_20261023/`：`before/` 是逐文件 SHA 核验的修改前源码和基线，`branches_before/` 是额外分支基线及计时脚本，`final/validation.json` 是最终验收，`code_audit.json` 是源码清理清单。历史结果及冻结 runtime 未改写。本轮没有重做 τ10 空间研究，也没有新增 τ200 精度或 50 分钟预算资格。

## 架构迁移验收（2026-10-22）

本轮验证针对 `src/+ipm_cr` 的新目录与接口，保留数值算术、完整阶段 CFL 和原 Poisson 验收门限。历史轨道与冻结运行未改写。

| 检查 | 结果 |
|---|---|
| 原网格/远区重分配，各 τ=0.03、2 步 | 初末完整状态、导数、速度、15975 点查询与原冻结基线逐字节相同 |
| 公开接口 | W/eta/rhs 无重复别名；mesh 为结构体；通用数组查询无需区域掩码 |
| 上下文 | 独立 solver/计划；错误 ctx、变更网格/背景/重构后复用旧 RHS 被拒绝；时间控制可独立改变 |
| v2 检查点 | τ=0.003→0.006 冷恢复和目录迁移差异为零；SHA 损坏、过期写游标被拒绝；STOP 保持指针；v1 只读导入 |
| 干净导出 | 无旧 solver 目录，从头编译；上述回归、接口、恢复和 native 检查全部通过 |
| 原生内核 | 五次多项式、常数消去、大常数平移、速度偏导和完整 CFL 检查通过 |
| 独立组件 | 谱 Poisson/独立方向加密、线性叠加、连接算子、解析背景、主尾、原规范检查通过 |
| 完整 RK4 | 固定空间 τ=0.4，20/40/80 个名义时间步，实际 22/42/82 步；W 时间观测阶 3.972，主要尺度/时钟约 3.96–3.97；X_shift 细层差分进入舍入区，未声明全部标量达阶；保留真实峰事件与全部阶段 CFL |
| 管理工具 | Python 编译、输入路径/复制竞态拒绝、prepare-only、完整组件套件及冻结管理器 τ=0.006 实际闭环执行通过 |

本轮在两个独立 MATLAB R2026a 进程中单线程对比旧快照和新包。网格 h=.03、基准 nr=385、ntheta=129；谱阶 64/80/20、requireAccuracy=true；固定径向重分配 factor=2、preserveRadius=4、fullRefinementRadius=20；RK4/CFL=.3/linear5/解析尾项保持相同。50 次完整 RHS 暖机后计时 40 次：

| 成本 | 旧 | 新 |
|---|---:|---:|
| 完整诊断 RHS 中位数 | 55.844 ms | 55.965 ms |
| τ0→0.1 的 8 个完整接受步总成本 | 1.907 s | 1.926 s |

步数、各步 dt、最大阶段 CFL=0.2862270358 及完整端点数值场一致。排除进程启动、初次构建、暖机、散列和存盘；短程成本约差 1%，未见明显性能回退。首轮暖机不足的结果另行保留，不用于稳态成本结论。

原始本地记录位于 `test/results/architecture_20261022/`（结果与二进制不入库）：`clean_acceptance/` 为干净构建/逐位回归/接口/恢复/native，`clean_suite/` 为 Python 全套实际执行，`time/` 为完整 RK4 时间缩步，`perf_steady/` 为充分暖机成本、脚本和完整端点。`manager_final/` 为冻结管理器的实际短程闭环；`before/` 是迁移前源码快照。

这些是架构等价性、事务、组件及短程成本证据。本轮没有重做 τ10 空间分辨率研究，也没有建立新架构 τ200、50 分钟预算或等误差加速资格。


## 发布入口整理（2026-10-24）

已验证的前批精简提交为 `6d63056`，其 source/data 与既有验收 SHA 全匹配，独立复核17项 proof和315处文件身份引用。该提交只认证维护等价；未把全局空间阶或长轨迹精度改为通过。

本次发布入口新增配置、Gaussian authority、接受态停止窗口、原生数据包及单一测试清单。统一 full 清单加两项显式历史回归共17检查实际通过，MATLAB总耗时55.48秒；源、测试和输入SHA前后相同。包含 native 停止/恢复/导出、候选恢复不推进、导出拒绝初值身份漂移，以及标准/重分配网格独立短轨迹逐位回归。CLI另有4种有效准备、33项无效输入拒绝及自定义authority冻结检查。当前12小时稳态实验仍保持其原冻结runtime，不用其结果认证新发布源码。

默认清单由 `tools/validation.json` 管理，README不再重复维护第二份测试表。历史比较明确传入数据，缺失标skipped；少量private工具只共享枚举、SHA与报告保存，独立制造解及误差公式不合并。本地完整证据位置及SHA见 EVIDENCE.json；干净克隆不包含这些历史大数据。

发布验收汇总：`results/release_20261024/validation.json`；完整维护清单、独立基线和构建receipt见同目录 `final_validation_01/`，CLI证据为 `cli_validation.json`。首两轮失败原样保留：过粗测试网格被活动域连接守卫拒绝；科学包JSON曾包含稀疏观测矩阵，现改为MAT完整保存、JSON标量摘要。最终17项请求检查通过；radial的完整空间资格仍false，未解析的舍入级差值不被伪作阶数。未重跑昂贵research/order/IVP，未增加稳态、全局空间或长期轨迹资格。
## 2026-11-03：真实中网格割线等值除零

冻结 `bc61655` 中网格从原初值演化，在内存接受态 `tau=12.479430743775856`、step3218之后的样条节点定位失败；最后持久态为 `tau=12.26748666645546`。保留原空间矩阵及该内存态的初始RHS进行单步重放，复现同一 `peak.m` 数组越界：第三次割线更新前后的函数值完全相同，均为 `2.3592239273284576e-16`，第四次分母为零，算得 `h=Inf`，随后RK阶段失败。原输入时间区间为 `0.0033627421368996437`；此前最后有限时间点为 `0.002874308598106941`。

最小修复在现有精确根终止条件中加入两次函数值完全相等的情形，在未定义的除法之前结束迭代，保留最后已计算的trial。没有近零阈值、备用求根、缩步重试、状态择优或回滚，也未修改峰值函数、Poisson、对流或remesh。

独立真实状态回归通过：从该内存态开始，首步正常定位原墙knot263，峰坐标缺陷为 `5.551115123125783e-17`；7步到达 `tau=12.5`，最大stage CFL为 `.2850036192704822`，完整F/Q、尺度和四字时钟均通过原接口检查，输入SHA前后相同。该回归明确使用冻结旧空间矩阵，只隔离候选时间代码，不能冒称新的空间方法长期轨道。

原默认linear5/FEM/RK4时间三层屏在相同问题、步长和阈值下通过：W阶 `3.9624346573`、q阶 `3.9551322852`、Q阶 `3.8109044361`；F差低于原分辨下限，不赋阶。正式封存为 `results/secant_zero_20261103_01/temporal_02`，源码和测试输入未改变。`temporal_01` 的数值测量完成，但因测试文件编辑与运行重叠触发输入变更检查，封存无效，保留为排程失败记录，不计入通过证据。

故障重放、修复重放和真实回归分别位于 `baseline_trace/REPLAY.json`、`PATCH_REPLAY.json`、`step_regression_01/step/step_api.json`。本次只修复完全相等函数值造成的未定义割线，不声称覆盖任意外推、近等值分母或全局分支切换；未重启原失败长实验。

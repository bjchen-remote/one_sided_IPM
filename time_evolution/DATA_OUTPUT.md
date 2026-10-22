# 原生 MATLAB 数据输出

保存结果的正式运行会直接生成 `data.mat` 和 `data_frames/`。接受态缓存已包含完整 RHS、流函数、速度和规范量；输出模块只读这些值，不重建 Poisson，不推进时间。旧 checkpoint 的 schema、配置和数值内容保持兼容。

```matlab
result = ipm.solve(options);                 % saveResults=true 时默认开启
data = ipm.output.loadData(result.metadata.datasetFile);
plot(data.tau,data.cr)
semilogy(data.tau,data.residual_inf)          % 全域残差最大值范数
data = ipm.output.selectTime(data,3.5);       % 最近保存帧，不相等会警告
mesh(data.x,data.y,data.rho_xt)
quiver(data.x,data.y,data.vx,data.vy)         % 选中时刻的合成输运速度
```

严格选时用 `selectTime(data,tau,'exact')`，按帧号用 `selectTime(data,k,'index')`。`ipm_select_time` 是同一函数的短入口。未选时的 `data.mat` 是标量历史与帧索引；`loadData` 自动读取最后一帧，普通 `load(...,'data')` 后也可显式选帧。

## 文件与结构

```text
run/
  data.mat                   普通 data 结构体：标量、诊断、帧索引、元数据
  data.mat.lock              操作系统文件锁；进程退出后锁自动释放
  data_frames/
    frame_step..._tau....mat  一个不可变 frame 结构体及 journal 元数据
  result...mat               原有 result 格式
  checkpoints/               原有 schema-4 恢复格式
```

`data.tau/t/step/cl/co/cr/Cx/Cy/Comega/Xshift/M/fullRateNorm/coreRateNorm` 是对齐的标量列向量。`tau` 为规范时间，`t` 为物理时间。原始公共量通过 `data.meta.sourceCommonNames` 映射到权威字段，其余量只在 `data.diagnostics.common` 保存一次。网格、规范、壁面记录分别位于 `data.diagnostics.mesh/gauge/wallCore`。质量掩码为 `trustedMask/continuousTrustedPrefix`；记录不会因质量标记而静默删除。

`data.frames(k)` 是包含时间、步数、历史行号、原生尺寸、文件名和 SHA256 的索引。`ipm.output.loadFrame(data,k)` 只读取该帧。帧内 `rho/w/rho_t/rho_xt/rho_y/psi/psiLow/u1/u2/transport_u1/transport_u2` 均为 `[y,x]`，每帧带自己的 `x/y` 与缩放因子，允许节点数改变。它不要求所有帧拼成一个三维数组。

`rho=R`，`w=R_X`，`rho_xt=∂τ(R_X)=acceptedRhoRate*Dx'`，保持计算坐标 X、Y 不变。它就是完整域收敛残差，不是相邻输出帧差分，也不是物理时间的 `∂t(rho_x)`。物理坐标和物理场在选帧时从缩放因子精确换算，不另外保存全部物理二维数组。流函数的低分量 `psiLow` 单独保存，选帧后名为 `psi_low`。

`data.residual_inf` 是与 `data.tau` 对齐的列向量，等于完整原生网格（含边界）上的 `max(abs(rho_xt),[],'all')`。它是绝对、无权最大范数，与原字段 `fullRateNorm` 数值相同；缺测保留 NaN。新维护版 writer 和打包产物保存此字段，`loadData` 或选帧接口也会为已有文件自动提供该字段，不改写原文件。已经冻结的 source24 runtime 保持原样，其原生文件中的范数仍保存为 `fullRateNorm`。

选帧后的 `data.vx/data.vy` 是已保存的合成输运速度，直接引用原生帧的 `transport_u1/transport_u2`；`data.vmag=hypot(data.vx,data.vy)`。它们位于计算坐标 X/Y，使用规范时间 τ，公式为 `Vx=U1+cx*X+cr`、`Vy=U2+cy*Y`；当前等向 source24 中 `cx=cy=cl`。流体速度仍为 `u1/u2`。原始二维合成速度已经逐帧保存，无需重复写数组；这些绘图字段从保存值读取，不重建内部中心化 RHS 的临时数组。

原配置核心宽度使用 `adaptiveLevels(end)`；当前生产值为 0.5，`wxCore/wyCore` 不标成 90% 宽度。完整壁面及二维帧足以在分析时求任意高度宽度、坍缩和尾部拟合。拟合与最终质量位于 `data.analysis`；不同规范、人工迭代和历史报告应保存为独立数据集。

## 保存策略与恢复

可选第三参数独立于 checkpoint 配置，不改变已有配置解析或签名合同：

```matlab
policy = struct('enabled',true,'file',fullfile(newRun,'data.mat'),'frameEvery',0.5);
result = ipm.solve(overrides,checkpoint,policy);
```

只允许 `enabled/file/frameEvery` 三个策略字段。默认启用取决于 `saveResults`，帧间隔取有限的 checkpoint 间隔，否则为 0.5τ。显式关闭用 `struct('enabled',false)`。保存初始/恢复状态、到期记录、记录到的 remesh 变化和终态。标量沿用求解器原记录频率，不额外改变数值程序的观察时刻。

帧先写同目录临时文件，再原子更名，最后原子发布 `data.mat`。中断留下的已完成帧可从 journal 恢复；临时文件不会当成完成帧。重复 `(tau,step)` 的密度、网格、缩放、时间及其它场必须逐位一致；只有 `rho_t/rho_xt` 和对应两个范数允许已实测的窄 RHS 重建舍入差（每字段 128 个机器精度单位）。原帧和首次范数不更改，差异保存在 `meta.replayChecks`。此规则只处理重建兼容性，不改变任何 PDE 接受或收敛阈值。运行 ID、数值配置与数值源码指纹不同就拒绝合并。同一文件只允许一个 writer，操作系统在进程崩溃后自动释放锁；锁文件可以保留。

续算到新目录时继承源 checkpoint 附近匹配的数据集索引；继承的帧可引用先前运行目录，不重复复制。只继承 checkpoint 历史内的帧，不接入未来分支。对开启原生输出以前的 checkpoint，旧标量历史仍完整保留，但未记录的残差或空间帧不会伪造；缺测残差为 NaN。输出失败会保留可保存的接受态 checkpoint 并报错，不宣称数据保存成功。

日常输出只重写标量目录，不重复写全量二维历史。标量目录仍随记录数增长；这次修改不改变原求解器累计 history 的内存策略。

已有核验过的 schema-2 packed 档案可在实际数值 runtime 下执行 `stream=ipm.output.dataset('seed',packedFile,policy); clear stream`，建立可追加 journal。该操作逐位复制已有原生帧并导入已记录残差，不复算 PDE；原档案保持只读。之后从同一轨迹 checkpoint 续算即可保留完整前史。

## 单文件归档

```matlab
ipm.output.packData('run/data.mat','archive.mat');
data = ipm.output.loadData('archive.mat');
```

归档把各帧 payload 放入 `data.frames(k).payload`，仍用同一选时接口，单独移动归档文件即可绘图。打包会使用与归档大小成比例的内存，作为显式交付操作，不放在时间推进循环中。移动持续运行的数据集时，应连同帧目录一起移动；跨续算目录引用也须保留，或先打包。

`tests/check_native_dataset.m` 验证真实演化与恢复的输出开关对照、范数、重复帧、journal 恢复、错误缓存/其它运行拒绝、锁、可变网格 I/O 和打包绘图。可变尺寸夹具只证明 I/O 合同，不声称新网格上的数值精度。

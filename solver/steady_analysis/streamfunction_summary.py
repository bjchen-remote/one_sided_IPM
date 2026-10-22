"""Render the streamfunction geometry/shape observations and their limitations."""
from pathlib import Path
import json
import os

os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/ipm_streamfunction_mpl')
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/streamfunction_geometry_20260928'


def shape_figure(data):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    a=np.load(OUT/'shape_fit.npz'); xi=a['xi']
    observed=a['sinh_observed_D_over_B']
    single=a['sinh_student_prediction_D_over_B']; rich=a['sinh_student_skew_even_prediction_D_over_B']
    fig,axes=plt.subplots(1,3,figsize=(15,4.6),layout='constrained')
    colors=plt.cm.viridis(np.linspace(.1,.9,7))
    for i,Y in enumerate(a['nativeY']):
        axes[0].plot(xi,observed[i],color=colors[i],alpha=.75,lw=1)
    axes[0].plot(xi,single,'--',color='#df7f0e',lw=2,label='1-parameter mother')
    axes[0].plot(xi,rich,color='#1f2937',lw=1.8,label='3-parameter mother')
    axes[0].set(xlabel=r'$\xi$',ylabel=r'$D/B$',title='Shared slope shape at seven heights');axes[0].legend(fontsize=8)
    for i in (1,3,5):
        axes[1].plot(xi,rich-observed[i],color=colors[i],label=f'Y={a["nativeY"][i]:g}')
    axes[1].axhline(0,color='#aaa',lw=.7)
    axes[1].set(xlabel=r'$\xi$',ylabel=r'$(D_{fit}-D_{obs})/B$',title='Held-out slope residual (3 parameters)');axes[1].legend(fontsize=8)
    for i in (1,3,5):
        g=data['geometry'][i]; scale=g['yInner']*g['B']*(g['L']+g['R'])/2
        axes[2].plot(xi,a['sinh_observed_psi'][i]/scale,color=colors[i],label=f'Y={g["Ynative"]:g}')
        axes[2].plot(xi,a['sinh_student_skew_even_prediction_psi'][i]/scale,'--',color=colors[i])
    axes[2].set(xlabel=r'$\xi$',ylabel=r'$\psi/[yB(L+R)/2]$',title='Streamfunction: observed / dashed fit');axes[2].legend(fontsize=8)
    fig.suptitle('Conditional shape fit | height geometry is measured at every height',fontsize=14)
    fig.savefig(OUT/'streamfunction_shape_fit.png',dpi=180);plt.close(fig)


def main():
    geom=json.loads((OUT/'geometry.json').read_text()); data=json.loads((OUT/'shape_fit.json').read_text())
    records=geom['records']; r20=records['native_tau20']
    lines=[
        '# 流函数的重点几何：零线、速度峰线与共同母形',
        '',
        '2026-09-28。已分析原生 τ≈10、15、20 三帧和四组条件有限箱候选。最明确的结果是：核心流函数横向单调，应该观测 ψ/y 的零线及其横向导数 u₂/y 的单峰几何；后者确实存在可压缩的共同母形。它的偏斜方向与 ρ_x 相反。当前拟合能很好重构 ψ 和速度，但二阶导数恢复 ρ_x 尚未通过相容性及正性检验。',
        '',
        '## 1. 数据与规范',
        '',
        '原生 checkpoint 不保存流函数；本轮用原求解器各作一次静态恢复，限制两线程，没有推进时间。τ20 恢复的 ψ/u₁/u₂ 与原 run_analysis.mat 逐元素完全一致。三帧恢复量保存于 source_frames.npz，内含导出脚本与来源。',
        '',
        r'项目采用 $u=(-\psi_y,\psi_x)$。中心化必须取 $\psi_c=\psi_{saved}+U_0y$，其中 $U_0=u_{1,saved}(0,0)$，才能得到 $u_{1,c}=u_{1,saved}-U_0$。条件有限箱根也各自重新中心化，不能直接使用保存的 ψ 做零线比较。',
        '',
        r'原生坐标与内坐标关系为 $X=a+\ell x,\ Y=\ell y$；$\psi_{native}=P\ell^2\psi_{inner}$。本文中未加 native 标记的公式使用内坐标，表格明确注明坐标。原生坐标是保存的重标度场坐标，不是未经缩放的物理坐标。',
        '',
        '这四个小箱候选均是冻结外部响应的条件根，不是已经认证的全半平面稳态；不同箱之间的差异作为观测结果保留。',
        '',
        '## 2. 应该拟合哪个量',
        '',
        r'定义 $H=\psi_c/y$，壁面连续延拓为 $H(x,0)=-u_1(x,0)$；再定义 $D=H_x=u_2/y$。准确恒等式为',
        '',
        r'\[u_1=-H-yH_y,\qquad u_2=yH_x,\qquad q=\rho_x=-2H_y-y(H_{xx}+H_{yy}).\]',
        '',
        r'因此 $H_y(x,0)=-q(x,0)/2$ 是骨架必须满足的壁面条件。不能把除以 y 当作一个任意归一化。',
        '',
        r'七组数据在共同核心 $|x|\le3,0\le y\le1$ 的密采样上均有 $D>0$。这说明 H 横向递增，在这个核心内没有截面峰谷，且 y>0 时没有中心化 Darcy 速度停滞点。这里不作全域拓扑判断。',
        '',
        '保留的几何是：H 的核心零线 Z(y)；D 的峰线 n(y)、峰高 B(y)、左右半宽 L(y),R(y)。其中 H=0 是 ψ=0 的内部流线，D 的峰线是每个高度上 u₂ 的峰线。ψ 的等值线属于中心化 Darcy 场；有效输运场 u+λ(x,y)+(v,0) 的散度为 2λ，不能把它的轨线与这些等值线混称。',
        '',
        f'![流函数几何]({OUT}/streamfunction_geometry.png)',
        '',
        '图中紫线为 ψ=0，绿虚线为 q 峰，粉线为 D 峰，橙虚线为 u₁=0。流函数零线在 q 峰右侧，D 峰却明显在左侧。',
        '',
        '## 3. 直接观测的几何',
        '',
        'τ20，下表全部使用原生坐标；L_D,R_D 是 D 峰到左右半高点的距离。',
        '',
        '| Y | ψ 零线 Z | q 峰 m | D 峰 n | B/A_q | L_D/R_D |',
        '|---:|---:|---:|---:|---:|---:|',
    ]
    for r in r20['nativeRows']:
        if r['Y']>=1-1e-9:
            n=r['native']; w=r['D']['widths']['0.5']
            lines.append(f"| {r['Y']:g} | {n['psiZero']:.5g} | {n['q']['m']:.5g} | {n['D']['m']:.5g} | {r['ratios']['DpeakOverQpeak']:.5f} | {w['left']/w['right']:.4f} |")
    lines += ['',
        '在 Y=4–256，D 的左半宽约是右半宽的 2.24–2.34 倍；q 则是右侧更宽。采用同一个 β=3/2 的双曲坐标，D 的偏斜参数约为 −0.60 至 −0.63，q 的偏斜参数为正。因此不能把两场的峰线、宽度和偏斜参数直接视为同一组数。',
        '',
        '同一原生高度 Y=64 的时间对照：',
        '',
        '| τ | Z/Y | n/Y | B/A_q | k_D（β=3/2） |',
        '|---:|---:|---:|---:|---:|']
    for tau in (10,15,20):
        r=next(v for v in records[f'native_tau{tau}']['nativeRows'] if abs(v['Y']-64)<1e-9)
        lines.append(f"| {tau} | {r['native']['psiZero']/64:.5f} | {r['native']['D']['m']/64:.5f} | {r['ratios']['DpeakOverQpeak']:.5f} | {r['D']['k']:.5f} |")
    lines += ['',
        'D 峰线的位置随时间变化相对较小；ψ 零线、振幅比与偏斜仍有明显变化。归一化母形较稳不意味着几何参数已经收敛。',
        '',f'![时间与高度变化]({OUT}/streamfunction_parameters.png)','',
        '四组条件根在内坐标 y=1 的比较：','',
        '| 条件根 | Z | q 峰 | D 峰 | B/A_q |','|---|---:|---:|---:|---:|']
    for name in ('L4_129x65','L4_257x129','L6_193x97','L8_257x129_continued'):
        r=next(v for v in records[name]['innerRows'] if v['y']==1.)
        lines.append(f"| {name} | {r['psiZero']:.5f} | {r['q']['m']:.5f} | {r['D']['m']:.5f} | {r['ratios']['DpeakOverQpeak']:.5f} |")
    lines += ['',
        '核心几何排序一致，但数值尚有显著箱与网格依赖。D 的左翼很宽，小箱内常测不到左半高点，缺失值按 null 保存，没有外推。个别箱边界值高于核心局部峰；程序保留核心局部峰并标记 endpointExceedsCorePeak，不能把它误称为全箱最大值。',
        '',
        'D 的几何主值由保存的 u₂/y 插值提取；另对 ψ 样条求导作位置敏感性对照。L8 在 y=1 的两种峰位相差约 5.8×10⁻⁵，峰高差约 1.2×10⁻⁷。多位峰位数字用于复算，不代表同等物理精度。',
        '', '## 4. 已拟合的横向母形', '',
        r'令 $x=n+WG_k(\xi)$，其中 $G_k=[\sinh(\beta\xi)+k(\cosh(\beta\xi)-1)]/\beta$、$\beta=3/2$，W,k 由 D 的左右半宽确定。拟合 $D/B=g(\xi)$。τ20 用 Y=4,16,64,256 训练，用 Y=8,32,128 留出，统一 $|\xi|\le1.5$。',
        '',
        '**每个高度（包括留出高度）的峰位、峰高、半宽和零线都来自实测；只留出了母形，不是高度函数预测。** 误差按均匀 ξ 采样，ψ 在各高度先除以 yB(L+R)/2 再汇总，防止大高度占据全部权重。',
        '',
        r'最简单的预设正母形为 $g_p(\xi)=[1+(2^{1/p}-1)\xi^2]^{-p}$，实测拟合 $p=1.7576$。它自动满足 $g(0)=1,g^{\prime}(0)=0,g(\pm1)=1/2$。',
        '',
        r'三参数对照在其上乘 $\exp[a\xi^3(\xi^2-1)/(1+\xi^2)^{5/2}+b\xi^2(\xi^2-1)/(1+\xi^2)^2]$，参数 $(p,a,b)=(0.96689,-0.13235,-0.30257)$。它保持上述标定点并在拟合区间内正、单峰；未据此声称全域尾部适定。',
        '', '| 母形 | 留出 D / u₂ 相对误差 | 留出 ψ 相对误差 |','|---|---:|---:|']
    for name,label in [('student','单参数'),('student_skew','双参数'),('student_skew_even','三参数')]:
        f=next(z for z in data['fits'] if z['coordinate']=='sinh' and z['family']==name)
        lines.append(f"| {label} | {100*f['shapeHoldout']['relativeL2']:.3f}% | {100*f['psiHoldoutHeightNormalized']['relativeL2']:.3f}% |")
    lines += ['',
        '此外测了均宽、左右分别缩放、双曲三种对齐方式及三个预设母形，共九组结果全部保留。不同坐标覆盖的物理横向带不同，不能只用其误差排名宣称某一坐标物理上最优。',
        '',f'![母形及积分拟合]({OUT}/streamfunction_shape_fit.png)','',
        '冻结 τ20 母形，不在 τ10/15 重新拟合；每个时刻仍使用其自己的实测几何，覆盖相同的七个原生高度：','',
        '| 母形 | τ | D / u₂ 误差 | ψ 误差 |','|---|---:|---:|---:|']
    for fr in data['temporalTransfer']['frames']:
        for m in fr['models']:
            label='单参数' if m['family']=='student' else '三参数'
            lines.append(f"| {label} | {fr['tau']} | {100*m['shapeAllHeights']['relativeL2']:.3f}% | {100*m['psiAllHeightsHeightNormalized']['relativeL2']:.3f}% |")
    lines += ['',
        '单参数母形的跨时刻迁移更稳；三参数主要提高 τ20 的当前精度。这里检验的是形状迁移，仍未预测几何随时间的演化。',
        '', '## 5. 导数验证：小流函数误差不等于正确的密度梯度', '',
        '利用相邻高度重新测得的几何，在固定 x 下对拟合 ψ 作五点高度差分，同时计算横向二阶导数。以下取 h/Y=0.005；完整 JSON 还保留 0.01、0.02 的步长对照。',
        '', '| 母形 | Y | u₁ 误差 | −Δψ 对 q 误差 | q_fit 最小值 | 负值点数 / 301 |',
        '|---|---:|---:|---:|---:|---:|']
    for key,label in [('singleParameterDerivativeAudit','单参数'),('derivativeAudit','三参数')]:
        for row in data[key]['rows']:
            s=row['stepChecks'][0]
            lines.append(f"| {label} | {row['Ynative']:g} | {100*s['u1']['relativeL2']:.3f}% | {100*s['q']['relativeL2']:.3f}% | {s['qfitMinimum']:.5g} | {s['qfitNegativeSamples']} |")
    lines += ['',
        'Y=128 保存 q 的最小值约 +8.66×10⁻⁴。单参数模型在三种步长下均出现负 q；三参数模型的负值较小且随步长改变符号，因此其二阶导数尚不稳定，不能认证正性。原始 ψ 样条在相同横向采样上的 q 误差仅约 0.029%–0.055%，远小于上述几个百分点。不能用 ψ 的高拟合精度替代二阶相容性检验。',
        '',
        '主几何脚本另在共同内区 |x|≤2、0.05≤y≤1 比较源场样条导数与独立保存速度；原生速度误差小于 8×10⁻⁷、q 误差约 4×10⁻⁵。壁面一侧样条误差较大，已单独保存。所有这些都是插值/离散一致性审核，不是连续解精度认证。',
        '', '## 6. 几何高度律与近壁约束', '',
        '作为有限区间描述，另对 τ20 原生 Y=4–256 的几何拟合 cY^p；使用同样的四个训练高度、三个留出高度：','',
        '| 几何 | c | p | 留出相对 L2 | 留出最大相对差 |','|---|---:|---:|---:|---:|']
    for key,label in [('psiZero','零线 Z'),('Dpeak','峰位 n'),('Damplitude','峰高 B'),('DleftHalfWidth','左半宽 L'),('DrightHalfWidth','右半宽 R')]:
        z=r20['powerFits']['fits'][key]
        lines.append(f"| {label} | {z['coefficient']:.6g} | {z['exponent']:.6f} | {100*z['heldoutRelativeL2']:.3f}% | {100*z['heldoutMaxRelative']:.3f}% |")
    lines += ['',
        '这些是局部经验幂，不是远场指数，也不能把这些高度拟合与上面的母形误差直接拼成完整模型精度：前面的母形检验仍使用实测高度几何。峰位 n 的最大相对差约 10%，尤其不适合被宣称为准确幂律。',
        '',
        r'壁面原点令 $s=H_x(0,0)$。由 Poisson 和无穿透，必须有 $\psi=sxy-q(0,0)y^2/2+O((|x|+y)^3)$，所以',
        '',
        r'\[Z_\psi^{\prime}(0)=\frac{q(0,0)}{2s},\qquad Z_{u_1}^{\prime}(0)=\frac{q(0,0)}s=2Z_\psi^{\prime}(0).\]',
        '',
        f"τ20 测得 s={r20['wallGeometry']['s']:.6f}、q(0,0)≈1，预测两条零线斜率分别为 {r20['wallGeometry']['predictedPsiZeroSlope']:.6f} 和 {r20['wallGeometry']['predictedU1ZeroSlope']:.6f}。这是近壁拟合可直接使用的约束。若进一步为严格稳态且使用现有峰规范，s=λ−μ；本轮没有把有限时刻数据当作严格稳态。",
        '',
        '## 7. 可保留的流函数骨架',
        '',
        r'将 $J(\xi;k)=\int_0^\xi g(t)G_k^{\prime}(t)dt$，$\xi_0(y)=G_k^{-1}((Z-n)/W)$，则观测支持的表达式是',
        '',
        r'\[\boxed{\psi(x,y)=yB(y)W(y)\big[J(\xi;k(y))-J(\xi_0(y);k(y))\big],\quad x=n(y)+W(y)G_{k(y)}(\xi).}\]',
        '',
        '这保留一条共同母形和五条几何函数 n、B、W、k、Z；它自动给出 u₂=yBg(ξ) 和指定的 ψ=0 零线。在这些函数壁面正则的前提下也自动满足无穿透。当前数据不支持把 Z 与 n 合并，也不支持把速度几何与 q 几何直接等同。',
        '',
        '后续筛选应让该骨架的 −Δψ 回生成实测 q 的峰线、峰高和偏斜，并首先强制 H_y(x,0)=−q(x,0)/2。只有再通过完整输运与全域速度选取，才能称为共同闭合骨架。本轮完成的是观测、有限区间拟合、导数审核与缺口定位。',
        '', '## 复现与文件', '',
        f'- [streamfunction_sources.py]({ROOT}/solver/steady_analysis/streamfunction_sources.py)：中心化原生与条件数据；[source_frames.npz]({OUT}/source_frames.npz) 内嵌静态导出脚本。',
        f'- [streamfunction_geometry.py]({ROOT}/solver/steady_analysis/streamfunction_geometry.py)：几何、近壁与源场导数审核；输出 [geometry.json]({OUT}/geometry.json)、两幅几何图。',
        f'- [streamfunction_shape_fit.py]({ROOT}/solver/steady_analysis/streamfunction_shape_fit.py)：九组固定族拟合、留出、时间迁移与两组导数审核；输出 [shape_fit.json]({OUT}/shape_fit.json) 和 [shape_fit.npz]({OUT}/shape_fit.npz)。',
        f'- [streamfunction_summary.py]({ROOT}/solver/steady_analysis/streamfunction_summary.py)：从 JSON/NPZ 生成本报告及母形图。',
        '',
        '运行顺序：geometry → shape_fit → summary。依赖 numpy/scipy/matplotlib；本轮使用已有 /private/tmp/ipm_steady_deps，限制 BLAS/OpenMP 两线程。原数据、求解器和旧报告未修改。',
    ]
    (OUT/'REPORT_ZH.md').write_text('\n'.join(lines)+'\n')
    shape_figure(data)
    print(OUT/'REPORT_ZH.md')


if __name__=='__main__':
    main()

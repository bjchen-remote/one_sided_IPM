"""Fresh inverse-H geometry from centred psi only; no inherited ridge model.

H=psi/y is only the exact regular wall quantity. All sampled heights are
positive. Horizontal profiles are never aligned to D peaks or fitted using
half widths, sinh coordinates, or any prescribed transverse mother shape.
"""
from pathlib import Path
import json
import numpy as np
from scipy.interpolate import RectBivariateSpline, PchipInterpolator
from scipy.optimize import brentq

from streamfunction_sources import load_native, load_conditional, CONDITIONAL_NAMES

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'runs/streamfunction_fresh_20260928'
HEIGHTS = np.array([4., 8., 16., 32., 64., 128., 256.])
CORE_HEIGHTS = np.array([.0625, .125, .25, .5, 1.])
ETA = np.linspace(-1., 1., 201)
TRAIN = np.array([True, False, True, False, True, False, True])


def relative(error, reference):
    return float(np.linalg.norm(error) / np.linalg.norm(reference))


class Observe:
    def __init__(self, profile):
        self.p = profile
        # No q or stored velocity is used in any observation below.
        self.psi = RectBivariateSpline(profile['y'], profile['x'], profile['psi'])

    def H(self, x, y):
        x = np.asarray(x, dtype=float)
        if y <= 0 or y > self.p['y'][-1] or np.min(x) < self.p['x'][0] or np.max(x) > self.p['x'][-1]:
            raise ValueError('Queries must remain inside saved positive-height support')
        return self.psi.ev(np.full_like(x, y), x) / y

    def row(self, y, radius, samples=16001):
        p = self.p
        xx = np.linspace(max(p['x'][0], -8*radius), min(p['x'][-1], 8*radius), samples)
        hh = self.H(xx, y)
        changes = np.flatnonzero(hh[:-1]*hh[1:] < 0)
        zeros = [brentq(lambda x: float(self.H(x, y)), xx[i], xx[i+1]) for i in changes]
        zeros.extend(xx[hh == 0].tolist())
        if not zeros:
            raise ValueError('No zero in the inspection window')
        z = float(min(zeros, key=abs))
        lo, hi = max(p['x'][0], z-4*radius), min(p['x'][-1], z+4*radius)
        support = np.linspace(lo, hi, samples)
        values = self.H(support, y)
        if not np.all(np.diff(values) > 0):
            raise ValueError('H is not monotone on this measured inversion window')
        hminus, hplus = float(self.H(z-radius, y)), float(self.H(z+radius, y))
        amplitude = (hplus-hminus)/2
        levels = amplitude*ETA
        if levels.min() < values[0] or levels.max() > values[-1]:
            raise ValueError('Requested level leaves the observed monotone branch')
        inverse = PchipInterpolator(values, support, extrapolate=False)
        displacement = (inverse(levels)-z)/radius
        return dict(yInner=float(y), nativeY=float(p['ell']*y), radiusInner=float(radius),
                    zeroInner=z, zeroNative=float(p['a']+p['ell']*z),
                    HscaleInner=amplitude, HscaleNative=float(p['P']*p['ell']*amplitude),
                    HAtProbeMinus=hminus, HAtProbePlus=hplus,
                    inversionSupportInner=[float(lo), float(hi)],
                    allZerosInInitialWindow=zeros, minSampledHIncrement=float(np.diff(values).min()),
                    normalizedInverseDisplacement=displacement.tolist())


def matrix(rows):
    return np.array([r['normalizedInverseDisplacement'] for r in rows])


def analyze(rows, train=None):
    data = matrix(rows)
    u, s, vh = np.linalg.svd(data, full_matrices=False)
    slopes = data@ETA/(ETA@ETA)
    linear = slopes[:, None]*ETA
    odd = (data-data[:, ::-1])/2
    result = dict(singularValues=s.tolist(),
                  rankRelativeL2={str(k): float(np.linalg.norm(s[k:])/np.linalg.norm(s)) for k in (1, 2, 3)},
                  bestLinearLevelModelRelativeL2=relative(data-linear, data),
                  evenPartRelativeL2=relative(data-odd, data),
                  learnedModes=vh[:3].tolist(), learnedCoefficients=(u[:, :3]*s[:3]).tolist())
    if train is not None:
        basis = np.linalg.svd(data[train], full_matrices=False)[2]
        result['trainingRows'] = np.flatnonzero(train).tolist()
        result['heldoutSubspaceReconstruction'] = {
            str(k): relative(data[~train]@basis[:k].T@basis[:k]-data[~train], data[~train]) for k in (1, 2, 3)}
    return result


def transfer(training_rows, rows, train=None):
    data = matrix(training_rows)
    if train is not None:
        data = data[train]
    modes = np.linalg.svd(data, full_matrices=False)[2]
    target = matrix(rows)
    return {str(k): relative(target@modes[:k].T@modes[:k]-target, target) for k in (1, 2, 3)}


def main():
    OUT.mkdir(exist_ok=True)
    natives = {tau: Observe(load_native(tau)) for tau in (10, 15, 20)}
    outer, core, sensitivity = {}, {}, {}
    for tau, observer in natives.items():
        ell = observer.p['ell']
        rows = [observer.row(Y/ell, Y/ell) for Y in HEIGHTS]
        outer[str(tau)] = dict(rows=rows, analysis=analyze(rows, TRAIN))
        rows = [observer.row(y, 1.) for y in CORE_HEIGHTS]
        core[f'native_tau{tau}'] = dict(rows=rows, analysis=analyze(rows))
    for name in CONDITIONAL_NAMES:
        observer = Observe(load_conditional(name))
        rows = [observer.row(y, 1.) for y in CORE_HEIGHTS]
        core[name] = dict(rows=rows, analysis=analyze(rows))
    obs = natives[20]
    for factor in (.5, 1., 2.):
        rows = [obs.row(Y/obs.p['ell'], factor*Y/obs.p['ell']) for Y in HEIGHTS]
        sensitivity[str(factor)] = dict(analysis=analyze(rows), rows=rows)
    outer_transfer = {tau: transfer(outer['20']['rows'], item['rows'], TRAIN) for tau, item in outer.items()}
    core_transfer = {name: transfer(core['native_tau20']['rows'], item['rows']) for name, item in core.items()}
    checks = []
    for observer, y, radius in [(natives[20], 256/natives[20].p['ell'], 256/natives[20].p['ell']),
                                (Observe(load_conditional(CONDITIONAL_NAMES[-1])), 1., 1.)]:
        fine, coarse = observer.row(y, radius), observer.row(y, radius, 8001)
        err = np.array(fine['normalizedInverseDisplacement'])-coarse['normalizedInverseDisplacement']
        checks.append(dict(source=observer.p['name'], maxNormalizedContourDifference=float(abs(err).max())))
    output = dict(protocol=dict(
        source='centred psi from streamfunction_sources; q, D peaks, widths and sinh fits are not used',
        H='psi/y at positive saved heights; exact regular variable, not a proposed shape',
        zero='H=0 root nearest x=0 in the stated inspection window',
        Hscale='S=[H(Z+R,y)-H(Z-R,y)]/2, a secant observable from psi',
        inverse='F_y(eta)=[X(S(y)*eta,y)-Z(y)]/R(y), eta in [-1,1]',
        probeRadius='outer: native R=Y; near wall: inner R=1; these are inspection windows, not fitted width laws',
        sensitivity='outer tau20 repeated with R=.5Y and R=2Y',
        model='SVD modes learned directly from inverse-contour samples; no analytic mother function',
        norms='equal weight per sampled height and eta; rank means linear algebra dimension after explicitly measured normalizations',
        validation='heldout and transfer scores project each full target row onto frozen learned modes; geometry, normalization and modal coefficients remain observed, not predicted',
        scope='finite-time saved fields and conditional finite-box roots; no global or steady-state certification',
        interpolation='bicubic psi, uniform support with no grid union, monotone PCHIP inverse; no extrapolation'),
        eta=ETA.tolist(), nativeHeights=HEIGHTS.tolist(), coreHeights=CORE_HEIGHTS.tolist(),
        outer=outer, core=core, probeRadiusSensitivity=sensitivity,
        frozenTau20OuterModeTransfer=outer_transfer, frozenNativeTau20CoreModeTransfer=core_transfer,
        interpolationChecks=checks)
    (OUT/'level_geometry.json').write_text(json.dumps(output, indent=2, allow_nan=False)+'\n')
    write_report(output)
    print(json.dumps(dict(outer={k:v['analysis']['rankRelativeL2'] for k,v in outer.items()},
                          outerTransfer=outer_transfer, coreTransfer=core_transfer,
                          interpolationChecks=checks), indent=2))


def pct(value):
    return f'{100*value:.3f}%'


def write_report(out):
    lines = ['# 流函数等值线的重新观测', '',
        '本轮只从原始中心化 ψ 取值出发；没有使用 rho_x 骨架、D 峰、左右半高宽或 sinh 坐标。', '',
        'H=ψ/y 只用于消除墙面的必然零因子。正高度上 H 的等值线就是同一截面内对应 ψ 等值线；跨高度固定 H 不等于固定 ψ。以下研究的是 H 等值线族，不能直接称作固定 ψ 的流线族。', '',
        '取零线 Z(y)，直接观测逆映射 X(H,y)。为比较可见范围，给定探针半径 R，用 ψ 的取值定义 S=[H(Z+R)−H(Z−R)]/2，然后比较 F_y(η)=[X(Sη,y)−Z]/R，η∈[−1,1]。外层检查 R=Y；近壁检查固定内坐标 R=1。R 是检查窗口，不是拟合所得宽度或新骨架假设。', '',
        '每个逆函数都从已存储的单调区间内插得到；不外推。SVD 的横向模态完全由这些曲线学习。', '',
        '## 外层等值线族：少量形变比单个固定母形更准确', '',
        '| τ | rank 1 误差 | rank 2 误差 | rank 2 留出子空间误差 | 最佳线性等值线模型误差 |',
        '|---|---:|---:|---:|---:|']
    for tau,item in out['outer'].items():
        a=item['analysis']
        lines.append(f"| {tau} | {pct(a['rankRelativeL2']['1'])} | {pct(a['rankRelativeL2']['2'])} | {pct(a['heldoutSubspaceReconstruction']['2'])} | {pct(a['bestLinearLevelModelRelativeL2'])} |")
    lines += ['', '原始高度为 Y=4,8,16,32,64,128,256。留出检验只用 Y=4,16,64,256 学习横向模态，再把其余整条已观测曲线投影到该子空间；它检验模态是否可复用，没有预测留出高度的 Z、S 或模态系数。', '',
              'tau20 的等值线明显不对称且非线性；低秩现象不是把一个微小线性邻域放大后的必然结果。rank 1 已较简洁，但数据保留系统形变，rank 2 显著改善。', '',
              '| 冻结 tau20 训练模态，投影到 | rank 1 | rank 2 |', '|---|---:|---:|']
    for tau,score in out['frozenTau20OuterModeTransfer'].items():
        lines.append(f"| τ={tau} | {pct(score['1'])} | {pct(score['2'])} |")
    lines += ['', '## 零线与窗口敏感性', '',
              '| τ | Z/Y，Y=4 | Z/Y，Y=256 |', '|---|---:|---:|']
    for tau,item in out['outer'].items():
        first,last=item['rows'][0],item['rows'][-1]
        lines.append(f"| {tau} | {first['zeroNative']/first['nativeY']:.6f} | {last['zeroNative']/last['nativeY']:.6f} |")
    lines += ['', 'Z/Y 随高度和时间变化，数据不支持把零线直接固定成一条跨时间通用直线。', '',
              '| tau20 探针 R/Y | rank 1 | rank 2 |', '|---|---:|---:|']
    for factor,item in out['probeRadiusSensitivity'].items():
        a=item['analysis']['rankRelativeL2'];lines.append(f"| {factor} | {pct(a['1'])} | {pct(a['2'])} |")
    lines += ['', '改变观察半径后，两个数据模态仍有效；误差会随范围变化，不能把有限范围塌缩解释为全局恒等式。', '',
              '## 近壁与条件根交叉检查', '',
              '内坐标 y=.0625,.125,.25,.5,1，固定探针 R=1。条件根保留其冻结外部源的限定。', '',
              '| 数据 | rank 1 | rank 2 | 投影到 native τ20 两模态 |', '|---|---:|---:|---:|']
    for name,item in out['core'].items():
        a=item['analysis']['rankRelativeL2'];score=out['frozenNativeTau20CoreModeTransfer'][name]['2']
        lines.append(f"| {name} | {pct(a['1'])} | {pct(a['2'])} | {pct(score)} |")
    lines += ['', '因此，目前直接观测支持的是“零线位置 + 直接读数尺度 + 两个数据学得的等值线形变模态”。这只是可检验的描述性候选，尚不是固定解析母形，也没有确定高度函数的闭合方程。', '',
              '数据反对：未检验就采用旧密度骨架；把流函数设成单峰；把固定母形、直零线或全局等值线自相似当成前提。这里没有排除其他、将来更优的非线性配准可能达到单模态；当前结果只比较明确记录的直接读数规范。', '',
              '## 复现与数值检查', '',
              '脚本：`solver/steady_analysis/streamfunction_fresh_levels.py`。完整数据、模态、系数、支持区间和规范保存在 `level_geometry.json`。', '',
              '把逆插值支持采样从 16001 减为 8001，所测最大归一等值线位移差：']
    for check in out['interpolationChecks']:
        lines.append(f"- {check['source']}: {check['maxNormalizedContourDifference']:.3e}")
    (OUT/'LEVEL_GEOMETRY.md').write_text('\n'.join(lines)+'\n')


if __name__ == '__main__':
    main()

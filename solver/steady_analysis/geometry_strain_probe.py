"""Read-only finite-patch strain audit of the measured sinh-profile geometry.

No PDE evolution, extrapolation, fitted tail, or global closure is used.
"""
from pathlib import Path
import json
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import CubicSpline, PchipInterpolator, RectBivariateSpline
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/geometry_loop_20260928"
BETA = 1.5
REGIONS = ["left_xi_lt_minus_1p5", "core", "right_xi_gt_1p5"]


def sum64(x):
    return float(np.sum(x, dtype=np.longdouble))


def forward(xi, geo):
    m, A, W, k = np.moveaxis(geo, -1, 0)
    return m + W * (np.sinh(BETA * xi) + k * (np.cosh(BETA * xi)-1)) / BETA


def inverse(x, geo):
    m, A, W, k = np.moveaxis(geo, -1, 0)
    return (np.arcsinh((BETA * (x-m)/W+k)/np.sqrt(1-k*k))
            - np.arctanh(k)) / BETA


def main():
    z = np.load(OUT / "tail_master_arrays.npz")
    xi, ys, rows, geometry = (z[k] for k in ("outer_xi", "heights", "rows", "geometry"))
    assert np.all(z["valid"]) and np.all(np.isfinite(rows))
    ref_mask = ys <= 1
    f_ref = np.mean(rows[ref_mask], axis=0)
    ref = CubicSpline(xi, f_ref, extrapolate=False)
    pos = ys > 0
    eta_measured = np.log(ys[pos])
    geo_cubic = CubicSpline(eta_measured, geometry[pos], axis=0, extrapolate=False)
    geo_pchip = PchipInterpolator(eta_measured, geometry[pos], axis=0, extrapolate=False)
    row_field = RectBivariateSpline(eta_measured, xi, rows[pos])
    native = loadmat(ROOT / "runs/steady_20260928_tau20/run_analysis.mat", simplify_cells=True)
    xx, yy, qq = (native[k] for k in ("xFull", "yFull", "rhoX"))
    field = RectBivariateSpline(yy, xx, qq)

    def exact_geometry(y):
        result = []
        for h in y:
            row = field(h, xx)[0]
            cs = CubicSpline(xx, row)
            roots = cs.derivative().roots(extrapolate=False)
            roots = roots[(roots > xx[0]) & (roots < xx[-1])]
            m = float(roots[np.argmax(cs(roots))])
            A = float(cs(m))
            half = CubicSpline(xx, row-A/2).roots(extrapolate=False)
            L = m-half[half < m][-1]
            R = half[half > m][0]-m
            W = (L+R)*BETA/(2*np.sinh(BETA))
            k = (R-L)/(R+L)/np.tanh(BETA/2)
            result.append([m, A, W, k])
        return np.asarray(result)

    def evaluate(name, source, geometry_method, rule, ny, ntheta=0, nxi=0):
        # Split at all saved height knots and reporting-bin boundaries.
        knots = np.unique(np.r_[eta_measured, np.log([1., 64.])])
        z_y, w_y = leggauss(ny)
        eta = (0.5*(knots[:-1,None]+knots[1:,None])
               + 0.5*np.diff(knots)[:,None]*z_y).ravel()
        ew = (0.5*np.diff(knots)[:,None]*w_y).ravel()
        y = np.exp(eta)
        if geometry_method == "cubic_log_y":
            geo = geo_cubic(eta)
        elif geometry_method == "pchip_log_y":
            geo = geo_pchip(eta)
        else:
            geo = exact_geometry(y)
        assert np.all(geo[:,1:3] > 0) and np.all(abs(geo[:,3]) < 1)

        all_integrals = []
        source_differences = []
        for j, (lo, hi) in enumerate([(-3.,-1.5),(-1.5,1.5),(1.5,3.)]):
            if rule == "theta_gauss":
                # x=y*tan(theta), dy=y*deta turns the kernel into
                # -(sin(2*theta)/pi) dtheta deta. Split physical x=0.
                lower = np.arctan(forward(np.full(len(y),lo), geo)/y)
                upper = np.arctan(forward(np.full(len(y),hi), geo)/y)
                cut = np.clip(np.zeros_like(lower), lower, upper)
                left = np.stack([lower,cut], axis=1)
                right = np.stack([cut,upper], axis=1)
                z_t, w_t = leggauss(ntheta)
                theta = 0.5*(left+right)[:,:,None]+0.5*(right-left)[:,:,None]*z_t
                tw = 0.5*(right-left)[:,:,None]*w_t
                x = y[:,None,None]*np.tan(theta)
                geo_b = np.broadcast_to(geo[:,None,None,:], x.shape+(4,))
                xiq = inverse(x, geo_b)
                weight = -np.sin(2*theta)/np.pi*tw*ew[:,None,None]
            else:
                # Independent composite quadrature in the original xi coordinate.
                bounds = np.linspace(lo,hi,round((hi-lo)/.01)+1)
                z_x,w_x = leggauss(nxi)
                xiq0 = (0.5*(bounds[:-1,None]+bounds[1:,None])
                        +0.5*np.diff(bounds)[:,None]*z_x).ravel()
                xw = (0.5*np.diff(bounds)[:,None]*w_x).ravel()
                xiq = np.broadcast_to(xiq0[None,:], (len(y),len(xiq0)))
                geo_b = np.broadcast_to(geo[:,None,:], xiq.shape+(4,))
                x = forward(xiq, geo_b)
                x_xi = geo[:,2,None]*(np.cosh(BETA*xiq)+geo[:,3,None]*np.sinh(BETA*xiq))
                yb = y[:,None]
                weight = -2/np.pi*x*yb/(x*x+yb*yb)**2*x_xi*yb*xw[None,:]*ew[:,None]
            assert np.min(xiq) >= -3-1e-10 and np.max(xiq) <= 3+1e-10
            xiq = np.clip(xiq,-3,3)
            yb = np.broadcast_to(y.reshape((-1,)+(1,)*(x.ndim-1)),x.shape)
            etab = np.broadcast_to(eta.reshape((-1,)+(1,)*(x.ndim-1)),x.shape)
            assert x.min() >= xx[0] and x.max() <= xx[-1]
            assert yb.min() >= yy[0] and yb.max() <= yy[-1]
            q_native = field.ev(yb.ravel(),x.ravel()).reshape(x.shape)
            amplitude = geo[:,1].reshape((-1,)+(1,)*(x.ndim-1))
            q_saved = amplitude*row_field.ev(etab.ravel(),xiq.ravel()).reshape(x.shape)
            q_true = q_native if source == "native_q" else q_saved
            q_comp = amplitude*ref(xiq)
            qs = {"S_true":q_true,"S_compressed":q_comp,"deltaS":q_true-q_comp}
            record = {k:sum64(weight*q) for k,q in qs.items()}
            record["physical_x_negative"] = {k:sum64(np.where(x<0,weight*q,0)) for k,q in qs.items()}
            record["physical_x_positive"] = {k:sum64(np.where(x>=0,weight*q,0)) for k,q in qs.items()}
            record["height_bins"] = {}
            for a,b in [(.025,1.),(1.,64.),(64.,1024.)]:
                sel = (yb>=a)&(yb<=b)
                record["height_bins"][f"{a:g}..{b:g}"] = {k:sum64(np.where(sel,weight*q,0)) for k,q in qs.items()}
            record["linearity_error"] = record["S_true"]-record["S_compressed"]-record["deltaS"]
            all_integrals.append(record)
            source_differences.append(sum64(weight*(q_native-q_saved)))
        totals = {k:sum(r[k] for r in all_integrals) for k in ["S_true","S_compressed","deltaS"]}
        totals["delta_over_true"] = totals["deltaS"]/totals["S_true"]
        totals["tail_deltaS"] = all_integrals[0]["deltaS"]+all_integrals[2]["deltaS"]
        totals["tail_delta_over_total_delta"] = totals["tail_deltaS"]/totals["deltaS"]
        totals["native_minus_saved_strain"] = sum(source_differences)
        return dict(name=name,source=source,geometry=geometry_method,rule=rule,
                    y_gauss_order=ny,theta_gauss_order=ntheta,xi_gauss_order=nxi,
                    y_nodes=len(y),regions=dict(zip(REGIONS,all_integrals)),total=totals)

    specs = [
        ("native_theta_coarse","native_q","cubic_log_y","theta_gauss",4,64,0),
        ("native_theta_fine","native_q","cubic_log_y","theta_gauss",8,192,0),
        ("native_xi_composite","native_q","cubic_log_y","xi_gauss",8,0,4),
        ("saved_profiles_theta_fine","saved_rows","cubic_log_y","theta_gauss",8,192,0),
        ("native_pchip_geometry","native_q","pchip_log_y","theta_gauss",8,192,0),
        ("native_remeasured_geometry","native_q","native_remeasured","theta_gauss",8,192,0),
    ]
    cases = []
    for args in specs:
        result = evaluate(*args)
        cases.append(result)
        print(json.dumps({"case":result["name"],**result["total"]}),flush=True)
    baseline = cases[1]
    sensitivity = {r["name"]:{k:r["total"][k]-baseline["total"][k]
                              for k in ["S_true","S_compressed","deltaS","tail_deltaS"]}
                   for r in cases if r is not baseline}
    report = dict(
        functional="S=-2/pi integral x*y/(x*x+y*y)^2 q dx dy",
        q_definition="q=rhoX in original tau20 coordinates; not peak-normalized inner coordinates",
        patch={"y":[.025,1024.],"xi":[-3.,3.],"beta":BETA,
               "map":"x=m+W*(sinh(beta*xi)+k*(cosh(beta*xi)-1))/beta"},
        reference={"definition":"unweighted mean of all measured rows with y<=1, including y=0",
                   "number_of_rows":int(ref_mask.sum()),"height_max":float(ys[ref_mask][-1]),
                   "xi_interpolation":"cubic spline on the saved 601 xi nodes"},
        interpretation=[
            "This is sensitivity of a finite-time geometric compression to a nonlocal functional; the far field has not converged.",
            "The correction is not an error lower bound for a limiting structure and does not exclude a single-profile limit; any need for tail modes requires cross-time checks.",
            "S_true and S_compressed are contributions from the same finite measured patch in each case.",
            "deltaS is integrated directly from q_true-q_compressed; no tail extrapolation is used.",
            "Left/right refer to xi intervals, not the sign of physical x; both physical signs are reported separately.",
            "The native field is the original saved-grid cubic interpolant, not an exact continuum solution.",
            "Remeasured geometry changes patch edges and compression slightly; it is a geometry-sensitivity check.",
            "No contribution from 0<=y<.025, y>1024, or |xi|>3 is included; this is not global Biot-Savart closure.",
            "Quadrature agreement does not establish mesh or domain convergence of the original trajectory."],
        baseline=baseline["name"],sensitivity_relative_to_baseline=sensitivity,cases=cases)
    (OUT/"strain_patch_checks.json").write_text(json.dumps(report,indent=2)+"\n")
    lines = ["# 已测补丁中的中心应变审核", "",
             "这是当前有限时刻几何压缩对非局部量的敏感性审核。远场尚未收敛，下述修正不是极限结构的误差下界，也不排除单剖面极限；候选是否需要独立尾模态须留待跨时间检验。", "",
             "只计算原始 τ=20 坐标下 $0.025\\le y\\le1024,-3\\le\\xi\\le3$ 的贡献。固定剖面是全部已测 $y\\le1$ 行（含壁面）的等权均值；它不是全域尾部模型。", "",
             "| 部分 | S_true | S_compressed | deltaS |", "|---|---:|---:|---:|"]
    for region,r in baseline["regions"].items():
        lines.append(f"| {region} | {r['S_true']:.10g} | {r['S_compressed']:.10g} | {r['deltaS']:.10g} |")
    r = baseline["total"]
    lines.append(f"| 总和 | {r['S_true']:.10g} | {r['S_compressed']:.10g} | {r['deltaS']:.10g} |")
    lines += ["", f"补丁总修正占其真实应变 {100*r['delta_over_true']:.5g}%；左右尾修正之和 {r['tail_deltaS']:.10g}，占总修正 {100*r['tail_delta_over_total_delta']:.5g}%。这里左/右以 ξ 定义，不能按 x 的符号替代。", "",
              "| 高度区间 | S_true | 总 deltaS | 两尾 deltaS |", "|---|---:|---:|---:|"]
    for height_bin in ["0.025..1","1..64","64..1024"]:
        parts = [b["height_bins"][height_bin] for b in baseline["regions"].values()]
        true_sum = sum(b["S_true"] for b in parts)
        delta_sum = sum(b["deltaS"] for b in parts)
        tail_sum = parts[0]["deltaS"]+parts[2]["deltaS"]
        lines.append(f"| {height_bin} | {true_sum:.10g} | {delta_sum:.10g} | {tail_sum:.10g} |")
    lines += ["",
              "稳定规则使用 $x=y\\tan\\theta,\\eta=\\log y$，将核及面积元准确变为 $-\\sin(2\\theta)d\\theta d\\eta/\\pi$，并在物理 x=0 处分段。另用原 ξ 网格每格四点 Gauss 求积作独立规则检查。高度积分在所有已测高度结点及 y=1,64 处分段。", "",
              "| 敏感性（相对主结果） | ΔS_true | ΔdeltaS | Δ尾修正 |", "|---|---:|---:|---:|"]
    for name,r in sensitivity.items():
        lines.append(f"| {name} | {r['S_true']:.4g} | {r['deltaS']:.4g} | {r['tail_deltaS']:.4g} |")
    lines += ["", "完整结果还按物理 x 正负和高度区间分解，见 strain_patch_checks.json。原场与保存剖面的插值差、几何插值差及求积差分别报告。", "",
              "没有计算 y<0.025、y>1024 或 |ξ|>3 的贡献，也没有验证原始轨道的网格/区域收敛。因此上述结果只度量已测形变对已测补丁应变的影响，不是完整速度闭合。"]
    (OUT/"STRAIN_PATCH.md").write_text("\n".join(lines)+"\n")


if __name__ == "__main__":
    main()

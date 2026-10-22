"""Test a renormalized geometric self-similarity closure for the tau=20 flow.

The constant drift and density source are removed exactly by

    z = x + v/lambda,        P = Q + v/mu.

Writing alpha=mu/lambda, p=1+alpha, t=log(r),

    P=r**alpha A(t,theta),   psi=r**p B(t,theta)

turns the far field into a fixed-point problem.  This script measures the
approach to that fixed point, tests the leading angular Poisson closure, and
checks the first nonlinear Mellin correction.  It is a read-only diagnostic;
it does not evolve the PDE or modify saved solutions.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/ipm_selfsimilar_mpl")

import numpy as np
from numpy.polynomial import chebyshev as cheb
from scipy.interpolate import RectBivariateSpline
from scipy.linalg import solve_banded

from streamfunction_sources import load_native


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/streamfunction_selfsimilar_closure_20260929"
NATIVE_RADII = 2.0 ** (np.arange(-16, 49) / 4)
THETA = (np.arange(401) + 0.5) * np.pi / 401


def relative(value: np.ndarray, reference: np.ndarray) -> float:
    return float(np.linalg.norm(value - reference) / np.linalg.norm(reference))


def normalized_spread(rows: np.ndarray) -> float:
    mean = rows.mean(axis=0)
    return relative(np.broadcast_to(mean, rows.shape), rows)


def chebyshev_fit(theta: np.ndarray, values: np.ndarray, degree: int):
    z = 2 * theta / np.pi - 1
    coefficients = cheb.chebfit(z, values, degree)

    def evaluate(order: int = 0) -> np.ndarray:
        derivative = cheb.chebder(coefficients, order)
        return cheb.chebval(z, derivative) * (2 / np.pi) ** order

    return coefficients, evaluate


def angular_poisson_from_density(theta: np.ndarray, A0: np.ndarray,
                                 alpha: float) -> np.ndarray:
    """Solve B0''+(1+alpha)^2 B0=-D_x[A0], B0(0)=B0(pi)=0."""
    _, afit = chebyshev_fit(theta, A0, 40)
    smooth = afit(0)
    derivative = afit(1)
    source = alpha * np.cos(theta) * smooth - np.sin(theta) * derivative
    h = theta[1] - theta[0]
    n = len(theta) - 2
    diagonal = np.full(n, -2 / h**2 + (1 + alpha) ** 2)
    off = np.full(n - 1, 1 / h**2)
    matrix = np.zeros((3, n))
    matrix[0, 1:] = off
    matrix[1] = diagonal
    matrix[2, :-1] = off
    B0 = np.zeros_like(theta)
    B0[1:-1] = solve_banded((1, 1), matrix, -source[1:-1])
    return B0


def sample(profile: dict, native_radii: np.ndarray, theta: np.ndarray,
           shifted: bool = True) -> dict:
    lam, mu, v = map(float, profile["rates"][:3])
    alpha = mu / lam
    radius = native_radii / profile["ell"]
    R, T = np.meshgrid(radius, theta, indexing="ij")
    center = -v / lam if shifted else 0.0
    x = center + R * np.cos(T)
    y = R * np.sin(T)
    splines = {key: RectBivariateSpline(profile["y"], profile["x"], profile[key])
               for key in ("Q", "q", "psi", "u1", "u2")}
    if x.min() < profile["x"][0] or x.max() > profile["x"][-1] or y.max() > profile["y"][-1]:
        raise ValueError("Polar samples leave the saved field")
    Q = splines["Q"].ev(y, x)
    q = splines["q"].ev(y, x)
    psi = splines["psi"].ev(y, x)
    u1 = splines["u1"].ev(y, x)
    u2 = splines["u2"].ev(y, x)
    Qy = splines["Q"].ev(y, x, dx=1)
    qx = splines["q"].ev(y, x, dy=1)
    qy = splines["q"].ev(y, x, dx=1)

    psi_r = np.cos(T) * u2 - np.sin(T) * u1
    psi_degree = np.sum(psi * R * psi_r, axis=1) / np.sum(psi**2, axis=1)
    psi_defect = np.linalg.norm(R * psi_r - psi_degree[:, None] * psi, axis=1) / np.linalg.norm(psi, axis=1)
    q_r = np.cos(T) * qx + np.sin(T) * qy
    q_degree = np.sum(q * R * q_r, axis=1) / np.sum(q**2, axis=1)
    q_defect = np.linalg.norm(R * q_r - q_degree[:, None] * q, axis=1) / np.linalg.norm(q, axis=1)

    shifted_density = Q + v / mu
    density_r = np.cos(T) * q + np.sin(T) * Qy
    dilation = lam * R * density_r - mu * shifted_density
    nonlinear = u1 * q + u2 * Qy
    full_transport = dilation + nonlinear
    transport_scale = np.linalg.norm(mu * shifted_density, axis=1)

    return dict(radius=radius, x=x, y=y, Q=Q, q=q, psi=psi,
                shiftedDensity=shifted_density, psiDegree=psi_degree,
                psiDefect=psi_defect, qDegree=q_degree, qDefect=q_defect,
                dilationRelative=np.linalg.norm(dilation, axis=1) / transport_scale,
                nonlinearRelative=np.linalg.norm(nonlinear, axis=1) / transport_scale,
                fullTransportRelative=np.linalg.norm(full_transport, axis=1) / transport_scale)


def outer_closure(profile: dict, observed: dict) -> tuple[dict, dict]:
    lam, mu, v = map(float, profile["rates"][:3])
    alpha = mu / lam
    power = 1 + alpha
    radius = observed["radius"]
    scaled_psi = observed["psi"] / radius[:, None] ** power
    scaled_q = observed["q"] * radius[:, None] ** (1 - alpha)

    mask = NATIVE_RADII >= 1024
    indices = np.flatnonzero(mask)
    train = indices[::2]
    holdout = indices[1::2]
    mean_scaled_psi = scaled_psi[train].mean(axis=0)
    coefficients, bfit = chebyshev_fit(THETA, mean_scaled_psi, 48)
    B0 = bfit(0)
    q0 = -(power**2 * B0 + bfit(2))
    angular_core = (THETA >= np.deg2rad(3)) & (THETA <= np.deg2rad(177))

    def closure_metrics(selected: np.ndarray) -> dict:
        return dict(
            psiFixedProfileRelativeL2=relative(
                np.broadcast_to(B0, scaled_psi[selected].shape), scaled_psi[selected]),
            qObservedCollapseRelativeL2=normalized_spread(scaled_q[selected][:, angular_core]),
            qFromLeadingPsiRelativeL2=relative(
                np.broadcast_to(q0[angular_core], scaled_q[selected][:, angular_core].shape),
                scaled_q[selected][:, angular_core]),
        )

    # Independently test the first nonlinear Mellin correction.  A0 is taken
    # from the largest observed density shells, then B0 is generated by the
    # angular Poisson equation rather than reusing the streamfunction fit.
    full_theta = np.linspace(0, np.pi, 401)
    full = sample(profile, NATIVE_RADII[NATIVE_RADII >= 256], full_theta, shifted=True)
    full_radius = full["radius"]
    full_scaled_density = full["shiftedDensity"] / full_radius[:, None] ** alpha
    far = NATIVE_RADII[NATIVE_RADII >= 256] >= 1024
    A0 = full_scaled_density[far].mean(axis=0)
    _, afit = chebyshev_fit(full_theta, A0, 40)
    A0 = afit(0)
    A0prime = afit(1)
    B0_from_A = angular_poisson_from_density(full_theta, A0, alpha)
    _, b_from_a_fit = chebyshev_fit(full_theta, B0_from_A, 50)
    B0prime = b_from_a_fit(1)
    A1 = (-alpha * B0prime * A0 + power * B0_from_A * A0prime) / (lam * (1 - alpha))
    beta2 = 2 * alpha - 1

    correction = {}
    for lower in (256, 1024):
        selected = NATIVE_RADII[NATIVE_RADII >= 256] >= lower
        actual = full["shiftedDensity"][selected]
        rr = full_radius[selected, None]
        leading = rr**alpha * A0
        corrected = leading + rr**beta2 * A1
        correction[str(lower)] = dict(
            nativeRadiusRange=[float(NATIVE_RADII[NATIVE_RADII >= 256][selected][0]),
                               float(NATIVE_RADII[NATIVE_RADII >= 256][selected][-1])],
            leadingDensityRelativeL2=relative(leading, actual),
            firstCorrectionDensityRelativeL2=relative(corrected, actual),
            correctionExponent=float(beta2),
        )

    result = dict(
        alpha=float(alpha), streamfunctionPower=float(power),
        dilationCenterInnerX=float(-v / lam),
        nativeRadiusRange=[1024.0, 4096.0],
        derivativeAngularRangeDegrees=[3.0, 177.0],
        trainingNativeRadii=NATIVE_RADII[train].tolist(),
        holdoutNativeRadii=NATIVE_RADII[holdout].tolist(),
        training=closure_metrics(train), holdout=closure_metrics(holdout),
        allOuter=closure_metrics(indices), firstMellinCorrection=correction,
        interpretation=("Leading streamfunction fixed profile is close on the measured outer shells; "
                        "derivative-level q closure is weaker, and the first formal nonlinear correction "
                        "does not improve this finite-time field."),
    )
    arrays = dict(B0=B0, q0=q0, scaledPsi=scaled_psi, scaledQ=scaled_q,
                  fullTheta=full_theta, densityA0=A0, streamfunctionB0FromDensity=B0_from_A,
                  densityA1=A1)
    return result, arrays


def angular_rank(observed: dict) -> dict:
    phi = observed["psi"] / np.sqrt(np.mean(observed["psi"]**2, axis=1))[:, None]
    singular = np.linalg.svd(phi, compute_uv=False)
    total = np.linalg.norm(singular)
    return dict(normalizedSingularValues=(singular[:8] / total).tolist(),
                rankRelativeErrors={str(rank): float(np.linalg.norm(singular[rank:]) / total)
                                    for rank in range(1, 7)})


def make_plot(profile: dict, core_observed: dict, shifted_observed: dict,
              closure: dict, arrays: dict) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lam, mu, _ = map(float, profile["rates"][:3])
    alpha = mu / lam
    power = 1 + alpha
    scaled_psi = arrays["scaledPsi"]
    scaled_q = arrays["scaledQ"]
    chosen = [int(np.argmin(abs(NATIVE_RADII - value))) for value in (256, 1024, 2048, 4096)]

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.facecolor": "white"})
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5), constrained_layout=True)
    ax = axes[0, 0]
    ax.semilogx(NATIVE_RADII, core_observed["psiDegree"], label=r"observed $p_\psi(R)$")
    ax.axhline(2, color="0.45", ls="--", label="regular core: 2")
    ax.axhline(power, color="#b24a3b", ls=":", lw=2,
               label=rf"outer fixed point: $1+\mu/\lambda={power:.4f}$")
    ax.set(xlabel="native radius R", ylabel="local radial degree",
           title="Inner and outer similarity exponents")
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    ax.loglog(NATIVE_RADII, core_observed["psiDefect"], label=r"$\psi$ angular drift")
    ax.loglog(NATIVE_RADII, core_observed["qDefect"], label=r"$q$ angular drift")
    ax.set(xlabel="native radius R", ylabel="relative homogeneity defect",
           title="The shape is slowly evolving, not exactly fixed")
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    for index in chosen:
        ax.plot(THETA * 180 / np.pi, scaled_psi[index], label=f"R={NATIVE_RADII[index]:g}")
    ax.plot(THETA * 180 / np.pi, arrays["B0"], "k--", lw=2, label="outer fixed profile")
    ax.set(xlabel="angle (degrees)", ylabel=rf"$\psi/r^{{{power:.3f}}}$",
           title="Outer streamfunction collapse")
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    for index in chosen:
        ax.plot(THETA * 180 / np.pi, scaled_q[index], lw=1,
                label=f"R={NATIVE_RADII[index]:g}")
    ax.plot(THETA * 180 / np.pi, arrays["q0"], "k--", lw=2,
            label=r"$-[B_0''+(1+\alpha)^2B_0]$")
    ax.set_xlim(3, 177)
    ax.set(xlabel="angle (degrees)", ylabel=rf"$r^{{1-{alpha:.3f}}}q$",
           title="Derivative-level Poisson closure")
    ax.legend(fontsize=8)
    for ax in axes.flat:
        ax.grid(alpha=0.16)
    fig.suptitle("Renormalized geometric self-similarity of the tau=20 streamfunction", fontsize=14)
    path = OUT / "selfsimilar_closure.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def write_report(profile: dict, closure: dict, rank: dict,
                 core_observed: dict, shifted_observed: dict) -> None:
    lam, mu, v = map(float, profile["rates"][:3])
    alpha = mu / lam
    power = 1 + alpha
    outer = closure["allOuter"]
    correction = closure["firstMellinCorrection"]
    selected = NATIVE_RADII >= 1024
    text = rf"""# 流函数几何自相似闭合试验

这次把图像中的自相似解释为一条**从核心二次固定点通向远场齐次固定点的重整化轨道**，而不是要求整个半平面只有一个幂律。

令

\[
z=x+\frac v\lambda,\qquad P=Q+\frac v\mu,
\qquad \alpha=\frac\mu\lambda={alpha:.9f},\qquad p_\infty=1+\alpha={power:.9f}.
\]

原稳态方程精确变成

\[
(u+\lambda(z,y))\cdot\nabla P=\mu P,\qquad P_z=-\Delta\psi.
\]

再写

\[
P=r^\alpha A(t,\theta),\qquad
\psi=r^{{1+\alpha}}B(t,\theta),\qquad t=\log r,
\]

得到精确的重整化方程

\[
-\big[(\partial_t+1+\alpha)^2+\partial_\theta^2\big]B
=\cos\theta(\partial_t+\alpha)A-\sin\theta A_\theta,
\]

\[
\lambda A_t+e^{{-(1-\alpha)t}}
\left[-B_\theta(A_t+\alpha A)+(B_t+(1+\alpha)B)A_\theta\right]=0.
\]

因此远场的固定点满足 `A_t=B_t=0`，自动给出流函数次数

\[
\boxed{{p_\infty=1+\mu/\lambda={power:.6f}}}.
\]

核心则由正则性强制为二次：

\[
\psi=r^2\left[(\lambda-\mu)\sin\theta\cos\theta-\tfrac12\sin^2\theta\right]+O(r^3).
\]

所以全局候选不是单幂，而是连接 `p=2` 与 `p={power:.6f}` 的异宿轨道。

## 数值检验

使用 τ=20 保存场。核心次数围绕实际壁面驻点测量，外层齐次塌缩围绕精确漂移中心 `x=-v/lambda={-v/lam:.6f}` 测量。实测流函数局部次数从最内层的 {core_observed['psiDegree'][0]:.6f} 降到最大保存半径的 {core_observed['psiDegree'][-1]:.6f}；其方向与理论远场次数一致，但尚未到达极限。这也说明内外固定点不应被强迫使用同一个有限半径中心。

在原生 `R=1024…4096` 上：

- `psi/r^p∞` 的固定角形误差为 {100*outer['psiFixedProfileRelativeL2']:.3f}%；
- `r^(1-alpha)q` 自身的角形漂移为 {100*outer['qObservedCollapseRelativeL2']:.3f}%；
- 在避开墙端数值二阶导数的 `3°…177°` 角带内，用固定流函数角形经 Poisson 算子生成 q，误差为 {100*outer['qFromLeadingPsiRelativeL2']:.3f}%。

整个半径区间的流函数角形矩阵并非严格秩一，但低秩很强：rank 1/2/3/4 的值误差依次为 {100*rank['rankRelativeErrors']['1']:.3f}%、{100*rank['rankRelativeErrors']['2']:.3f}%、{100*rank['rankRelativeErrors']['3']:.3f}%、{100*rank['rankRelativeErrors']['4']:.3f}%。这支持“少数角模态沿 log r 演化”，不支持“一个固定角形覆盖全部半径”。

## 第一项非线性修正

远场展开

\[
P\sim r^\alpha A_0+r^{{2\alpha-1}}A_1+\cdots
\]

中的第一项修正由

\[
A_1=\frac{{-\alpha B_0'A_0+(1+\alpha)B_0A_0'}}{{\lambda(1-\alpha)}}
\]

唯一生成。本次直接用最大外层提取 `A0` 后，加入这项并没有改善有限时刻数据：`R=1024…4096` 的密度相对误差从 {100*correction['1024']['leadingDensityRelativeL2']:.3f}% 变为 {100*correction['1024']['firstCorrectionDensityRelativeL2']:.3f}%；`R=256…4096` 从 {100*correction['256']['leadingDensityRelativeL2']:.3f}% 变为 {100*correction['256']['firstCorrectionDensityRelativeL2']:.3f}%。因此目前数据只支持接近远场固定点，尚不支持把前两项渐近级数当作已经闭合的 τ=20 全局解。

更直接的障碍是：在最大半径处，纯伸缩输运缺陷 `||lambda r P_r-mu P||/||mu P||` 仍为 {100*shifted_observed['dilationRelative'][-1]:.2f}%，而自诱导非线性项只有 {100*shifted_observed['nonlinearRelative'][-1]:.2f}%。这说明 τ=20 外层的流函数几何已经接近齐次形状，但密度输运尚未到稳态固定点。

## 结论

可保留的闭合结构是

\[
\boxed{{\psi=r^{{1+\mu/\lambda}}B(\log r,\theta)}}
\]

并让少数角模态系数随 `log r` 形成从核心二次态到远场固定态的轨道。它同时解释原点有界和外层约 1.38 的表观次数。当前 τ=20 数据在流函数值层面明显支持该结构，但二阶 Poisson 与密度输运层面尚未闭合；不能把它宣布为已求得的稳态解。
"""
    (OUT / "REPORT_ZH.md").write_text(text)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    profile = load_native(20)
    core_observed = sample(profile, NATIVE_RADII, THETA, shifted=False)
    shifted_observed = sample(profile, NATIVE_RADII, THETA, shifted=True)
    closure, arrays = outer_closure(profile, shifted_observed)
    rank = angular_rank(core_observed)
    result = dict(
        source=profile["source"], tau=profile["tau"],
        rates=dict(lambda_=float(profile["rates"][0]), mu=float(profile["rates"][1]),
                   v=float(profile["rates"][2])),
        exactReduction=dict(
            shiftedCoordinates="z=x+v/lambda, P=Q+v/mu",
            alpha="mu/lambda", streamfunctionPower="1+mu/lambda",
            poisson="-[(d_t+1+alpha)^2+d_theta^2]B=cos(theta)(d_t+alpha)A-sin(theta)A_theta",
            transport="lambda A_t+exp(-(1-alpha)t)[-B_theta(A_t+alpha A)+(B_t+(1+alpha)B)A_theta]=0"),
        outerClosure=closure, angularRank=rank,
        radialDiagnostics=dict(nativeRadii=NATIVE_RADII.tolist(),
            coreCenter=dict(psiDegree=core_observed["psiDegree"].tolist(),
                            psiDefect=core_observed["psiDefect"].tolist(),
                            qDegree=core_observed["qDegree"].tolist(),
                            qDefect=core_observed["qDefect"].tolist()),
            dilationCenter=dict(dilationRelative=shifted_observed["dilationRelative"].tolist(),
                                nonlinearRelative=shifted_observed["nonlinearRelative"].tolist(),
                                fullTransportRelative=shifted_observed["fullTransportRelative"].tolist())),
        limitations=["tau=20 is a finite-time field, not a certified steady state",
                     "R=1024..4096 is a measured outer band, not an infinite-radius limit",
                     "value-level angular low rank does not certify derivative-level closure",
                     "the pure-y part of shifted density is invisible to P_x and must be fixed by transport/matching"])
    (OUT / "closure.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    np.savez_compressed(OUT / "closure_data.npz", nativeRadii=NATIVE_RADII,
                        theta=THETA, psiDegree=core_observed["psiDegree"],
                        psiDefect=core_observed["psiDefect"], qDegree=core_observed["qDegree"],
                        qDefect=core_observed["qDefect"], **arrays)
    make_plot(profile, core_observed, shifted_observed, closure, arrays)
    write_report(profile, closure, rank, core_observed, shifted_observed)
    print(json.dumps(dict(alpha=closure["alpha"], power=closure["streamfunctionPower"],
        outer=closure["allOuter"], correction=closure["firstMellinCorrection"],
        rank=rank["rankRelativeErrors"],
        lastTransport=dict(dilation=float(shifted_observed["dilationRelative"][-1]),
                           nonlinear=float(shifted_observed["nonlinearRelative"][-1]),
                           full=float(shifted_observed["fullTransportRelative"][-1])),
        output=str(OUT)), indent=2))


if __name__ == "__main__":
    main()

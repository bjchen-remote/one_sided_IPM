#!/usr/bin/env python3
"""Draw the actual transport and spectral meshes exported by export_mesh.m.

Requires NumPy, SciPy, Matplotlib and the bundled mesh rendering settings.
Every exported mesh line is
retained; a near-field figure clips lines to its native-coordinate window
rather than thinning the mesh. Captions and SHA-256 provenance are written
alongside PNG, SVG and PDF figures.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, Patch, Rectangle
from matplotlib import patheffects
from scipy.io import loadmat


STYLE_PROFILE = Path(__file__).resolve().parent/"styles"/"mesh.json"


def load_profile(path: Path):
    """Load the portable profile; its palette order defines mesh semantics."""
    profile = json.loads(path.read_text())
    if profile.get("name") != "mesh" or len(profile.get("colors", [])) < 3:
        raise ValueError("Expected the bundled mesh settings with three or more colors.")
    return profile


PROFILE = load_profile(STYLE_PROFILE)
STYLE = {**PROFILE["rc_params"], "svg.hashsalt": "ipm_cr_mesh_v1"}
# The first three categorical colors retain their mesh roles across panels.
# Neutral display colors come from named rc tokens; no plotted grid is added.
COLORS = dict(zip(("core", "polar", "ghost"), PROFILE["colors"][:3]))
COLORS.update(overlap=STYLE["grid.color"], panel=STYLE["ytick.color"],
              ink=STYLE["text.color"], paper=STYLE["axes.facecolor"])


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def vector(value) -> np.ndarray:
    return np.asarray(value, dtype=float).reshape(-1)


def physical_axes(ax, xlim, ylim) -> None:
    ax.set(xlim=xlim, ylim=ylim, xlabel="$X$", ylabel="$Y$")
    ax.set_aspect("equal", adjustable="box")


def visible_nodes(X, Y, ax, expand=False):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    visible = (X >= xmin) & (X <= xmax) & (Y >= ymin) & (Y <= ymax)
    if expand:
        # Include adjacent nodes so a clipped line reaches the frame exactly.
        grown = visible.copy()
        grown[1:] |= visible[:-1]
        grown[:-1] |= visible[1:]
        grown[:, 1:] |= visible[:, :-1]
        grown[:, :-1] |= visible[:, 1:]
        visible = grown
    return visible


def grid_lines(ax, X, Y, mask=None, *, color, width=0.4, alpha=0.7,
               linestyle="solid", zorder=2, clip=True, touch=False):
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    if X.shape != Y.shape or X.ndim != 2:
        raise ValueError("Grid coordinate arrays must be matching matrices.")
    mask = np.ones_like(X, dtype=bool) if mask is None else np.asarray(mask, bool)
    if mask.shape != X.shape:
        raise ValueError("The active/ghost mask must match its coordinate grid.")
    if clip:
        mask = mask & visible_nodes(X, Y, ax, expand=True)
    segments = []
    # Every row and column is included; no stride or mesh-line subsampling.
    for xx, yy, mm in ((X, Y, mask), (X.T, Y.T, mask.T)):
        for xrow, yrow, mrow in zip(xx, yy, mm):
            # Active edges have two active endpoints. A ghost edge has at
            # least one ghost endpoint, including the edge to the nearest
            # active node; that node itself remains active in the marker plot.
            edge_mask = (mrow[:-1] | mrow[1:]) if touch else (mrow[:-1] & mrow[1:])
            ids = np.flatnonzero(edge_mask)
            for run in np.split(ids, np.flatnonzero(np.diff(ids) > 1) + 1):
                if len(run):
                    run = np.append(run, run[-1]+1)
                    segments.append(np.column_stack((xrow[run], yrow[run])))
    collection = LineCollection(segments, colors=color, linewidths=width,
                                alpha=alpha, linestyles=linestyle, zorder=zorder)
    ax.add_collection(collection)
    return collection


def nodes(ax, grid, mask, *, color, size=7, hollow=False, zorder=5):
    X, Y = np.asarray(grid["X"]), np.asarray(grid["Y"])
    keep = np.asarray(mask, bool) & visible_nodes(X, Y, ax)
    return ax.scatter(X[keep], Y[keep], s=size, marker="o", linewidths=0.55,
                      facecolors=COLORS["paper"] if hollow else color,
                      edgecolors=color, zorder=zorder)


def arc(ax, radius, *, color, width=1.0, linestyle="solid", alpha=1, zorder=4):
    t = np.linspace(0, np.pi, 1001)
    return ax.plot(radius*np.cos(t), radius*np.sin(t), color=color,
                   lw=width, ls=linestyle, alpha=alpha, zorder=zorder)[0]


def box_line(ax, box, *, color, width=1.0, linestyle="solid", zorder=4):
    left, right, low, high = vector(box)
    return ax.plot([left, left, right, right], [low, high, high, low],
                   color=color, lw=width, ls=linestyle, zorder=zorder)[0]


def overlap_fill(ax, box, radius, *, alpha=0.85):
    left, right, low, high = vector(box)
    x = np.linspace(left, right, 1501)
    bottom = np.maximum(low, np.sqrt(np.maximum(0, radius*radius-x*x)))
    ax.fill_between(x, bottom, high, where=bottom <= high,
                    color=COLORS["overlap"], alpha=alpha, linewidth=0, zorder=0.4)


def transport_lines(ax, m, *, markers=False, faint=False):
    core, polar = m["transport"]["core"], m["transport"]["polar"]
    overlap_fill(ax, m["bounds"]["coreActiveBox"], polar["innerActiveRadius"],
                 alpha=0.60 if faint else 0.85)
    for grid, category in ((polar, "polar"), (core, "core")):
        grid_lines(ax, grid["X"], grid["Y"], grid["active"],
                   color=COLORS[category],
                   width=0.21 if faint else (0.38 if markers else 0.27),
                   alpha=0.43 if faint else (0.70 if markers else 0.56))
    for grid in (polar, core):
        grid_lines(ax, grid["X"], grid["Y"], grid["ghost"],
                   color=COLORS["ghost"], width=0.48 if markers else 0.37,
                   linestyle=(0, (2.2, 2.2)), alpha=0.8, zorder=3, touch=True)
    if markers:
        for grid, category in ((polar, "polar"), (core, "core")):
            nodes(ax, grid, grid["active"], color=COLORS[category], size=7)
        for grid in (polar, core):
            nodes(ax, grid, grid["ghost"], color=COLORS["ghost"], size=12, hollow=True)
    box_line(ax, m["bounds"]["coreActiveBox"], color=COLORS["core"], width=1.05)
    arc(ax, polar["innerActiveRadius"], color=COLORS["polar"], width=1.05)
    ax.axhline(0, color=COLORS["ink"], lw=0.75, zorder=7)


def transport_legend(*, markers=False, overlap=True):
    handles = [
        Line2D([], [], color=COLORS["core"], lw=1, label="Core",
               marker="o" if markers else None, markersize=3),
        Line2D([], [], color=COLORS["polar"], lw=1, label="Annulus",
               marker="o" if markers else None, markersize=3),
        Line2D([], [], color=COLORS["ghost"], lw=1, ls="--", label="Ghost",
               marker="o" if markers else None, markerfacecolor=COLORS["paper"], markersize=4),
    ]
    if overlap:
        handles.append(Patch(facecolor=COLORS["overlap"], alpha=0.85,
                             edgecolor="none", label="Overlap"))
    return handles


def collar_view(m):
    """The same real Cartesian ghost window used in the connection figure."""
    p = m["transport"]["polar"]
    h = float(np.diff(vector(m["transport"]["core"]["y"]))[0])
    span = max(0.30, float(p["innerActiveRadius"]-p["innerHaloRadius"])+6*h)
    high = vector(m["bounds"]["coreBox"])[3]
    return (0.34, 0.34+span), (high-span+0.035, high+0.035)


def near_figure(m):
    # Unequal panels allocate space to geometry rather than equal-size cards.
    # Fixed axes are the sole layout mechanism for this linked overview/crop.
    fig = plt.figure(figsize=(8.6, 3.9), layout="none")
    ax = fig.add_axes((0.06, 0.19, 0.65, 0.715))
    detail = fig.add_axes((0.78, 0.51, 0.18, 0.397))
    extent = max(abs(vector(m["bounds"]["coreBox"])))*1.7
    physical_axes(ax, (-extent, extent), (0, extent))
    transport_lines(ax, m, faint=True)
    box_line(ax, m["bounds"]["coreBox"], color=COLORS["ghost"],
             linestyle=(0, (2.2, 2.2)), width=0.55)
    ax.set_xticks((-2, -1, 0, 1, 2))
    ax.set_yticks((0, 1, 2))
    xlim, ylim = collar_view(m)
    x0, y0 = xlim[0], ylim[0]
    span = xlim[1]-x0
    physical_axes(detail, xlim, ylim)
    transport_lines(detail, m, markers=True)
    detail.set(xlabel="", ylabel="")
    detail.locator_params(axis="both", nbins=3)
    detail.tick_params(labelsize=7.5, length=2, pad=2)
    for spine in detail.spines.values():
        spine.set_visible(True)
        spine.set_color(COLORS["ghost"])
        spine.set_linewidth(0.6)
    ax.add_patch(Rectangle((x0, y0), span, span, fill=False,
                          edgecolor=COLORS["ink"], lw=0.8, zorder=9))
    fig.add_artist(ConnectionPatch(
        xyA=(x0+span, y0+span), coordsA=ax.transData,
        xyB=(x0, y0+span), coordsB=detail.transData,
        color=COLORS["ghost"], lw=0.65, zorder=1, clip_on=False))
    fig.text(0.06, 0.945, "a", fontsize=10, weight="bold")
    fig.text(0.085, 0.945, "Near field", fontsize=9.5)
    fig.text(0.78, 0.945, "b", fontsize=10, weight="bold")
    fig.text(0.805, 0.945, "Ghost collar", fontsize=9.5)
    fig.legend(handles=transport_legend(markers=True), loc="upper left",
               bbox_to_anchor=(0.777, 0.435), ncol=1, fontsize=8.5,
               handlelength=1.7, labelspacing=0.6, borderaxespad=0)
    return fig


def join_figure(m):
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 4.4), layout="constrained")
    # Leave a real top margin: some Matplotlib releases undercount titles
    # when fixed-aspect axes share an outside figure legend.
    fig.get_layout_engine().set(rect=(0, 0, 1, 0.94))
    p = m["transport"]["polar"]
    h = float(np.diff(vector(m["transport"]["core"]["y"]))[0])
    r0, ra = float(p["innerHaloRadius"]), float(p["innerActiveRadius"])
    span = max(0.30, (ra-r0)+6*h)
    xmid = (r0+ra)/2
    physical_axes(axes[0], (xmid-span/2, xmid+span/2), (0, span))
    physical_axes(axes[1], *collar_view(m))
    for ax in axes:
        transport_lines(ax, m, markers=True)
        ax.locator_params(axis="both", nbins=5)
    axes[0].set_title("(a) Annular collar", loc="left", pad=9)
    axes[1].set_title("(b) Cartesian collar", loc="left", pad=9)
    fig.legend(handles=transport_legend(markers=True), loc="outside lower center",
               ncol=4, handlelength=1.6, columnspacing=1.5)
    return fig


def warp(X, Y):
    R = np.hypot(X, Y)
    factor = np.divide(np.log10(1+R), R, out=np.full_like(R, 1/np.log(10)), where=R != 0)
    return X*factor, Y*factor


def global_figure(m):
    fig = plt.figure(figsize=(8.6, 4.8), layout="none")
    ax = fig.add_axes((0.065, 0.15, 0.915, 0.81))
    rmax = float(m["bounds"]["outer"])
    limit = np.log10(1+rmax)
    ax.set(xlim=(-limit*1.025, limit*1.025), ylim=(0, limit*1.04),
           xlabel=r"$\xi$", ylabel=r"$\upsilon$")
    ax.set_aspect("equal", adjustable="box")
    for category in ("polar", "core"):
        grid = m["transport"][category]
        X, Y = warp(np.asarray(grid["X"]), np.asarray(grid["Y"]))
        grid_lines(ax, X, Y, grid["active"], color=COLORS[category],
                   width=0.18 if category == "polar" else 0.30,
                   alpha=0.40 if category == "polar" else 0.75, clip=False)
        grid_lines(ax, X, Y, grid["ghost"], color=COLORS["ghost"],
                   width=0.30, alpha=0.75, linestyle="dashed", clip=False, touch=True)
    # Outline the actual exported outer row, not an invented coarser ring.
    polar = m["transport"]["polar"]
    X, Y = warp(np.asarray(polar["X"])[-1], np.asarray(polar["Y"])[-1])
    ax.plot(X, Y, color=COLORS["polar"], lw=0.85, alpha=0.9, zorder=4)
    # Three short leaders locate exact native radii; they are annotations,
    # not replacement grid rows. Text halos avoid large opaque label boxes.
    for exponent, angle, offset in ((2, 135, (-27, 13)),
                                     (6, 125, (-28, 16)),
                                     (10, 102, (-24, 11))):
        radius = np.log10(1+10.0**exponent)
        theta = np.deg2rad(angle)
        label = ax.annotate(
            rf"$r=10^{{{exponent}}}$", xy=(radius*np.cos(theta), radius*np.sin(theta)),
            xytext=offset, textcoords="offset points", fontsize=8.5,
            ha="right", va="bottom", color=COLORS["ink"], zorder=7,
            arrowprops={"arrowstyle": "-", "lw": 0.6, "color": COLORS["panel"]})
        label.set_path_effects([patheffects.withStroke(linewidth=2.5, foreground=COLORS["paper"])])
    ax.text(0.99, 0.99, "$s=\\log_{10}(1+r)$\n$(\\xi,\\upsilon)=s(\\cos\\theta,\\sin\\theta)$",
            transform=ax.transAxes, ha="right", va="top", fontsize=8.5,
            linespacing=1.55, color=COLORS["panel"])
    ax.set_xticks((-10, -5, 0, 5, 10))
    ax.set_yticks(np.arange(0, 11, 2))
    fig.legend(handles=transport_legend(overlap=False), loc="lower center",
               bbox_to_anchor=(0.53, 0.015), ncol=3, handlelength=1.8,
               columnspacing=2.0, fontsize=8.5)
    return fig


def poisson_figure(m):
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.5), layout="constrained",
                             gridspec_kw={"width_ratios": (1.3, 1)})
    core, polar = m["poisson"]["core"], m["poisson"]["polar"]
    a = max(abs(vector(m["bounds"]["spectralCoreBox"])))
    physical_axes(axes[0], (-1.65*a, 1.65*a), (0, 1.65*a))
    grid_lines(axes[0], core["X"], core["Y"], color=COLORS["core"], width=0.26, alpha=0.62)
    grid_lines(axes[0], polar["X"], polar["Y"], color=COLORS["polar"], width=0.23, alpha=0.50)
    box_line(axes[0], m["bounds"]["spectralCoreBox"], color=COLORS["core"])
    arc(axes[0], m["bounds"]["spectralInner"], color=COLORS["polar"])
    for edge in np.exp(vector(polar["logEdges"]))[1:-1]:
        if edge <= np.hypot(1.65*a, 1.65*a):
            arc(axes[0], edge, color=COLORS["panel"], width=0.70, alpha=0.85,
                linestyle=(0, (4, 2)))
    axes[0].axhline(0, color=COLORS["ink"], lw=0.75, zorder=7)
    axes[0].set_xticks((-1.5, 0, 1.5))
    axes[0].set_yticks((0, 0.75, 1.5))
    theta, t = vector(polar["theta"]), vector(polar["t"])
    T, TH = np.meshgrid(t, theta, indexing="ij")
    axes[1].set(xlim=(0, np.pi), ylim=(min(t), max(t)), xlabel=r"$\theta$", ylabel=r"$t=\log r$")
    # Equal panel heights, while only the native-coordinate panel is isometric.
    axes[1].set_box_aspect(0.65)
    grid_lines(axes[1], TH, T, color=COLORS["polar"], width=0.20, alpha=0.50, clip=False)
    for edge in vector(polar["logEdges"]):
        axes[1].axhline(edge, color=COLORS["panel"], lw=0.6, ls=(0, (4, 2)), alpha=0.8, zorder=3)
    axes[1].set_xticks((0, np.pi/2, np.pi), ("$0$", "$\\pi/2$", "$\\pi$"))
    axes[1].set_yticks((0, 5, 10, 15, 20))
    for ax, label in zip(axes, ("(a) Native near field", "(b) Full log-polar annulus")):
        ax.set_title(label, loc="left", pad=9)
    handles = [Line2D([], [], color=COLORS["core"], lw=1, label="Chebyshev core"),
               Line2D([], [], color=COLORS["polar"], lw=1, label="Spectral annulus"),
               Line2D([], [], color=COLORS["panel"], lw=1, ls=(0, (4, 2)), label="Panel boundary")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3,
               handlelength=1.8, columnspacing=2)
    return fig


def captions(m):
    count, bounds = m["counts"], m["bounds"]
    extent = max(abs(vector(bounds["coreBox"])))*1.7
    ra = float(bounds["polarInnerActive"])
    p = m["poisson"]["polar"]
    tail_density = float(m["radial"]["actualTailDensityFactor"])
    return {
        "mesh-near": {
            "caption": f"Transport mesh in the solver's native (X,Y) coordinates with equal axis scaling in both panels. "
            f"The Cartesian core contains {count['coreActive']:,} active nodes; the annulus "
            f"contains {count['polarActive']:,}. Blue and coral lines mark their active "
            f"interiors. Both evolve in the shaded overlap. Gray dashed lines are the "
            f"three-layer interpolation collars. The active annulus starts at r={ra:.6f}; "
            "the dark line is the wall Y=0. Panel (a) retains the original overview window; "
            "the outlined ROI is magnified in (b), using the same Cartesian collar window as mesh-join (b). "
            "Filled colored points are active nodes and hollow gray points are ghosts. All mesh lines "
            "intersecting either window are shown; no donor links are inferred. The outer domain extends to r=10^10.",
            "caption_zh": "含时输运网格在求解器原生 (X,Y) 坐标中的近场图，两窗均横纵等比例。蓝色为张量核心的活动节点网格，"
            "珊瑚色为极坐标环带的活动节点网格；浅灰蓝底是两套活动自由度独立推进的重叠区。"
            "灰虚线为三层插值 ghost，黑线是壁面 Y=0。(a) 保留原概览窗口，框选区域在 (b) 等比例放大；"
            "彩色实心点为活动点，灰色空心点为 ghost。放大窗与 mesh-join (b) 相同；未画 donor 连线，未对窗口内网格线抽样。",
            "coordinates": "native, unwarped X,Y; equal aspect in both windows; no time-dependent coordinate conversion",
            "overview_window_xy": [-extent, extent, 0, extent],
            "detail_window_xy": [*collar_view(m)[0], *collar_view(m)[1]],
        },
        "mesh-join": {
            "caption": "Two native-coordinate magnifications of the real interpolation "
            "collars. (a) The three inward annular ghost rows receive values from active "
            "Cartesian donors. (b) The three top Cartesian ghost rows receive values from "
            "active annular donors. Filled blue/coral circles are active nodes; hollow "
            "gray circles are interpolation ghosts. Both axes use equal native-coordinate scaling "
            "within each panel. All visible nodes and lines are retained; the two panels "
            "show different windows. The 7-by-7 tensor interpolation stencils are not drawn.",
            "caption_zh": "原生 (X,Y) 坐标中真实接缝的两个等比例局部窗口。(a) 环带内侧三层 ghost 从活动张量核心取值；"
            "(b) 张量核心顶部三层 ghost 从活动环带取值。蓝/珊瑚实心点是活动节点，灰空心点是插值 ghost。"
            "显示窗口内全部节点和网格线；两窗覆盖不同区域，7×7 张量插值模板未画出。",
            "coordinates": "native X,Y; equal aspect within each of two windows",
        },
        "mesh-global": {
            "caption": f"The full transport mesh under the display-only radial transform "
            "r=sqrt(X^2+Y^2), s=log10(1+r), (xi,upsilon)=s(cos(theta),sin(theta)). "
            f"This is a mapped view. All {count['radialAnnulus']} retained annular "
            f"radial rows and {count['angularTransport']} angular columns are drawn, reaching "
            "r=10^10. The small blue patch is the same Cartesian core under the same "
            "transform. Three fine leaders identify exact native-coordinate radii 10^2, 10^6 and 10^10; "
            "these labels are geometric references, not replacement mesh rows. The stronger outer contour "
            "is the actual last exported radial row. Redistribution preserves "
            f"r<=4 and uses a configured tail density factor of 2 (realized {tail_density:.6f}), "
            "with full refinement from r=20; no density "
            "field or time trajectory was computed to make this figure.",
            "caption_zh": "全域输运网格使用仅用于显示的径向变换 s=log10(1+r)，"
            "(ξ,υ)=s(cosθ,sinθ)，其中 r=√(X²+Y²)，不代表原生坐标等比例。完整保留环带所有径向行和角向列，外边界 r=10^10。"
            f"蓝色小片为同一变换下的张量核心。径向加密保留 r≤4，r≥20 的尾部加密配置为 2 倍，"
            f"实际密度比 {tail_density:.6f}；三条细引线标出原生半径 10²、10⁶、10¹⁰ 的位置，"
            "不替代实际网格行。较强外轮廓取自实际最外层径向网格行；本图未推进物理场。",
            "coordinates": "xi=log10(1+r)*cos(theta), upsilon=log10(1+r)*sin(theta); display transform only",
        },
        "mesh-poisson": {
            "caption": "The independent spectral Poisson discretization. (a) Native "
            "near field, with equal axis scaling: Chebyshev tensor core in blue and "
            "spectral annulus in coral. (b) The complete annulus in (theta,t=log r) "
            f"coordinates, showing all {int(p['nPanels'])} radial panels of degree "
            f"{int(p['degree'])} and angular degree {int(p['angularDegree'])}. Graphite dashed lines "
            "are actual panel boundaries; adjoining panels share their endpoint node, "
            "where derivative continuity is imposed. Every exported grid line is "
            "retained. These are distinct from the finite-difference transport grids.",
            "caption_zh": "独立的谱 Poisson 网格。(a) 等比例原生 (X,Y) 近场：蓝色 Chebyshev 张量核心，"
            "珊瑚色谱环带；(b) (θ,t=log r) 坐标下的完整环带。石墨色虚线为实际径向分片边界，"
            "相邻片共享端点并施加导数连续条件。未抽样网格线；此网格与输运差分网格相互独立。",
            "coordinates": "panel a: native X,Y, equal aspect; panel b: theta,t=natural log(r)",
        },
    }


def save(fig, stem: Path):
    files = []
    for suffix in ("png", "svg", "pdf"):
        path = stem.with_suffix("."+suffix)
        metadata = {"Creator": "ipm_cr tools/plot_mesh.py"}
        if suffix == "pdf":
            metadata.update({"CreationDate": None, "ModDate": None})
        if suffix == "svg":
            metadata["Date"] = None
        # savefig restores a None layout engine through rcParams. Keep manual
        # figures manual across PNG/SVG/PDF; existing constrained engines stay.
        with plt.rc_context({"figure.constrained_layout.use": False,
                             "figure.autolayout": False}):
            fig.savefig(path, dpi=250, metadata=metadata)
        files.append({"file": path.name, "sha256": digest(path), "bytes": path.stat().st_size})
    plt.close(fig)
    return files


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="mesh.mat from examples/export_mesh.m")
    parser.add_argument("--output", type=Path, required=True, help="Directory for figures and provenance")
    args = parser.parse_args()
    if not args.input.is_file():
        parser.error(f"Missing mesh export: {args.input}")
    exported = loadmat(args.input, simplify_cells=True)
    if "m" not in exported or exported["m"].get("format") != "ipm_cr_mesh_v1":
        parser.error("Expected m.format=ipm_cr_mesh_v1 from export_mesh.m")
    m = exported["m"]
    input_hash = digest(args.input)
    source_info = args.input.with_suffix(".json")
    source = json.loads(source_info.read_text()) if source_info.exists() else None
    if source is not None and source.get("sha256") != input_hash:
        raise ValueError("MAT export SHA-256 does not match its JSON provenance.")
    args.output.mkdir(parents=True, exist_ok=True)
    records = captions(m)
    constructors = {"mesh-near": near_figure, "mesh-join": join_figure,
                    "mesh-global": global_figure, "mesh-poisson": poisson_figure}
    with plt.rc_context(STYLE):
        for name, constructor in constructors.items():
            fig = constructor(m)
            records[name]["size_inches"] = list(fig.get_size_inches())
            records[name]["files"] = save(fig, args.output/name)
    provenance = {
        "format": "ipm_cr_mesh_figures_v1",
        "input": args.input.name,
        "input_sha256": input_hash,
        "script": "tools/plot_mesh.py",
        "script_sha256": digest(Path(__file__)),
        "style": PROFILE["name"],
        "style_profile": "tools/styles/mesh.json",
        "style_profile_sha256": digest(STYLE_PROFILE),
        "style_details": {
            "palette_roles": COLORS,
            "axes": "white background, thin left/bottom axes, sparse outward ticks",
            "typography": "portable sans-serif; bundled profile controls font sizes",
            "distinction": "solid active nodes; hollow ghost nodes; dashed collars and panel joins",
            "grid_lines": "all true mesh lines retained, thin strokes for dense meshes",
            "composition": "near: original overview with a linked equal-aspect real collar crop; global: full mapped mesh, three radius leaders and the actual outer contour",
        },
        "line_subsampling": False,
        "no_time_integration": True,
        "no_physical_field_data": True,
        "matplotlib": matplotlib.__version__,
        "source": source,
        "figures": records,
    }
    (args.output/"figures.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False)+"\n")
    caption_text = "\n\n".join(f"{name}\n{record['caption']}\n{record['caption_zh']}" for name, record in records.items())
    (args.output/"captions.txt").write_text(caption_text+"\n")
    print(json.dumps({"output": str(args.output), "figures": list(records),
                      "input_sha256": input_hash}, ensure_ascii=False))


if __name__ == "__main__":
    main()

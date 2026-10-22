"""Read-only, consistently centred streamfunction observations.

Native checkpoints store rho, not flow.  The tau10/15 archive contains one
static restoreCheckpoint evaluation per saved state, with no time evolution.
Native tau20 uses the original run_analysis.mat export.  Conditional roots
remain frozen-exterior finite-box observations, not full half-plane solutions.

All returned fields use (ny,nx) layout and inner coordinates.  The sign is
u=(-psi_y,psi_x): centring u1 by subtracting U0 ADDS U0*y to psi.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.interpolate import CubicSpline

from fixed_point import fd_matrix

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "runs/steady_20260928_tau20"
ARCHIVE = ROOT / "runs/streamfunction_geometry_20260928/source_frames.npz"
CONDITIONAL_NAMES = ("L4_129x65", "L4_257x129", "L6_193x97", "L8_257x129_continued")


def load_native(tau: int) -> dict:
    """Return centred inner psi/u and the saved Q/q for tau=10,15,20.

    Keys: x,y,Q,q,psi,u1,u2,psi_saved,a,P,ell,rates,metadata,source,gauge,
    name,kind,U0. psi_saved is the uncentred streamfunction in inner units.
    U0 is the removed native horizontal speed (before dividing by P*ell).
    """
    tau = int(tau)
    if tau not in (10, 15, 20):
        raise ValueError("Only the observed tau10/15/20 frames are supported")
    frame_path = ROOT / f"runs/geometry_loop_20260928/frame_tau{tau}.mat"
    frame = loadmat(frame_path, simplify_cells=True)
    a, P, ell = (float(frame[k]) for k in ("a", "P", "ell"))
    if tau == 20:
        raw_path = SOURCE / "run_analysis.mat"
        raw = loadmat(raw_path, simplify_cells=True)
        meta = dict(raw["metadata"])
        psi_abs, ua, va = (np.asarray(raw[k], float) for k in ("psi", "u1", "u2"))
        U0 = float(meta["uAtPeak"])
        rates = np.asarray(meta["innerRates"], float)
        y_native = np.asarray(raw["yFull"], float)
        metadata = dict(meta, restoration="original_saved_analysis_export")
    else:
        if not ARCHIVE.exists():
            raise FileNotFoundError(f"Static checkpoint flow archive is not yet available: {ARCHIVE}")
        raw_path = ARCHIVE
        with np.load(ARCHIVE, allow_pickle=False) as raw:
            prefix = f"tau{tau}_"
            psi_abs = raw[prefix + "psi_abs"]
            ua, va = raw[prefix + "u1_abs"], raw[prefix + "u2_abs"]
            U0 = float(raw[prefix + "uAtPeak"])
            y_native = raw[prefix + "yFull"]
            rates = raw[prefix + "rates"]
            metadata = json.loads(str(raw[prefix + "metadata_json"]))
    psi_saved = psi_abs / (P * ell**2)
    psi = (psi_abs + U0 * y_native[:, None]) / (P * ell**2)
    x, y = np.asarray(frame["x"], float), np.asarray(frame["y"], float)
    assert psi.shape == (len(y), len(x)) == np.asarray(frame["q"]).shape
    gauge = "psi=(psi_abs+uAtPeak*Y)/(P*ell^2); u1=(u1_abs-uAtPeak)/(P*ell); u2=u2_abs/(P*ell)"
    metadata.update(sourceCheckpoint=str(frame["sourceCheckpoint"]), frame=str(frame_path),
                    a=a, P=P, ell=ell, U0=U0, gauge=gauge,
                    velocityDerivative="native saved/reconstructed solver derivative operators")
    return dict(name=f"native_tau{tau}", kind="finite_time_native", tau=float(frame["tau"]),
                x=x, y=y, Q=np.asarray(frame["Q"], float), q=np.asarray(frame["q"], float),
                psi=psi, psi_saved=psi_saved, u1=(ua-U0)/(P*ell), u2=va/(P*ell),
                a=a, P=P, ell=ell, U0=U0, rates=rates,
                source=str(raw_path), gauge=gauge, metadata=metadata)


def load_conditional(name: str) -> dict:
    """Return a finite-box conditional root with its own origin speed removed."""
    path = SOURCE / name / "state.npz"
    config_path = path.with_name("config.json")
    config = json.loads(config_path.read_text()) if config_path.exists() else {}
    stencil = int(config.get("stencil", 5))
    with np.load(path, allow_pickle=False) as raw:
        x, y, psi_saved = (np.asarray(raw[k], float) for k in ("x", "y", "psi"))
        Q, q, rates = (np.asarray(raw[k], float) for k in ("Q", "qx", "rates"))
    ua = -fd_matrix(y, 1, stencil) @ psi_saved
    va = psi_saved @ fd_matrix(x, 1, stencil).T
    va[0] = 0.0
    U0 = float(CubicSpline(x, ua[0])(0.0))
    psi = psi_saved + U0 * y[:, None]
    gauge = "psi=psi_saved+U0*y, U0=(-Dy@psi_saved)(x=0,y=0); u1=-Dy@psi_saved-U0"
    metadata = dict(config=config, derivativeStencil=stencil, U0=U0, gauge=gauge,
                    scope="conditional frozen-exterior finite-box root; no global BS certification",
                    a=0.0, P=1.0, ell=1.0,
                    normalization="coordinates already in the root's inner normalization")
    return dict(name=name, kind="conditional_frozen_exterior_root", x=x, y=y,
                Q=Q, q=q, psi=psi, psi_saved=psi_saved, u1=ua-U0, u2=va,
                a=0.0, P=1.0, ell=1.0, U0=U0, rates=rates,
                source=str(path), gauge=gauge, metadata=metadata)


def psi_over_y(profile: dict) -> np.ndarray:
    """Continuous wall extension of centred psi/y; no epsilon denominator.

    At the wall psi/y=psi_y=-u1. The saved solver derivative defines the
    wall trace; independent spline/finite-difference checks should report
    their mismatch instead of silently changing the centring convention.
    """
    y = np.asarray(profile["y"])
    out = np.empty_like(profile["psi"])
    out[0] = -profile["u1"][0]
    out[1:] = profile["psi"][1:] / y[1:, None]
    return out


def profiles():
    for tau in (10, 15, 20):
        yield load_native(tau)
    for name in CONDITIONAL_NAMES:
        yield load_conditional(name)

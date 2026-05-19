#!/usr/bin/env python
"""
composite_patch_tests.py

First-stage composite benchmark/sanity tests for:

  comp_dkmq24.CompDKMQ24
  comp_dkmt18.CompDKMT18

Purpose
-------
This is not a full structural composite benchmark yet.  It validates the most
important laminate mechanics before plate/shell examples:

  1. isotropic laminate reproduces isotropic base element stiffness
  1b. one-ply isotropic laminate also reproduces isotropic base element stiffness
  2. symmetric laminate has B ~= 0
  3. unsymmetric laminate has B != 0
  4. ABD energy for constant membrane/bending fields is correct
  5. unsymmetric stiffness matrix remains symmetric

Run
---
  python composite_patch_tests.py
  python composite_patch_tests.py --json composite_patch_summary.json

Interpretation
--------------
PASS here means the element can consume A/B/D/As laminate stiffness correctly at
single-element/constant-strain level.  It does NOT yet validate ply stress
recovery, failure criteria, buckling Kg, thermal loads, or Nastran PCOMP syntax.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import json
import math
import sys
from dataclasses import asdict
from typing import Any

import numpy as np

# Sandbox fallback: some uploaded workspaces store core.py as core(2).py.
HERE = Path(__file__).resolve().parent
if "core" not in sys.modules and not (HERE / "core.py").exists() and (HERE / "core(2).py").exists():
    import importlib.util
    spec = importlib.util.spec_from_file_location("core", HERE / "core(2).py")
    core_mod = importlib.util.module_from_spec(spec)
    sys.modules["core"] = core_mod
    spec.loader.exec_module(core_mod)

try:
    from core import Node
except Exception:
    from core_sparse import Node

from laminate_utils import Laminate
from comp_dkmq24 import CompDKMQ24
from comp_dkmt18 import CompDKMT18
from dkmq24_katili import DKMQ24
from dkmt18_maknun_snorm import DKMT18MaknunSNORM


def make_node(nid: int, x: float, y: float, z: float = 0.0):
    return Node(nid, x, y, z)


def maxabs(a: np.ndarray) -> float:
    return float(np.max(np.abs(np.asarray(a, dtype=float))))


def relerr(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    den = max(1.0, float(np.max(np.abs(b))))
    return maxabs(a - b) / den


def passfail(ok: bool) -> str:
    return "PASS" if ok else "FAIL"


def q4_nodes():
    return [
        make_node(1, 0.0, 0.0, 0.0),
        make_node(2, 1.0, 0.0, 0.0),
        make_node(3, 1.0, 1.0, 0.0),
        make_node(4, 0.0, 1.0, 0.0),
    ]


def t3_nodes():
    return [
        make_node(1, 0.0, 0.0, 0.0),
        make_node(2, 1.0, 0.0, 0.0),
        make_node(3, 0.0, 1.0, 0.0),
    ]


def sample_laminates():
    mat = dict(E1=135e9, E2=10e9, nu12=0.28, G12=5e9, G13=5e9, G23=3.8e9, t=0.000125)
    sym_crossply = Laminate.from_plies([
        dict(**mat, theta=0),
        dict(**mat, theta=90),
        dict(**mat, theta=90),
        dict(**mat, theta=0),
    ])

    unsym_crossply = Laminate.from_plies([
        dict(**mat, theta=0),
        dict(**mat, theta=90),
    ])

    angle_symmetric = Laminate.from_plies([
        dict(**mat, theta=45),
        dict(**mat, theta=-45),
        dict(**mat, theta=-45),
        dict(**mat, theta=45),
    ])

    return {
        "sym_crossply_0_90_90_0": sym_crossply,
        "unsym_crossply_0_90": unsym_crossply,
        "sym_angle_45_-45_-45_45": angle_symmetric,
    }


def abd_energy(area: float, lam: Laminate, eps: np.ndarray, kap: np.ndarray) -> float:
    x = np.r_[eps, kap]
    return 0.5 * area * float(x @ lam.ABD @ x)


def shear_energy(area: float, lam: Laminate, gam: np.ndarray) -> float:
    return 0.5 * area * float(gam @ lam.As @ gam)


def run_tests(args) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    def add(name, metric, value, tol, ok, note=""):
        rows.append({
            "test": name,
            "metric": metric,
            "value": float(value) if isinstance(value, (int, float, np.floating)) and math.isfinite(float(value)) else value,
            "tolerance": tol,
            "status": passfail(bool(ok)),
            "note": note,
        })

    E = 210e9
    nu = 0.3
    h = 0.01
    lam_iso = Laminate.isotropic(E, nu, h)
    G = E / (2.0 * (1.0 + nu))
    lam_oneply_iso = Laminate.from_plies([
        dict(E1=E, E2=E, nu12=nu, G12=G, G13=G, G23=G, theta=0, t=h),
    ])

    qnodes = q4_nodes()
    tnodes = t3_nodes()

    # 1. Isotropic reproduction
    q_iso = DKMQ24(1, qnodes, E, nu, h, use_mindlin=True).k_global()
    q_comp = CompDKMQ24(1, qnodes, lam_iso, use_mindlin=True).k_global()
    err_q = relerr(q_comp, q_iso)
    add("isotropic_reproduction_dkmq24", "relative_matrix_error", err_q, args.iso_tol, err_q <= args.iso_tol)
    add("isotropic_symmetry_dkmq24", "max_abs_K_minus_KT", maxabs(q_comp - q_comp.T), args.sym_tol, maxabs(q_comp - q_comp.T) <= args.sym_tol)

    t_iso = DKMT18MaknunSNORM(1, tnodes, E, nu, h).k_global()
    t_comp = CompDKMT18(1, tnodes, lam_iso).k_global()
    err_t = relerr(t_comp, t_iso)
    add("isotropic_reproduction_dkmt18", "relative_matrix_error", err_t, args.iso_tol, err_t <= args.iso_tol)
    add("isotropic_symmetry_dkmt18", "max_abs_K_minus_KT", maxabs(t_comp - t_comp.T), args.sym_tol, maxabs(t_comp - t_comp.T) <= args.sym_tol)

    q_comp_oneply = CompDKMQ24(11, qnodes, lam_oneply_iso, use_mindlin=True).k_global()
    err_q_oneply = relerr(q_comp_oneply, q_iso)
    add("oneply_iso_reproduction_dkmq24", "relative_matrix_error", err_q_oneply, args.iso_tol, err_q_oneply <= args.iso_tol)
    add("oneply_iso_symmetry_dkmq24", "max_abs_K_minus_KT", maxabs(q_comp_oneply - q_comp_oneply.T), args.sym_tol, maxabs(q_comp_oneply - q_comp_oneply.T) <= args.sym_tol)

    t_comp_oneply = CompDKMT18(11, tnodes, lam_oneply_iso).k_global()
    err_t_oneply = relerr(t_comp_oneply, t_iso)
    add("oneply_iso_reproduction_dkmt18", "relative_matrix_error", err_t_oneply, args.iso_tol, err_t_oneply <= args.iso_tol)
    add("oneply_iso_symmetry_dkmt18", "max_abs_K_minus_KT", maxabs(t_comp_oneply - t_comp_oneply.T), args.sym_tol, maxabs(t_comp_oneply - t_comp_oneply.T) <= args.sym_tol)

    # 2. Laminate B checks
    laminates = sample_laminates()
    for lname, lam in laminates.items():
        bmax = maxabs(lam.B)
        dmax = max(1.0, maxabs(lam.D))
        normalized_b = bmax * lam.h / dmax
        if "unsym" in lname:
            ok = bmax > args.unsym_b_min
            add(f"{lname}_B_nonzero", "max_abs_B", bmax, f"> {args.unsym_b_min:g}", ok)
        else:
            ok = normalized_b <= args.b_zero_tol
            add(f"{lname}_B_zero", "max_abs_B_times_h_over_maxD", normalized_b, args.b_zero_tol, ok)

    # 3. ABD energy algebra, independent of element B matrices.
    # Area = unit square for Q4; right triangle area=0.5 for T3.
    lam_u = laminates["unsym_crossply_0_90"]
    area_q = 1.0
    area_t = 0.5

    eps_cases = {
        "eps_x": np.array([1.0e-4, 0.0, 0.0]),
        "eps_y": np.array([0.0, 1.0e-4, 0.0]),
        "gam_xy": np.array([0.0, 0.0, 1.0e-4]),
        "mixed_membrane": np.array([1.0e-4, -0.5e-4, 0.25e-4]),
    }
    kap_cases = {
        "kap_x": np.array([2.0e-3, 0.0, 0.0]),
        "kap_y": np.array([0.0, -1.5e-3, 0.0]),
        "kap_xy": np.array([0.0, 0.0, 1.0e-3]),
        "mixed_bending": np.array([2.0e-3, -1.0e-3, 0.5e-3]),
    }

    for cname, eps in eps_cases.items():
        zero = np.zeros(3)
        U = abd_energy(area_q, lam_u, eps, zero)
        Uref = 0.5 * area_q * float(eps @ lam_u.A @ eps)
        err = abs(U - Uref) / max(1.0, abs(Uref))
        add(f"abd_membrane_energy_{cname}", "relative_error", err, args.energy_tol, err <= args.energy_tol)

    for cname, kap in kap_cases.items():
        zero = np.zeros(3)
        U = abd_energy(area_q, lam_u, zero, kap)
        Uref = 0.5 * area_q * float(kap @ lam_u.D @ kap)
        err = abs(U - Uref) / max(1.0, abs(Uref))
        add(f"abd_bending_energy_{cname}", "relative_error", err, args.energy_tol, err <= args.energy_tol)

    # Coupled energy includes eps.B.kap twice under 0.5*[eps,kap] ABD [eps,kap].
    eps = np.array([1.0e-4, -0.25e-4, 0.1e-4])
    kap = np.array([2.0e-3, -1.25e-3, 0.4e-3])
    U = abd_energy(area_q, lam_u, eps, kap)
    Uref = 0.5 * area_q * (
        float(eps @ lam_u.A @ eps)
        + 2.0 * float(eps @ lam_u.B @ kap)
        + float(kap @ lam_u.D @ kap)
    )
    err = abs(U - Uref) / max(1.0, abs(Uref))
    add("abd_membrane_bending_coupling_energy", "relative_error", err, args.energy_tol, err <= args.energy_tol)

    # 4. Unsymmetric element stiffness symmetry
    q_unsym = CompDKMQ24(2, qnodes, lam_u, use_mindlin=True).k_global()
    t_unsym = CompDKMT18(2, tnodes, lam_u).k_global()
    sq = maxabs(q_unsym - q_unsym.T)
    st = maxabs(t_unsym - t_unsym.T)
    add("unsymmetric_laminate_stiffness_symmetry_dkmq24", "max_abs_K_minus_KT", sq, args.sym_tol, sq <= args.sym_tol)
    add("unsymmetric_laminate_stiffness_symmetry_dkmt18", "max_abs_K_minus_KT", st, args.sym_tol, st <= args.sym_tol)

    # 5. Shear As positive and nonzero
    evals = np.linalg.eigvalsh(lam_u.As)
    min_ev = float(np.min(evals))
    add("laminate_shear_As_spd", "min_eigenvalue_As", min_ev, "> 0", min_ev > 0.0)

    return rows


def print_table(rows: list[dict[str, Any]]) -> None:
    print("Composite patch / laminate sanity tests")
    print("=" * 78)
    for r in rows:
        print(f"{r['status']:4s}  {r['test']:<52s} {r['metric']:<28s} {r['value']}")
    print("=" * 78)
    nfail = sum(1 for r in rows if r["status"] != "PASS")
    print(f"TOTAL: {len(rows)} tests, {nfail} failed")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iso-tol", type=float, default=1.0e-10)
    ap.add_argument("--sym-tol", type=float, default=1.0e-8)
    ap.add_argument("--energy-tol", type=float, default=1.0e-12)
    ap.add_argument("--b-zero-tol", type=float, default=1.0e-10)
    ap.add_argument("--unsym-b-min", type=float, default=1.0e-12)
    ap.add_argument("--json", default="")
    ap.add_argument("--csv", default="composite_patch_summary.csv")
    args = ap.parse_args()

    rows = run_tests(args)
    print_table(rows)

    try:
        import pandas as pd
        pd.DataFrame(rows).to_csv(args.csv, index=False)
        print(f"CSV: {args.csv}")
    except Exception:
        pass

    if args.json:
        Path(args.json).write_text(json.dumps(rows, indent=2), encoding="utf-8")
        print(f"JSON: {args.json}")

    if any(r["status"] != "PASS" for r in rows):
        raise SystemExit(1)


if __name__ == "__main__":
    main()

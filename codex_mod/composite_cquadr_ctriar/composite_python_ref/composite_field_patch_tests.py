#!/usr/bin/env python
"""
composite_field_patch_tests.py

Multi-element composite field patch tests for:

  comp_dkmq24
  comp_dkmt18

This assembles a real FE mesh, applies an exact prescribed displacement field,
and compares

    U_FE = 0.5 * u^T K u

against

    U_REF = 0.5 * area * [eps0,kappa]^T ABD [eps0,kappa]

Fields:
  ux = ex*x + 0.5*gxy*y
  uy = ey*y + 0.5*gxy*x
  w  = -0.5*kx*x^2 -0.5*ky*y^2 -0.5*kxy*x*y
  rx = dw/dy
  ry = -dw/dx
  rz = 0

Run:
  python composite_field_patch_tests.py
  python composite_field_patch_tests.py --case membrane
  python composite_field_patch_tests.py --case bending
  python composite_field_patch_tests.py --case coupled
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
import numpy as np

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

DOF = 6


def make_node(nid, x, y, z=0.0):
    return Node(nid, x, y, z)


def nid_ij(i, j, n):
    return j * (n + 1) + i + 1


def build_mesh(a, n):
    nodes = {}
    for j in range(n + 1):
        y = a * j / n
        for i in range(n + 1):
            x = a * i / n
            nid = nid_ij(i, j, n)
            nodes[nid] = make_node(nid, x, y, 0.0)

    quads = []
    tris = []
    eid = 1
    for j in range(n):
        for i in range(n):
            n1 = nid_ij(i, j, n)
            n2 = nid_ij(i + 1, j, n)
            n3 = nid_ij(i + 1, j + 1, n)
            n4 = nid_ij(i, j + 1, n)
            quads.append((eid, [n1, n2, n3, n4]))
            tris.append((2 * eid - 1, [n1, n2, n3]))
            tris.append((2 * eid, [n1, n3, n4]))
            eid += 1
    return nodes, quads, tris


def make_laminate(name, h=0.01):
    name = name.lower()
    if name == "isotropic":
        return Laminate.isotropic(E=70e9, nu=0.3, h=h)

    base = dict(E1=135e9, E2=10e9, nu12=0.28, G12=5e9, G13=5e9, G23=3.8e9)

    if name == "sym_crossply":
        t = h / 4.0
        return Laminate.from_plies([
            dict(**base, theta=0, t=t),
            dict(**base, theta=90, t=t),
            dict(**base, theta=90, t=t),
            dict(**base, theta=0, t=t),
        ])

    if name == "unsym_crossply":
        t = h / 2.0
        return Laminate.from_plies([
            dict(**base, theta=0, t=t),
            dict(**base, theta=90, t=t),
        ])

    if name == "sym_angle":
        t = h / 4.0
        return Laminate.from_plies([
            dict(**base, theta=45, t=t),
            dict(**base, theta=-45, t=t),
            dict(**base, theta=-45, t=t),
            dict(**base, theta=45, t=t),
        ])

    raise ValueError(f"unknown laminate {name}")


def make_elements(element, nodes, quads, tris, laminate):
    if element == "comp_dkmq24":
        return [CompDKMQ24(eid, [nodes[i] for i in conn], laminate) for eid, conn in quads]

    if element == "comp_dkmt18":
        normals = {nid: np.array([0.0, 0.0, 1.0], dtype=float) for nid in nodes}
        return [CompDKMT18(eid, [nodes[i] for i in conn], laminate, nodal_normals=normals) for eid, conn in tris]

    raise ValueError(f"unknown element {element}")


def assemble_K(nodes, elems):
    ndof = DOF * len(nodes)
    K = np.zeros((ndof, ndof), dtype=float)
    for elem in elems:
        ke = elem.k_global()
        edofs = []
        for node in elem.nodes:
            base = DOF * (node.nid - 1)
            edofs.extend([base + k for k in range(DOF)])
        for a, ia in enumerate(edofs):
            for b, ib in enumerate(edofs):
                K[ia, ib] += ke[a, b]
    return 0.5 * (K + K.T)


def displacement_vector(nodes, eps, kap):
    U = np.zeros(DOF * len(nodes), dtype=float)
    ex, ey, gxy = eps
    kx, ky, kxy = kap

    for nid, node in nodes.items():
        x, y, _ = node.coords
        ux = ex * x + 0.5 * gxy * y
        uy = ey * y + 0.5 * gxy * x

        w = -0.5 * kx * x * x - 0.5 * ky * y * y - 0.5 * kxy * x * y
        dwdx = -kx * x - 0.5 * kxy * y
        dwdy = -ky * y - 0.5 * kxy * x

        rx = dwdy
        ry = -dwdx
        rz = 0.0

        base = DOF * (nid - 1)
        U[base + 0] = ux
        U[base + 1] = uy
        U[base + 2] = w
        U[base + 3] = rx
        U[base + 4] = ry
        U[base + 5] = rz
    return U


def abd_energy(area, lam, eps, kap):
    x = np.r_[eps, kap]
    return 0.5 * area * float(x @ lam.ABD @ x)


def run_case(element, laminate_name, n, a, h, eps, kap, tol):
    lam = make_laminate(laminate_name, h=h)
    nodes, quads, tris = build_mesh(a, n)
    elems = make_elements(element, nodes, quads, tris, lam)
    K = assemble_K(nodes, elems)
    U = displacement_vector(nodes, eps, kap)
    Ufe = 0.5 * float(U @ K @ U)
    Uref = abd_energy(a * a, lam, eps, kap)
    rel = abs(Ufe - Uref) / max(1.0, abs(Uref))
    return {
        "element": element,
        "laminate": laminate_name,
        "n": n,
        "nodes": len(nodes),
        "elements": len(elems),
        "Ufe": Ufe,
        "Uref": Uref,
        "rel_error": rel,
        "Bmax": float(np.max(np.abs(lam.B))),
        "status": "PASS" if rel <= tol else "FAIL",
    }


def parse_ns(s):
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def vec3_csv(s):
    vals = [float(x.strip()) for x in s.split(",")]
    if len(vals) != 3:
        raise ValueError("expected three comma-separated numbers")
    return np.array(vals, dtype=float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--elements", nargs="*", default=["comp_dkmq24", "comp_dkmt18"])
    ap.add_argument("--laminates", nargs="*", default=["isotropic", "sym_crossply", "unsym_crossply", "sym_angle"])
    ap.add_argument("--ns", default="1,2,4")
    ap.add_argument("--a", type=float, default=1.0)
    ap.add_argument("--h", type=float, default=0.01)
    ap.add_argument("--tol", type=float, default=1.0e-6)
    ap.add_argument("--case", choices=["membrane", "bending", "coupled"], default="coupled")
    ap.add_argument("--eps", default="")
    ap.add_argument("--kap", default="")
    ap.add_argument("--csv", default="composite_field_patch_summary.csv")
    args = ap.parse_args()

    if args.eps:
        eps = vec3_csv(args.eps)
    elif args.case == "bending":
        eps = np.zeros(3)
    else:
        eps = np.array([1.0e-4, -0.5e-4, 0.25e-4])

    if args.kap:
        kap = vec3_csv(args.kap)
    elif args.case == "membrane":
        kap = np.zeros(3)
    else:
        kap = np.array([2.0e-3, -1.0e-3, 0.5e-3])

    rows = []
    print("element, laminate, case, n, nodes, elements, Ufe, Uref, rel_error, Bmax, status")
    for element in args.elements:
        for laminate in args.laminates:
            for n in parse_ns(args.ns):
                try:
                    r = run_case(element, laminate, n, args.a, args.h, eps, kap, args.tol)
                    rows.append({**r, "case": args.case})
                    print(
                        f"{r['element']}, {r['laminate']}, {args.case}, {r['n']}, {r['nodes']}, {r['elements']}, "
                        f"{r['Ufe']:.12e}, {r['Uref']:.12e}, {r['rel_error']:.12e}, "
                        f"{r['Bmax']:.12e}, {r['status']}"
                    )
                except Exception as exc:
                    print(f"{element}, {laminate}, {args.case}, {n}, NaN, NaN, NaN, NaN, NaN, NaN, ERROR: {type(exc).__name__}: {exc}")

    try:
        import pandas as pd
        pd.DataFrame(rows).to_csv(args.csv, index=False)
        print(f"CSV: {args.csv}")
    except Exception:
        pass

    if any(r.get("status") != "PASS" for r in rows):
        raise SystemExit(1)


if __name__ == "__main__":
    main()

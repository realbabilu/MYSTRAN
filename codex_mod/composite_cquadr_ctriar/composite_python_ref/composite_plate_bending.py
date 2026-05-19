#!/usr/bin/env python
"""
composite_plate_bending.py

Composite laminated square plate bending benchmark.

Elements:
  comp_dkmq24
  comp_dkmt18

Laminates:
  isotropic
  oneply_iso
  sym_crossply     [0/90/90/0]
  unsym_crossply   [0/90]
  sym_angle        [45/-45/-45/45]

Problem:
  square plate a x a
  simply supported boundary, transverse displacement w=0 on all edges
  rotations free
  minimal in-plane anchors to remove rigid body modes
  uniform transverse pressure q in global -Z

Purpose:
  This is a structural convergence benchmark, not an exact closed-form
  validation yet.  The important composite checks before this are in
  composite_patch_tests.py.  Here we check convergence and element-to-element
  consistency for laminated plates, including unsymmetric B coupling.

Run examples:
  python composite_plate_bending.py --element comp_dkmq24 --laminate sym_crossply --n 16
  python composite_plate_bending.py --element comp_dkmt18 --laminate sym_crossply --n 16
  python composite_plate_bending.py --all --sweep --ns 4,8,16,32
  python composite_plate_bending.py --all --laminate unsym_crossply --sweep --ns 4,8,16,32 > composite_plate_unsym.txt

Output columns:
  element, laminate, n, nodes, elements, center_w, corner_ux, corner_uy, max_abs_w
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Iterable

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
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

try:
    from core import Node
except Exception:
    from core_sparse import Node

from laminate_utils import Laminate
from comp_dkmq24 import CompDKMQ24
from comp_dkmt18 import CompDKMT18


DOF_PER_NODE = 6


def make_node(nid: int, x: float, y: float, z: float = 0.0):
    return Node(nid, x, y, z)


def make_laminate(name: str, h: float = 0.001) -> Laminate:
    """
    Return laminate in SI-like units.

    CFRP-ish default:
      E1=135 GPa, E2=10 GPa, nu12=0.28,
      G12=5 GPa, G13=5 GPa, G23=3.8 GPa.
    """
    name = name.lower()
    if name == "isotropic":
        return Laminate.isotropic(E=70e9, nu=0.3, h=h)

    if name == "oneply_iso":
        E = 70e9
        nu = 0.3
        G = E / (2.0 * (1.0 + nu))
        return Laminate.from_plies([
            dict(E1=E, E2=E, nu12=nu, G12=G, G13=G, G23=G, theta=0, t=h),
        ])

    # Use ply thickness so total h is close to requested h.
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

    raise ValueError(f"unknown laminate: {name}")


def nid_ij(i: int, j: int, n: int) -> int:
    return j * (n + 1) + i + 1


def build_mesh(a: float, n: int):
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
            # split along diagonal n1-n3
            tris.append((2 * eid - 1, [n1, n2, n3]))
            tris.append((2 * eid, [n1, n3, n4]))
            eid += 1

    return nodes, quads, tris


def polygon_area_xy(coords: np.ndarray) -> float:
    x = coords[:, 0]
    y = coords[:, 1]
    return 0.5 * abs(float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def make_elements(element: str, nodes: dict[int, Node], quads, tris, laminate: Laminate):
    element = element.lower()
    elems = []

    if element == "comp_dkmq24":
        for eid, conn in quads:
            elems.append(CompDKMQ24(eid, [nodes[i] for i in conn], laminate))
        return elems

    if element == "comp_dkmt18":
        normal_dict = {nid: np.array([0.0, 0.0, 1.0], dtype=float) for nid in nodes}
        for eid, conn in tris:
            elems.append(CompDKMT18(eid, [nodes[i] for i in conn], laminate, nodal_normals=normal_dict))
        return elems

    raise ValueError(f"unknown element: {element}")


def assemble(nodes: dict[int, Node], elems, pressure: float):
    ndof = DOF_PER_NODE * len(nodes)
    F = np.zeros(ndof, dtype=float)

    if HAVE_SCIPY:
        rows = []
        cols = []
        vals = []
    else:
        K = np.zeros((ndof, ndof), dtype=float)

    for elem in elems:
        ke = elem.k_global()
        conn = [n.nid for n in elem.nodes]
        edofs = []
        for nid in conn:
            base = DOF_PER_NODE * (nid - 1)
            edofs.extend([base + k for k in range(DOF_PER_NODE)])

        # stiffness
        if HAVE_SCIPY:
            for a, ia in enumerate(edofs):
                for b, ib in enumerate(edofs):
                    v = ke[a, b]
                    if v != 0.0:
                        rows.append(ia)
                        cols.append(ib)
                        vals.append(v)
        else:
            for a, ia in enumerate(edofs):
                for b, ib in enumerate(edofs):
                    K[ia, ib] += ke[a, b]

        # simple consistent nodal pressure by element area share
        coords = np.array([nodes[nid].coords for nid in conn], dtype=float)
        area = polygon_area_xy(coords)
        fnode = -pressure * area / len(conn)
        for nid in conn:
            F[DOF_PER_NODE * (nid - 1) + 2] += fnode

    if HAVE_SCIPY:
        K = sp.coo_matrix((vals, (rows, cols)), shape=(ndof, ndof)).tocsr()
        K = 0.5 * (K + K.T)
    else:
        K = 0.5 * (K + K.T)

    return K, F


def boundary_conditions(a: float, n: int, nodes: dict[int, Node]):
    """
    Simply supported:
      w=0 on all edges.

    Minimal in-plane anchors:
      corner (0,0): ux=uy=0
      corner (a,0): uy=0

    This lets unsymmetric laminate extension-bending coupling develop while
    preventing in-plane rigid body modes.
    """
    bcs = {}

    def add(nid, dof, val=0.0):
        bcs[DOF_PER_NODE * (nid - 1) + dof] = val

    # w = 0 on all boundary nodes
    for j in range(n + 1):
        for i in range(n + 1):
            if i == 0 or i == n or j == 0 or j == n:
                add(nid_ij(i, j, n), 2, 0.0)

    # in-plane anchors
    add(nid_ij(0, 0, n), 0, 0.0)
    add(nid_ij(0, 0, n), 1, 0.0)
    add(nid_ij(n, 0, n), 1, 0.0)

    return bcs


def solve_system(K, F, bcs: dict[int, float]):
    ndof = len(F)
    fixed = np.array(sorted(bcs.keys()), dtype=int)
    fixed_vals = np.array([bcs[i] for i in fixed], dtype=float)
    allidx = np.arange(ndof, dtype=int)
    free = np.setdiff1d(allidx, fixed)

    U = np.zeros(ndof, dtype=float)
    U[fixed] = fixed_vals

    if HAVE_SCIPY and sp.issparse(K):
        rhs = F[free] - K[free][:, fixed] @ U[fixed]
        U[free] = spla.spsolve(K[free][:, free], rhs)
    else:
        rhs = F[free] - K[np.ix_(free, fixed)] @ U[fixed]
        U[free] = np.linalg.solve(K[np.ix_(free, free)], rhs)

    return U


def run_case(element: str, laminate_name: str, n: int, a: float, h: float, pressure: float, verbose: bool = True):
    laminate = make_laminate(laminate_name, h=h)
    nodes, quads, tris = build_mesh(a, n)
    elems = make_elements(element, nodes, quads, tris, laminate)
    K, F = assemble(nodes, elems, pressure)
    bcs = boundary_conditions(a, n, nodes)
    U = solve_system(K, F, bcs)

    center_nid = nid_ij(n // 2, n // 2, n) if n % 2 == 0 else None
    if center_nid is None:
        # nearest node to center
        target = np.array([a / 2.0, a / 2.0, 0.0])
        center_nid = min(nodes, key=lambda nid: float(np.linalg.norm(nodes[nid].coords - target)))

    center_w = U[DOF_PER_NODE * (center_nid - 1) + 2]
    w_all = U[2::DOF_PER_NODE]
    max_abs_w = float(np.max(np.abs(w_all)))

    corner = nid_ij(n, n, n)
    corner_ux = U[DOF_PER_NODE * (corner - 1) + 0]
    corner_uy = U[DOF_PER_NODE * (corner - 1) + 1]

    if verbose:
        print("Composite square plate bending")
        print(f"element       : {element}")
        print(f"laminate      : {laminate_name}")
        print(f"n             : {n}")
        print(f"nodes         : {len(nodes)}")
        print(f"elements      : {len(elems)}")
        print(f"a, h, pressure: {a}, {h}, {pressure}")
        print(f"laminate Bmax : {float(np.max(np.abs(laminate.B))):.6e}")
        print(f"center node   : {center_nid}")
        print(f"center_w      : {center_w:.12e}")
        print(f"max_abs_w     : {max_abs_w:.12e}")
        print(f"corner_ux     : {corner_ux:.12e}")
        print(f"corner_uy     : {corner_uy:.12e}")

    return {
        "element": element,
        "laminate": laminate_name,
        "n": n,
        "nodes": len(nodes),
        "elements": len(elems),
        "center_w": center_w,
        "max_abs_w": max_abs_w,
        "corner_ux": corner_ux,
        "corner_uy": corner_uy,
        "Bmax": float(np.max(np.abs(laminate.B))),
    }


def parse_ns(s: str):
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--element", choices=["comp_dkmq24", "comp_dkmt18"], default="comp_dkmq24")
    ap.add_argument("--laminate", choices=["isotropic", "oneply_iso", "sym_crossply", "unsym_crossply", "sym_angle"], default="sym_crossply")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--sweep", action="store_true")
    ap.add_argument("--ns", default="4,8,16,32")
    ap.add_argument("--n", type=int, default=16)
    ap.add_argument("--a", type=float, default=1.0)
    ap.add_argument("--h", type=float, default=0.01)
    ap.add_argument("--pressure", type=float, default=1000.0)
    ap.add_argument("--csv", default="")
    args = ap.parse_args()

    elems = ["comp_dkmq24", "comp_dkmt18"] if args.all else [args.element]
    ns = parse_ns(args.ns) if args.sweep else [args.n]

    results = []

    if args.sweep or args.all:
        print("element, laminate, n, nodes, elements, center_w, max_abs_w, corner_ux, corner_uy, Bmax")
        for elem in elems:
            for n in ns:
                try:
                    r = run_case(elem, args.laminate, n, args.a, args.h, args.pressure, verbose=False)
                    results.append(r)
                    print(
                        f"{r['element']}, {r['laminate']}, {r['n']}, {r['nodes']}, {r['elements']}, "
                        f"{r['center_w']:.12e}, {r['max_abs_w']:.12e}, "
                        f"{r['corner_ux']:.12e}, {r['corner_uy']:.12e}, {r['Bmax']:.12e}"
                    )
                except Exception as exc:
                    print(f"{elem}, {args.laminate}, {n}, NaN, NaN, NaN, NaN, NaN, NaN, ERROR: {type(exc).__name__}: {exc}")
    else:
        r = run_case(args.element, args.laminate, args.n, args.a, args.h, args.pressure, verbose=True)
        results.append(r)

    if args.csv:
        try:
            import pandas as pd
            pd.DataFrame(results).to_csv(args.csv, index=False)
            print(f"CSV: {args.csv}")
        except Exception as exc:
            print(f"WARNING: failed to write csv: {exc}")


if __name__ == "__main__":
    main()

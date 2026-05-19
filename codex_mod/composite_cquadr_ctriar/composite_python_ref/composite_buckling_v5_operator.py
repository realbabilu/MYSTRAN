#!/usr/bin/env python
"""
composite_buckling_v5_operator.py

Composite plate buckling benchmark v5 operator.

Compared to composite_eigen_buckling_v2.py buckling:
  v2: KG uses user-given Nx, Ny, Nxy directly.
  v3: KG uses laminate prebuckling resultants recovered from prescribed
      membrane/curvature fields:

        N = A eps0 + B kappa
        M = B eps0 + D kappa

Then solves:
        K phi = lambda KG phi

Interpretation
--------------
If eps_pre is prescribed, lambda is a multiplier on that prescribed
prebuckling strain/resultant state.

Example:
  --eps -1e-5,0,0 means eps_x = -1e-5 compression.
  lambda_crit = 120 means critical eps_x ~= -120e-5.

This is still a benchmark-level plate buckling KG, not a full MYSTRAN
prebuckling static solve with element-by-element stress recovery.

V5 operator method:
  KG may be singular and may act only on transverse w DOFs.  Do not delete
  rotations and do not form dense Schur complements.  Instead solve the
  reciprocal eigenproblem:
      OP(x) = K^{-1} KG x
      OP phi = mu phi
      lambda = 1 / mu
  This keeps all free DOFs coupled through K and is the production-oriented
  route for sparse/banded solvers.
"""

from __future__ import annotations

import argparse
import math
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
    import scipy.linalg as la
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


DOF = 6


def make_node(nid, x, y, z=0.0):
    return Node(nid, x, y, z)


def nid_ij(i, j, n):
    return j * (n + 1) + i + 1


def make_laminate(name, h=0.01):
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


def make_elements(element, nodes, quads, tris, laminate):
    if element == "comp_dkmq24":
        return [CompDKMQ24(eid, [nodes[i] for i in conn], laminate) for eid, conn in quads], quads

    if element == "comp_dkmt18":
        normals = {nid: np.array([0.0, 0.0, 1.0]) for nid in nodes}
        return [CompDKMT18(eid, [nodes[i] for i in conn], laminate, nodal_normals=normals) for eid, conn in tris], tris

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


def q4_dN_nat(r, s):
    return 0.25 * np.array([
        [-(1-s), +(1-s), +(1+s), -(1+s)],
        [-(1-r), -(1+r), +(1+r), +(1-r)],
    ], dtype=float)


def q4_geom_stiffness(coords, Nmat):
    gp = 1.0 / math.sqrt(3.0)
    Kww = np.zeros((4, 4), dtype=float)
    xy = coords[:, :2]
    for r in (-gp, gp):
        for s in (-gp, gp):
            dN = q4_dN_nat(r, s)
            J = dN @ xy
            detJ = float(np.linalg.det(J))
            dNxy = np.linalg.inv(J) @ dN
            Kww += detJ * (dNxy.T @ Nmat @ dNxy)
    return 0.5 * (Kww + Kww.T)


def tri_geom_stiffness(coords, Nmat):
    xy = coords[:, :2]
    x1, y1 = xy[0]
    x2, y2 = xy[1]
    x3, y3 = xy[2]
    A2 = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
    area = 0.5 * abs(A2)
    if abs(A2) < 1e-30:
        raise ValueError("zero-area triangle")
    b = np.array([y2-y3, y3-y1, y1-y2], dtype=float) / A2
    c = np.array([x3-x2, x1-x3, x2-x1], dtype=float) / A2
    G = np.vstack([b, c])
    return area * (G.T @ Nmat @ G)


def laminate_resultants(lam: Laminate, eps: np.ndarray, kap: np.ndarray):
    eps = np.asarray(eps, dtype=float).reshape(3)
    kap = np.asarray(kap, dtype=float).reshape(3)
    N = lam.A @ eps + lam.B @ kap
    M = lam.B @ eps + lam.D @ kap
    return N, M


def assemble_KG_from_resultants(nodes, conns, N):
    """
    N = [Nxx, Nyy, Nxy]. Compression convention:
    If eps_x is negative compression, Nxx is negative.
    KG for K = lambda KG should be positive for compressive state,
    therefore use -N in geometric stiffness.
    """
    ndof = DOF * len(nodes)
    KG = np.zeros((ndof, ndof), dtype=float)

    Nxx, Nyy, Nxy = N
    Nmat = -np.array([[Nxx, Nxy], [Nxy, Nyy]], dtype=float)

    for eid, conn in conns:
        coords = np.array([nodes[i].coords for i in conn], dtype=float)
        if len(conn) == 4:
            kgww = q4_geom_stiffness(coords, Nmat)
        else:
            kgww = tri_geom_stiffness(coords, Nmat)

        for a, na in enumerate(conn):
            ia = DOF * (na - 1) + 2
            for b, nb in enumerate(conn):
                ib = DOF * (nb - 1) + 2
                KG[ia, ib] += kgww[a, b]

    return 0.5 * (KG + KG.T)


def boundary_dofs(n):
    fixed = set()

    def add(nid, dof):
        fixed.add(DOF * (nid - 1) + dof)

    for j in range(n + 1):
        for i in range(n + 1):
            if i == 0 or i == n or j == 0 or j == n:
                add(nid_ij(i, j, n), 2)

    add(nid_ij(0, 0, n), 0)
    add(nid_ij(0, 0, n), 1)
    add(nid_ij(n, 0, n), 1)

    return sorted(fixed)


def solve_buckling(K, KG, fixed, nmodes):
    """
    Buckling by reciprocal operator:

        K phi = lambda KG phi
        K^{-1} KG phi = mu phi
        lambda = 1 / mu

    This avoids dense Schur condensation and keeps rotations/in-plane DOFs
    coupled through the full stiffness solve.  It is the closest Python
    analogue to a MYSTRAN/ARPACK reverse-communication implementation:

        y = KG*x
        z = solve(K,y)
        return z
    """
    ndof = K.shape[0]
    fixed = np.array(fixed, dtype=int)
    free = np.setdiff1d(np.arange(ndof), fixed)

    Kf = K[np.ix_(free, free)]
    KGf = KG[np.ix_(free, free)]

    nfree = Kf.shape[0]
    k_ask = min(max(nmodes + 6, 12), max(1, nfree - 2))

    if HAVE_SCIPY and nfree > 4:
        # Sparse factorization of dense arrays converted to CSC.
        # In MYSTRAN this maps to banded/sparse factor K once.
        Ksp = sp.csc_matrix(Kf)
        KGsp = sp.csr_matrix(KGf)
        solveK = spla.factorized(Ksp)

        def matvec(x):
            return solveK(KGsp @ x)

        OP = spla.LinearOperator((nfree, nfree), matvec=matvec, dtype=float)

        try:
            # OP is generally nonsymmetric in Euclidean inner product.
            vals = spla.eigs(OP, k=k_ask, which="LR", return_eigenvectors=False, tol=1e-9, maxiter=max(2000, 20*nfree))
        except Exception:
            vals = spla.eigs(OP, k=min(nmodes + 2, max(1, nfree - 2)), which="LM", return_eigenvectors=False, tol=1e-8, maxiter=max(3000, 30*nfree))

        vals = np.asarray(vals)
        vals = vals[np.isfinite(vals)]
        vals = vals[np.abs(vals.imag) <= 1e-6 * np.maximum(1.0, np.abs(vals.real))]
        mu = vals.real
    else:
        # Dense fallback for small/debug use.
        vals = np.linalg.eigvals(np.linalg.solve(Kf, KGf))
        vals = np.asarray(vals)
        vals = vals[np.isfinite(vals)]
        vals = vals[np.abs(vals.imag) <= 1e-6 * np.maximum(1.0, np.abs(vals.real))]
        mu = vals.real

    # Positive mu corresponds to positive buckling lambda.
    mu = mu[mu > 1e-14]
    if mu.size == 0:
        return np.array([], dtype=float)

    # Largest positive mu gives smallest positive lambda.
    lam = 1.0 / mu
    lam = lam[np.isfinite(lam)]
    lam = lam[lam > 1e-10]
    lam.sort()
    return lam[:nmodes]


def run_case(element, laminate_name, n, a, h, eps, kap, nmodes):
    lam = make_laminate(laminate_name, h=h)
    nodes, quads, tris = build_mesh(a, n)
    elems, conns = make_elements(element, nodes, quads, tris, lam)
    K = assemble_K(nodes, elems)
    N, M = laminate_resultants(lam, eps, kap)
    KG = assemble_KG_from_resultants(nodes, conns, N)
    fixed = boundary_dofs(n)
    eig = solve_buckling(K, KG, fixed, nmodes)
    return eig, N, M


def parse_ns(s):
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def vec3(s):
    vals = [float(x.strip()) for x in s.split(",")]
    if len(vals) != 3:
        raise ValueError("expected 3 comma-separated values")
    return np.array(vals, dtype=float)


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
    ap.add_argument("--eps", default="-1e-5,0,0")
    ap.add_argument("--kap", default="0,0,0")
    ap.add_argument("--nmodes", type=int, default=6)
    args = ap.parse_args()

    elems = ["comp_dkmq24", "comp_dkmt18"] if args.all else [args.element]
    ns = parse_ns(args.ns) if args.sweep else [args.n]
    eps = vec3(args.eps)
    kap = vec3(args.kap)

    print("analysis, element, laminate, n, mode, lambda, eps_multiplier, Nxx, Nyy, Nxy, Mxx, Myy, Mxy")
    for elem in elems:
        for n in ns:
            try:
                eig, N, M = run_case(elem, args.laminate, n, args.a, args.h, eps, kap, args.nmodes)
                for i, lamval in enumerate(eig, start=1):
                    print(
                        f"buckling_v5_operator, {elem}, {args.laminate}, {n}, {i}, "
                        f"{lamval:.12e}, {lamval:.12e}, "
                        f"{N[0]:.12e}, {N[1]:.12e}, {N[2]:.12e}, "
                        f"{M[0]:.12e}, {M[1]:.12e}, {M[2]:.12e}"
                    )
            except Exception as exc:
                print(f"buckling_v5_operator, {elem}, {args.laminate}, {n}, NaN, NaN, ERROR: {type(exc).__name__}: {exc}")


if __name__ == "__main__":
    main()

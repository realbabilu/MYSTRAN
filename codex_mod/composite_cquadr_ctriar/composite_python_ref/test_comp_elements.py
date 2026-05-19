#!/usr/bin/env python
"""
test_comp_elements.py

Quick checks for composite DKMT18/DKMQ24 wrappers.

Checks:
  1. Isotropic laminate reproduces isotropic DKMQ24/DKMT18 stiffness.
  2. One-ply isotropic laminate also reproduces isotropic DKMQ24/DKMT18 stiffness.
  2. Symmetric laminate has B approximately zero.
  3. Unsymmetric laminate has B nonzero and element stiffness remains symmetric.

Run:
  python test_comp_elements.py
"""

from __future__ import annotations

# Sandbox fallback: some uploaded workspaces store core.py as core(2).py.
from pathlib import Path as _Path
import sys as _sys
_HERE = _Path(__file__).resolve().parent
if "core" not in _sys.modules and not (_HERE / "core.py").exists() and (_HERE / "core(2).py").exists():
    import importlib.util as _importlib_util
    _spec = _importlib_util.spec_from_file_location("core", _HERE / "core(2).py")
    _core_mod = _importlib_util.module_from_spec(_spec)
    _sys.modules["core"] = _core_mod
    _spec.loader.exec_module(_core_mod)

import numpy as np

try:
    from core import Node
except Exception:
    from core_sparse import Node

from laminate_utils import Laminate
from comp_dkmq24 import CompDKMQ24
from comp_dkmt18 import CompDKMT18
from dkmq24_katili import DKMQ24
from dkmt18_maknun_snorm import DKMT18MaknunSNORM


def make_node(nid, x, y, z=0.0):
    try:
        return Node(nid, x, y, z)
    except TypeError:
        n = Node(nid, x, y, z)
        return n


def maxabs(a):
    return float(np.max(np.abs(a)))


def relerr(a, b):
    den = max(1.0, float(np.max(np.abs(b))))
    return maxabs(a - b) / den


def main():
    E = 210e9
    nu = 0.3
    h = 0.01
    lam_iso = Laminate.isotropic(E, nu, h)
    G = E / (2.0 * (1.0 + nu))
    lam_oneply_iso = Laminate.from_plies([
        dict(E1=E, E2=E, nu12=nu, G12=G, G13=G, G23=G, theta=0, t=h),
    ])

    # Q4 unit square
    qnodes = [
        make_node(1, 0.0, 0.0, 0.0),
        make_node(2, 1.0, 0.0, 0.0),
        make_node(3, 1.0, 1.0, 0.0),
        make_node(4, 0.0, 1.0, 0.0),
    ]
    q_iso = DKMQ24(1, qnodes, E, nu, h, use_mindlin=True).k_global()
    q_comp = CompDKMQ24(1, qnodes, lam_iso, use_mindlin=True).k_global()
    print("DKMQ24 isotropic reproduction relerr:", f"{relerr(q_comp, q_iso):.3e}")
    print("DKMQ24 comp symmetry:", f"{maxabs(q_comp-q_comp.T):.3e}")
    q_comp_oneply = CompDKMQ24(11, qnodes, lam_oneply_iso, use_mindlin=True).k_global()
    print("DKMQ24 oneply_iso reproduction relerr:", f"{relerr(q_comp_oneply, q_iso):.3e}")

    # T3 right triangle
    tnodes = [
        make_node(1, 0.0, 0.0, 0.0),
        make_node(2, 1.0, 0.0, 0.0),
        make_node(3, 0.0, 1.0, 0.0),
    ]
    t_iso = DKMT18MaknunSNORM(1, tnodes, E, nu, h).k_global()
    t_comp = CompDKMT18(1, tnodes, lam_iso).k_global()
    print("DKMT18 isotropic reproduction relerr:", f"{relerr(t_comp, t_iso):.3e}")
    print("DKMT18 comp symmetry:", f"{maxabs(t_comp-t_comp.T):.3e}")
    t_comp_oneply = CompDKMT18(11, tnodes, lam_oneply_iso).k_global()
    print("DKMT18 oneply_iso reproduction relerr:", f"{relerr(t_comp_oneply, t_iso):.3e}")

    # Symmetric and unsymmetric laminate checks
    mat = dict(E1=135e9, E2=10e9, nu12=0.28, G12=5e9, G13=5e9, G23=3.8e9, t=0.000125)
    lam_sym = Laminate.from_plies([
        dict(**mat, theta=0),
        dict(**mat, theta=90),
        dict(**mat, theta=90),
        dict(**mat, theta=0),
    ])
    lam_unsym = Laminate.from_plies([
        dict(**mat, theta=0),
        dict(**mat, theta=90),
    ])

    print("Symmetric laminate max |B|:", f"{maxabs(lam_sym.B):.3e}")
    print("Unsymmetric laminate max |B|:", f"{maxabs(lam_unsym.B):.3e}")

    q_unsym = CompDKMQ24(2, qnodes, lam_unsym).k_global()
    t_unsym = CompDKMT18(2, tnodes, lam_unsym).k_global()
    print("DKMQ24 unsym stiffness symmetry:", f"{maxabs(q_unsym-q_unsym.T):.3e}")
    print("DKMT18 unsym stiffness symmetry:", f"{maxabs(t_unsym-t_unsym.T):.3e}")


if __name__ == "__main__":
    main()

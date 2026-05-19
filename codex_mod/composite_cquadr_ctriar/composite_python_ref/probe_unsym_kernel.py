"""
probe_unsym_kernel.py

Small diagnostic to compare symmetric vs unsymmetric composite laminate inputs
for the Python composite shell kernels.  This helps localize whether the
remaining unsymmetric buckling mismatch is likely driven by:

1. preload/resultant state,
2. laminate ABD data, or
3. the element stiffness path under B-coupling.
"""

from __future__ import annotations

import numpy as np

from composite_buckling_v5_operator import make_laminate
from comp_dkmq24 import CompDKMQ24
from comp_dkmt18 import CompDKMT18
from core import Node


def frob(x: np.ndarray) -> float:
    return float(np.sqrt(np.sum(x * x)))


def maxabs(x: np.ndarray) -> float:
    return float(np.max(np.abs(x)))


def make_nodes_quad():
    return [
        Node(1, 0.0, 0.0, 0.0),
        Node(2, 1.0, 0.0, 0.0),
        Node(3, 1.0, 1.0, 0.0),
        Node(4, 0.0, 1.0, 0.0),
    ]


def make_nodes_tri():
    return [
        Node(1, 0.0, 0.0, 0.0),
        Node(2, 1.0, 0.0, 0.0),
        Node(3, 0.0, 1.0, 0.0),
    ]


def print_mat(name: str, M: np.ndarray):
    print(name)
    with np.printoptions(precision=6, suppress=False):
        print(M)


def main():
    lam_sym = make_laminate("sym_crossply")
    lam_unsym = make_laminate("unsym_crossply")

    print("Laminate comparison")
    print(f"sym   h = {lam_sym.h:.12e}")
    print(f"unsym h = {lam_unsym.h:.12e}")
    print(f"||A_sym - A_unsym||_F = {frob(lam_sym.A - lam_unsym.A):.12e}")
    print(f"||B_sym - B_unsym||_F = {frob(lam_sym.B - lam_unsym.B):.12e}")
    print(f"||D_sym - D_unsym||_F = {frob(lam_sym.D - lam_unsym.D):.12e}")
    print(f"max|B_sym|   = {maxabs(lam_sym.B):.12e}")
    print(f"max|B_unsym| = {maxabs(lam_unsym.B):.12e}")
    print()

    print_mat("A_sym", lam_sym.A)
    print_mat("A_unsym", lam_unsym.A)
    print_mat("B_unsym", lam_unsym.B)
    print_mat("D_sym", lam_sym.D)
    print_mat("D_unsym", lam_unsym.D)
    print()

    q_sym = CompDKMQ24(1, make_nodes_quad(), lam_sym)
    q_unsym = CompDKMQ24(1, make_nodes_quad(), lam_unsym)
    t_sym = CompDKMT18(1, make_nodes_tri(), lam_sym)
    t_unsym = CompDKMT18(1, make_nodes_tri(), lam_unsym)

    Kq_sym = q_sym.k_global()
    Kq_unsym = q_unsym.k_global()
    Kt_sym = t_sym.k_global()
    Kt_unsym = t_unsym.k_global()

    print("Element stiffness comparison")
    print(f"DKMQ24 ||K_sym - K_unsym||_F = {frob(Kq_sym - Kq_unsym):.12e}")
    print(f"DKMT18 ||K_sym - K_unsym||_F = {frob(Kt_sym - Kt_unsym):.12e}")
    print(f"DKMQ24 max|K_sym - K_unsym| = {maxabs(Kq_sym - Kq_unsym):.12e}")
    print(f"DKMT18 max|K_sym - K_unsym| = {maxabs(Kt_sym - Kt_unsym):.12e}")
    print(f"DKMQ24 ||K_unsym||_F        = {frob(Kq_unsym):.12e}")
    print(f"DKMT18 ||K_unsym||_F        = {frob(Kt_unsym):.12e}")


if __name__ == "__main__":
    main()

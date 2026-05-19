import sys
from pathlib import Path

import numpy as np

REF_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(REF_DIR))

from composite_buckling_v5_operator import (  # noqa: E402
    build_mesh,
    assemble_K,
    make_elements,
    make_laminate,
    laminate_resultants,
    assemble_KG_from_resultants,
    solve_buckling,
)


def fixed_edge_dofs(nodes, a):
    fixed = []
    for i, nd in enumerate(nodes.values()):
        x = nd.x
        y = nd.y
        if abs(x) < 1e-12 or abs(x - a) < 1e-12 or abs(y) < 1e-12 or abs(y - a) < 1e-12:
            for d in range(6):
                fixed.append(6 * i + d)
    return fixed


def main():
    if len(sys.argv) != 3:
        print("usage: py -3 probe_composite_buckling_model.py <comp_dkmq24|comp_dkmt18> <sym_crossply|unsym_crossply>")
        raise SystemExit(2)
    elem_name = sys.argv[1]
    lam_name = sys.argv[2]
    n = 2
    a = 1.0
    h = 0.01
    eps = np.array([-1.0e-5, 0.0, 0.0], dtype=float)
    kap = np.array([0.0, 0.0, 0.0], dtype=float)

    lam = make_laminate(lam_name, h)
    nodes, quads, tris = build_mesh(a, n)
    elems, conns = make_elements(elem_name, nodes, quads, tris, lam)
    K = assemble_K(nodes, elems)

    N, M = laminate_resultants(lam, eps, kap)
    KG = assemble_KG_from_resultants(nodes, conns, N)
    fixed = fixed_edge_dofs(nodes, a)
    eig = solve_buckling(K, KG, fixed, 6)

    print(f"element={elem_name}")
    print(f"laminate={lam_name}")
    print(f"K_norm={np.linalg.norm(K):.12e}")
    print(f"KG_norm={np.linalg.norm(KG):.12e}")
    print(f"N={N}")
    print(f"M={M}")
    print("eig=" + ",".join(f"{v:.12e}" for v in eig))


if __name__ == "__main__":
    main()

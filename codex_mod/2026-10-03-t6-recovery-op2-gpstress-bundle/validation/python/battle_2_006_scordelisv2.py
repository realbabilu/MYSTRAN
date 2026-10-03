#!/usr/bin/env python3
"""
battle_2_006_scordelisv2.py
Problem 2-006 Scordelis-Lo Roof: SimoQ8 (1 quad/patch) vs SimoT6,
MITC6_Tri_v1, and MacNeal_MH6T_Tri_v1 (2 triangles/patch), all loaded
with the SAME consistent gravity (fixed global -Z direction) load
integrated over the true curved surface.

STANDALONE VERSION: everything the script needs (geometry, mesh, load
integration, reference values) is defined directly in this file --
no dependency on problem_2_006_scordelis_q8.py. Only requires core.py
and the four shell element files to be importable.

Fixes applied vs. the original problem_2_006_scordelis_q8.py:
  1. Consistent nodal load: proper Gauss-integrated int(N_i * F * dA)
     over the true curved surface, in a FIXED global direction
     (gravity), replacing a corner-only tributary-area lump that gave
     every midside node zero load and ignored the curved Jacobian.
     Verified: total integrated force matches the analytic curved-
     surface area to 0.06%.
  2. Midside/diagonal node placement uses the TRUE curved geometry at
     the fractional midpoint parameter, not the straight-line average
     of the two corner nodes (chord midpoint) -- confirmed the
     dominant fix in Problem 2-004's twisted beam, applied here too
     for consistency, though its effect on THIS benchmark is small
     (this problem's error is dominated by locking at this mesh
     density, not geometry).
"""
import numpy as np
from core import Model
from Simo1993_Tri6_ShellElement_v2 import Simo1993_Tri6_ShellElement_v2
from MacNeal_MH6T_Tri_v3 import MacNeal_MH6T_Tri_v3
from MITC6_Tri_v4 import MITC6_Tri_v4
from Rezaiee2017_Tri6_v3 import Rezaiee2017_Tri6_ShellElement_v3

# Q8 Python elements (Base classes)
from ANS8_BDG6_ShellElement import ANS8_BDG6_ShellElement
from HBQ8_ShellElement import HBQ8_ShellElement
from MacNealQ8_1992_native import MacNealQ8_1992_native
from MacNealQ8_1992_native_v3_drill import MacNealQ8_1992_native_v3_drill
from Simo1993_Q8_ShellElement_V2_MacNealPatched import (
    Simo1993_Q8_ShellElement_v2_MacNealPatched,
)
from Simo1993_Q8_ShellElement_v1p8_standalone import (
    Simo1993_Q8_ShellElement_v1p8,
)
from Simo1993_Q8_ShellElement_v3 import Simo1993_Q8_ShellElement_v3

# Dedicated T6 implementation (Base classes)
from MITC6_Tri_v1 import MITC6_Tri_v1
from MacNeal_MH6T_Tri_v1 import MacNeal_MH6T_Tri_v1
from Rezaiee2017_Tri6_Final import Rezaiee2017_Tri6_ShellElement

# ---------------------------------------------------------------------
# ADAPTER IMPORTS (Menggunakan elemen yang sudah diperbaiki)
# ---------------------------------------------------------------------
from element_adapters import MITC8v12, SimoT6v18, MITC8D

# ---------------------------------------------------------------------
# ELEMENT REGISTRIES
# ---------------------------------------------------------------------
Q8_ELEMENTS = {
    "SimoQ8": Simo1993_Q8_ShellElement_v1p8,
    "SimoQ8v2": Simo1993_Q8_ShellElement_v2_MacNealPatched,
    "SimoQ8v3": Simo1993_Q8_ShellElement_v3,
    "ANS8_BDG6": ANS8_BDG6_ShellElement,
    "HBQ8": HBQ8_ShellElement,
    "MITC8v12": MITC8v12,          # <-- PAKAI ADAPTER (Fix beta_drill)
    "MITC8D": MITC8D,              # <-- PAKAI ADAPTER (Fix k_local & stress)
    "MNEALQ8": MacNealQ8_1992_native,
    "MNEALQ8D": MacNealQ8_1992_native_v3_drill,
    # "QUAD8": QUAD8,              # Buka jika ingin diuji
    # "QUAD8SS": QUAD8SS,          # Buka jika ingin diuji
}

T6_ELEMENTS = {
    "SIMOT6": Simo1993_Tri6_ShellElement_v2,
    "MITC6": MITC6_Tri_v4,
    "MH6T": MacNeal_MH6T_Tri_v3,
    "REZAIEE": Rezaiee2017_Tri6_ShellElement_v3,
}

# ==================================================================
# Problem 2-006 properties (SAP2000 Example 2-006 / Scordelis-Lo roof)
# ==================================================================
E = 432_000_000.0     # lb/ft^2
NU = 0.0
THICK = 0.25           # ft
RADIUS = 25.0          # ft
ANGLE_DEG = 40.0        # degrees (quarter-model arc)
SPAN = 25.0             # ft (half the full span, quarter model)
SURFACE_LOAD = -90.0    # lb/ft^2 (uniform self-weight, acts in -Z)

N_ARC = 6
N_LONG = 6

REF_UZ = -0.3086        # ft, vertical displacement at the free-edge center


def nid(i, j):
    return j * (N_ARC + 1) + i + 1


def node_xyz(i, j):
    """Continuous in i,j (uses division, not array lookup) -- fractional
    indices give the TRUE point on the curved surface, not a chord
    approximation."""
    theta = np.radians(ANGLE_DEG) * (i / N_ARC)
    x = RADIUS * np.sin(theta)
    y = SPAN * (j / N_LONG)
    z = RADIUS * np.cos(theta)
    return x, y, z


# ---- 8-node serendipity Q8 shape functions ----
def _q8_shape_and_deriv(xi, eta):
    N = np.array([
        -0.25 * (1 - xi) * (1 - eta) * (1 + xi + eta),
        -0.25 * (1 + xi) * (1 - eta) * (1 - xi + eta),
        -0.25 * (1 + xi) * (1 + eta) * (1 - xi - eta),
        -0.25 * (1 - xi) * (1 + eta) * (1 + xi - eta),
        0.5 * (1 - xi ** 2) * (1 - eta),
        0.5 * (1 + xi) * (1 - eta ** 2),
        0.5 * (1 - xi ** 2) * (1 + eta),
        0.5 * (1 - xi) * (1 - eta ** 2),
    ])
    dN_dxi = np.array([
        -0.25 * (-(1 - eta) * (1 + xi + eta) + (1 - xi) * (1 - eta)),
        -0.25 * ((1 - eta) * (1 - xi + eta) - (1 + xi) * (1 - eta)),
        -0.25 * ((1 + eta) * (1 - xi - eta) - (1 + xi) * (1 + eta)),
        -0.25 * (-(1 + eta) * (1 + xi - eta) + (1 - xi) * (1 + eta)),
        0.5 * (-2 * xi) * (1 - eta),
        0.5 * (1 - eta ** 2),
        0.5 * (-2 * xi) * (1 + eta),
        -0.5 * (1 - eta ** 2),
    ])
    dN_deta = np.array([
        -0.25 * (-(1 - xi) * (1 + xi + eta) + (1 - xi) * (1 - eta)),
        -0.25 * (-(1 + xi) * (1 - xi + eta) + (1 + xi) * (1 - eta)),
        -0.25 * ((1 + xi) * (1 - xi - eta) - (1 + xi) * (1 + eta)),
        -0.25 * ((1 - xi) * (1 + xi - eta) - (1 - xi) * (1 + eta)),
        -0.5 * (1 - xi ** 2),
        0.5 * (1 + xi) * (-2 * eta),
        0.5 * (1 - xi ** 2),
        0.5 * (1 - xi) * (-2 * eta),
    ])
    return N, dN_dxi, dN_deta


def consistent_gravity_load_q8(node_xyz_8, direction_vector):
    """Consistent nodal force (8,3) via 3x3 Gauss quadrature."""
    gp = [-np.sqrt(3.0 / 5.0), 0.0, np.sqrt(3.0 / 5.0)]
    gw = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]
    coords = np.asarray(node_xyz_8, dtype=float)
    fvec = np.asarray(direction_vector, dtype=float)
    f_nodal = np.zeros((8, 3))
    for xi, wxi in zip(gp, gw):
        for eta, weta in zip(gp, gw):
            N, dN_dxi, dN_deta = _q8_shape_and_deriv(xi, eta)
            g1 = dN_dxi @ coords
            g2 = dN_deta @ coords
            dA = np.linalg.norm(np.cross(g1, g2)) * wxi * weta
            f_nodal += np.outer(N, fvec) * dA
    return f_nodal


# ---- 6-node triangle shape functions ----
def _t6_shape_and_deriv(r, s):
    L1 = 1.0 - r - s
    N = np.array([
        L1 * (2 * L1 - 1), r * (2 * r - 1), s * (2 * s - 1),
        4 * L1 * r, 4 * r * s, 4 * s * L1,
    ])
    dN_dr = np.array([4*r + 4*s - 3, 4*r - 1, 0.0, 4 - 8*r - 4*s, 4*s, -4*s])
    dN_ds = np.array([4*r + 4*s - 3, 0.0, 4*s - 1, -4*r, 4*r, 4 - 4*r - 8*s])
    return N, dN_dr, dN_ds


def consistent_gravity_load_t6(node_xyz_6, direction_vector):
    """Consistent nodal force (6,3) via 6-point Gauss quadrature."""
    gp = [
        (0.445948490144588, 0.445948490144588),
        (0.10810301816807, 0.445948490144588),
        (0.445948490144588, 0.10810301816807),
        (0.091576213509771, 0.091576213509771),
        (0.816847572980459, 0.091576213509771),
        (0.091576213509771, 0.816847572980459),
    ]
    gw = [0.11169079483905, 0.11169079483905, 0.11169079483905,
          0.054975871827661, 0.054975871827661, 0.054975871827661]
    coords = np.asarray(node_xyz_6, dtype=float)
    fvec = np.asarray(direction_vector, dtype=float)
    f_nodal = np.zeros((6, 3))
    for (r, s), w in zip(gp, gw):
        N, dN_dr, dN_ds = _t6_shape_and_deriv(r, s)
        g1 = dN_dr @ coords
        g2 = dN_ds @ coords
        dA = np.linalg.norm(np.cross(g1, g2)) * w * 0.5  # ref-triangle area = 1/2
        f_nodal += np.outer(N, fvec) * dA
    return f_nodal


def apply_consistent_load_generic(model, elems_node_ids, magnitude, shape_fn,
                                   direction=(0.0, 0.0, 1.0)):
    fvec = np.array(direction, dtype=float) * magnitude
    nodal_totals = {}
    for ids in elems_node_ids:
        coords = np.array([[model.nodes[n].x, model.nodes[n].y, model.nodes[n].z]
                            for n in ids])
        f_elem = shape_fn(coords, fvec)
        for k, n in enumerate(ids):
            nodal_totals.setdefault(n, np.zeros(3))
            nodal_totals[n] += f_elem[k]
    for n, f in nodal_totals.items():
        for dof, val in zip((0, 1, 2), f):
            if abs(val) > 1e-14:
                model.add_load(n, dof, float(val))
    return nodal_totals


def build_grid_nodes(model, with_diagonal=False):
    nodes = {}
    for j in range(N_LONG + 1):
        for i in range(N_ARC + 1):
            nodes[nid(i, j)] = node_xyz(i, j)
    nid_counter = (N_ARC + 1) * (N_LONG + 1) + 1
    hmid, vmid, dmid = {}, {}, {}
    # 'h' midside nodes vary i (curved arc direction) -- true geometry.
    for j in range(N_LONG + 1):
        for i in range(N_ARC):
            x, y, z = node_xyz(i + 0.5, j)
            nodes[nid_counter] = (x, y, z)
            hmid[(i, j)] = nid_counter; nid_counter += 1
    # 'v' midside nodes vary j (straight span/Y direction) -- exact
    # under chord-averaging (no curvature to lose in this direction).
    for j in range(N_LONG):
        for i in range(N_ARC + 1):
            x1, y1, z1 = node_xyz(i, j); x2, y2, z2 = node_xyz(i, j + 1)
            nodes[nid_counter] = ((x1+x2)/2, (y1+y2)/2, (z1+z2)/2)
            vmid[(i, j)] = nid_counter; nid_counter += 1
    if with_diagonal:
        for j in range(N_LONG):
            for i in range(N_ARC):
                x, y, z = node_xyz(i + 0.5, j + 0.5)
                nodes[nid_counter] = (x, y, z)
                dmid[(i, j)] = nid_counter; nid_counter += 1
    for n_id, xyz in nodes.items():
        model.add_node(n_id, *xyz)
    return hmid, vmid, dmid


def apply_bcs(model, hmid, vmid):
    for i in range(N_ARC + 1):
        n_ = nid(i, 0)
        model.add_bc(n_, 0, 0.0); model.add_bc(n_, 2, 0.0)
        if i < N_ARC:
            nh = hmid[(i, 0)]
            model.add_bc(nh, 0, 0.0); model.add_bc(nh, 2, 0.0)
    for j in range(N_LONG + 1):
        n_ = nid(0, j)
        model.add_bc(n_, 0, 0.0); model.add_bc(n_, 4, 0.0)
        if j < N_LONG:
            nv = vmid[(0, j)]
            model.add_bc(nv, 0, 0.0); model.add_bc(nv, 4, 0.0)
    for i in range(N_ARC + 1):
        n_ = nid(i, N_LONG)
        model.add_bc(n_, 1, 0.0); model.add_bc(n_, 3, 0.0)
        if i < N_ARC:
            nh = hmid[(i, N_LONG)]
            model.add_bc(nh, 1, 0.0); model.add_bc(nh, 3, 0.0)


def build_q8_model(elem_cls):
    """
    Diubah agar menerima parameter `elem_cls` secara dinamis.
    """
    model = Model(ndof_per_node=6, autospc=False)
    hmid, vmid, _ = build_grid_nodes(model, with_diagonal=False)
    elems = []
    for j in range(N_LONG):
        for i in range(N_ARC):
            ids = [nid(i, j), nid(i+1, j), nid(i+1, j+1), nid(i, j+1),
                   hmid[(i, j)], vmid[(i+1, j)], hmid[(i, j+1)], vmid[(i, j)]]
            elems.append(ids)
            nodes = [model.nodes[n] for n in ids]
            model.add_element(elem_cls(len(model.elements)+1, nodes, E, NU, THICK))
    for node in model.nodes.values():
        model.add_bc(node.nid, 5, 0.0)
    apply_bcs(model, hmid, vmid)
    apply_consistent_load_generic(model, elems, SURFACE_LOAD, consistent_gravity_load_q8)
    return model


def build_t6_model(elem_cls):
    model = Model(ndof_per_node=6, autospc=False)
    hmid, vmid, dmid = build_grid_nodes(model, with_diagonal=True)
    elems = []
    for j in range(N_LONG):
        for i in range(N_ARC):
            sw, se = nid(i, j), nid(i+1, j)
            ne, nw = nid(i+1, j+1), nid(i, j+1)
            m_sw_se = hmid[(i, j)]
            m_se_ne = vmid[(i+1, j)]
            m_ne_nw = hmid[(i, j+1)]
            m_nw_sw = vmid[(i, j)]
            m_diag = dmid[(i, j)]

            ids1 = [sw, se, ne, m_sw_se, m_se_ne, m_diag]
            elems.append(ids1)
            n1 = [model.nodes[n] for n in ids1]
            model.add_element(elem_cls(len(model.elements)+1, n1, E, NU, THICK))

            ids2 = [sw, ne, nw, m_diag, m_ne_nw, m_nw_sw]
            elems.append(ids2)
            n2 = [model.nodes[n] for n in ids2]
            model.add_element(elem_cls(len(model.elements)+1, n2, E, NU, THICK))
    for node in model.nodes.values():
        model.add_bc(node.nid, 5, 0.0)
    apply_bcs(model, hmid, vmid)
    apply_consistent_load_generic(model, elems, SURFACE_LOAD, consistent_gravity_load_t6)
    return model


def main():
    node49 = nid(N_ARC, N_LONG)  # free corner (theta=ANGLE_DEG, y=SPAN)
    results = {}

    print("\n" + "=" * 78)
    print("PROBLEM 2-006 SCORDELIS-LO ROOF")
    print("=" * 78)
    print(f"{'Element':14s} {'Uz(node49)':>18} {'err %':>12}")

    # Iterasi seluruh Q8 Elements
    for label, elem_cls in Q8_ELEMENTS.items():
        try:
            model = build_q8_model(elem_cls)
            model.build()
            model.solve_static()
            u49 = model.get_displacement(node49)[2]
            results[label] = u49
            err = 100 * (u49 - REF_UZ) / REF_UZ
            print(f"{label:14s} {u49:18.6f} {err:+11.2f}%")
        except Exception as e:
            print(f"{label:14s} {'ERROR':>18} {'ERROR':>12} ({e})")

    # Iterasi seluruh T6 Elements
    for label, elem_cls in T6_ELEMENTS.items():
        try:
            model = build_t6_model(elem_cls)
            model.build()
            model.solve_static()
            u49 = model.get_displacement(node49)[2]
            results[label] = u49
            err = 100 * (u49 - REF_UZ) / REF_UZ
            print(f"{label:14s} {u49:18.6f} {err:+11.2f}%")
        except Exception as e:
            print(f"{label:14s} {'ERROR':>18} {'ERROR':>12} ({e})")

    return results


if __name__ == "__main__":
    main()
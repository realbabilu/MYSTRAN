#!/usr/bin/env python3
"""
Problem 2-002 — Q8/T6 benchmark with MSC/Nastran reference runs.

Two-stage workflow:

1) FIRST RUN
   - Run the Python Q8/T6 elements.
   - Create MSC/Nastran decks in ./working_nastran/
   - Do NOT launch Nastran automatically.
   - Plot Python results and any Nastran F06 files that already exist.

2) SOLVE THE DECKS
   Run the generated .dat files with MSC/Nastran.

3) SECOND RUN
   - Run this script again.
   - Existing .F06/.f06 files in ./working_nastran/ are detected.
   - Their tip responses are parsed and plotted together with Python results.

Nastran elements:
   CQUAD8  = one quadratic quad per parent rectangle
   CTRIA6  = two quadratic triangles per parent rectangle

The CTRIA6 diagonal midside GRID is shared by the two triangles.

This benchmark uses the independent beam-theory reference values, NOT
SAP2000 values, as the primary reference.
"""

from pathlib import Path
import csv
import sys
import numpy as np
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
LINEAR_DIR = HERE.parent / "linear"
if str(LINEAR_DIR) not in sys.path:
    sys.path.insert(0, str(LINEAR_DIR))

from core import Model
from Simo1993_Tri6_ShellElement_v2 import Simo1993_Tri6_ShellElement_v2
from MacNeal_MH6T_Tri_v3 import MacNeal_MH6T_Tri_v3
from MITC6_Tri_v4 import MITC6_Tri_v4
from Rezaiee2017_Tri6_v3 import Rezaiee2017_Tri6_ShellElement_v3

# Q8 Python elements
from ANS8_BDG6_ShellElement import ANS8_BDG6_ShellElement
from HBQ8_ShellElement import HBQ8_ShellElement
from MITC8_ShellElement_v1p2 import MITC8_ShellElement_v1p2
from MITC8D_ShellElement import MITC8D_ShellElement

from MacNealQ8_1992_native import MacNealQ8_1992_native
from MacNealQ8_1992_native_v3_drill import MacNealQ8_1992_native_v3_drill

from KikuchiMacNeal_Q8_ShellElement_v1 import KikuchiMacNeal_Q8_ShellElement_v1 
from KikuchiMacNeal_ANS8_v1 import KikuchiMacNeal_ANS8_v1
from KikuchiMacNeal_HBQ8_v1 import KikuchiMacNeal_HBQ8_v1

from Simo1993_Q8_ShellElement_v1p8_standalone import (
    Simo1993_Q8_ShellElement_v1p8,
)


# Dedicated T6 implementation
from element_adapters import MITC8v12, SimoT6v18

# ---------------------------------------------------------------------
# Problem 2-002
# ---------------------------------------------------------------------
LENGTH = 6.0
HEIGHT = 0.2
THICK = 0.1
E = 1.0e7
NU = 0.3

# Independent beam-theory values used by the benchmark.
REF = {
    1: ("UX", 3.000e-5),
    2: ("UZ", 1.081e-1),
    3: ("UY", 4.321e-1),
    4: ("UY", 3.410e-3),
    5: ("UX", 9.000e-4),
    6: ("RZ", 3.600e-2),
}
REF_DOF = {1: 0, 2: 2, 3: 1, 4: 1, 5: 0, 6: 5}

# For LC4/LC5 the two tip values have opposite signs.
USE_ABS = {4: True, 5: True}

# Several meshes are useful, but the fine end can be reduced if desired.
NX_LIST = [1, 2, 4, 8, 16, 24, 32, 48, 64]

NASTRAN_DIR = HERE / "working_nastran"
NASTRAN_DIR.mkdir(parents=True, exist_ok=True)
MYSTRAN_DIR = HERE / "working_mystran"
MYSTRAN_DIR.mkdir(parents=True, exist_ok=True)

PLOT_ALL_PATH = "prob_2_002_q8_t6_nastran.png"
PLOT_MYSTRAN_VS_NASTRAN_PATH = "prob_2_002_q8_t6_mystran_Vs_nastran.png"
PLOT_PYTHON_VS_NASTRAN_PATH = "prob_2_002_q8_t6_python_Vs_nastran.png"
CSV_PATH = "prob_2_002_q8_t6_nastran.csv"

Q8_ELEMENTS = {


    "SimoQ8": Simo1993_Q8_ShellElement_v1p8,
    "SimoQ8_Kikuchi": KikuchiMacNeal_Q8_ShellElement_v1 ,

    "ANS8_BDG6": ANS8_BDG6_ShellElement,
    "ANS8_KIKU":KikuchiMacNeal_ANS8_v1,
    "HBQ8": HBQ8_ShellElement,
    "HBQ8_KIKU": KikuchiMacNeal_HBQ8_v1,

    "MITC8": MITC8_ShellElement_v1p2,
    "MITC8v12": MITC8v12,
    "MITC8D": MITC8D_ShellElement,
    "MNEALQ8": MacNealQ8_1992_native,
    "MNEALQ8D": MacNealQ8_1992_native_v3_drill,



}

T6_ELEMENTS = {
    "SIMOT6": Simo1993_Tri6_ShellElement_v2,
    "MITC6": MITC6_Tri_v4,
    "MH6T": MacNeal_MH6T_Tri_v3,
    "REZAIEE": Rezaiee2017_Tri6_ShellElement_v3,
}

MYSTRAN_FORMULATIONS = [
    ("cquad8", "QUAD8TYP", "SIMOQ8", "MYSTRAN CQUAD8 SIMOQ8", "cquad8_simoq8"),
    ("cquad8", "QUAD8TYP", "ANS8BDG6", "MYSTRAN CQUAD8 ANS8_BDG6", "cquad8_ans8bdg6"),
    ("cquad8", "QUAD8TYP", "HBQ8", "MYSTRAN CQUAD8 HBQ8", "cquad8_hbq8"),
    ("cquad8", "QUAD8TYP", "MITC8", "MYSTRAN CQUAD8 MITC8", "cquad8_mitc8"),
    ("cquad8", "QUAD8TYP", "MITC8D", "MYSTRAN CQUAD8 MITC8D", "cquad8_mitc8d"),
    ("cquad8", "QUAD8TYP", "MACQ8D", "MYSTRAN CQUAD8 MACQ8D", "cquad8_macq8d"),
    ("ctria6", "TRIA6TYP", "SIMOT6", "MYSTRAN CTRIA6 SIMOT6", "ctria6_simot6"),
    ("ctria6", "TRIA6TYP", "MITC6", "MYSTRAN CTRIA6 MITC6", "ctria6_mitc6"),
    ("ctria6", "TRIA6TYP", "MH6T", "MYSTRAN CTRIA6 MH6T", "ctria6_mh6t"),
    ("ctria6", "TRIA6TYP", "REZAIEE", "MYSTRAN CTRIA6 REZAIEE", "ctria6_rezaiee"),
]


# ---------------------------------------------------------------------
# Common mesh
# ---------------------------------------------------------------------
def build_parent_grid(nx, diagonal=False):
    """
    Common 6 x 0.2 parent strip.

    Corner nodes:
        bottom j=0, top j=1

    Q8 midsides:
        hmid = longitudinal edge midsides
        vmid = vertical edge midsides

    T6 only:
        dmid = shared midside of SW-NE diagonal.
    """
    model = Model(ndof_per_node=6, autospc=False)

    dx = LENGTH / nx
    dz = HEIGHT

    corner = {}
    hmid = {}
    vmid = {}
    dmid = {}

    nid = 1

    for j in range(2):
        for i in range(nx + 1):
            model.add_node(nid, i * dx, 0.0, j * dz)
            corner[(i, j)] = nid
            nid += 1

    for j in range(2):
        for i in range(nx):
            model.add_node(nid, (i + 0.5) * dx, 0.0, j * dz)
            hmid[(i, j)] = nid
            nid += 1

    for i in range(nx + 1):
        model.add_node(nid, i * dx, 0.0, 0.5 * dz)
        vmid[i] = nid
        nid += 1

    if diagonal:
        for i in range(nx):
            model.add_node(nid, (i + 0.5) * dx, 0.0, 0.5 * dz)
            dmid[i] = nid
            nid += 1

    return model, corner, hmid, vmid, dmid


def build_q8_old(nx, shell_class):
    model, corner, hmid, vmid, _ = build_parent_grid(nx, False)

    eid = 1
    for i in range(nx):
        nodes = [
            model.nodes[corner[(i, 0)]],
            model.nodes[corner[(i + 1, 0)]],
            model.nodes[corner[(i + 1, 1)]],
            model.nodes[corner[(i, 1)]],
            model.nodes[hmid[(i, 0)]],
            model.nodes[vmid[i + 1]],
            model.nodes[hmid[(i, 1)]],
            model.nodes[vmid[i]],
        ]
        model.add_element(shell_class(eid, nodes, E, NU, THICK))
        eid += 1

    return (
        model,
        corner[(0, 0)],
        corner[(0, 1)],
        corner[(nx, 0)],
        corner[(nx, 1)],
    )

def build_q8(nx, shell_class):
    # Diubah dari autospc=False menjadi autospc=True
    model = Model(ndof_per_node=6, autospc=True)

    dx = LENGTH / nx
    dz = HEIGHT

    corner = {}
    hmid = {}
    vmid = {}
    dmid = {}

    nid = 1

    for j in range(2):
        for i in range(nx + 1):
            model.add_node(nid, i * dx, 0.0, j * dz)
            corner[(i, j)] = nid
            nid += 1

    for j in range(2):
        for i in range(nx):
            model.add_node(nid, (i + 0.5) * dx, 0.0, j * dz)
            hmid[(i, j)] = nid
            nid += 1

    for i in range(nx + 1):
        model.add_node(nid, i * dx, 0.0, 0.5 * dz)
        vmid[i] = nid
        nid += 1

    eid = 1
    for i in range(nx):
        nodes = [
            model.nodes[corner[(i, 0)]],
            model.nodes[corner[(i + 1, 0)]],
            model.nodes[corner[(i + 1, 1)]],
            model.nodes[corner[(i, 1)]],
            model.nodes[hmid[(i, 0)]],
            model.nodes[vmid[i + 1]],
            model.nodes[hmid[(i, 1)]],
            model.nodes[vmid[i]],
        ]
        model.add_element(shell_class(eid, nodes, E, NU, THICK))
        eid += 1

    return (
        model,
        corner[(0, 0)],
        corner[(0, 1)],
        corner[(nx, 0)],
        corner[(nx, 1)],
    )


def build_t6(nx, shell_class):
    """
    Two T6 elements per parent rectangle.

    T6-1: SW-SE-NE
    T6-2: SW-NE-NW

    The SW-NE diagonal midside GRID is shared.
    """
    model, corner, hmid, vmid, dmid = build_parent_grid(nx, True)

    eid = 1
    for i in range(nx):
        sw = corner[(i, 0)]
        se = corner[(i + 1, 0)]
        ne = corner[(i + 1, 1)]
        nw = corner[(i, 1)]

        m_sw_se = hmid[(i, 0)]
        m_se_ne = vmid[i + 1]
        m_sw_ne = dmid[i]
        m_ne_nw = hmid[(i, 1)]
        m_nw_sw = vmid[i]

        nodes1 = [
            model.nodes[sw], model.nodes[se], model.nodes[ne],
            model.nodes[m_sw_se], model.nodes[m_se_ne],
            model.nodes[m_sw_ne],
        ]
        model.add_element(shell_class(eid, nodes1, E, NU, THICK))
        eid += 1

        nodes2 = [
            model.nodes[sw], model.nodes[ne], model.nodes[nw],
            model.nodes[m_sw_ne], model.nodes[m_ne_nw],
            model.nodes[m_nw_sw],
        ]
        model.add_element(shell_class(eid, nodes2, E, NU, THICK))
        eid += 1

    return (
        model,
        corner[(0, 0)],
        corner[(0, 1)],
        corner[(nx, 0)],
        corner[(nx, 1)],
    )


# ---------------------------------------------------------------------
# Loads / BCs
# ---------------------------------------------------------------------
def apply_load_case(model, fixed_bottom, fixed_top, tip_bottom, tip_top, lc):
    # Exact support pattern used in the benchmark:
    # bottom: UX UY UZ RZ
    # top:    UX UY RZ
    for dof in [0, 1, 2, 5]:
        model.add_bc(fixed_bottom, dof, 0.0)
    for dof in [0, 1, 5]:
        model.add_bc(fixed_top, dof, 0.0)

    if lc == 1:
        model.add_load(tip_bottom, 0, +0.5)
        model.add_load(tip_top,    0, +0.5)

    elif lc == 2:
        model.add_load(fixed_top,  2, -0.5)
        model.add_load(tip_bottom, 2, +0.5)
        model.add_load(tip_top,    2, +0.5)

    elif lc == 3:
        model.add_load(tip_bottom, 1, +0.5)
        model.add_load(tip_top,    1, +0.5)

    elif lc == 4:
        model.add_load(tip_bottom, 1, -5.0)
        model.add_load(tip_top,    1, +5.0)

    elif lc == 5:
        model.add_load(tip_bottom, 0, -5.0)
        model.add_load(tip_top,    0, +5.0)

    elif lc == 6:
        model.add_load(tip_bottom, 5, +0.5)
        model.add_load(tip_top,    5, +0.5)

    else:
        raise ValueError(f"Unknown load case {lc}")


def solve_one(nx, kind, shell_class, lc):
    if kind == "Q8":
        model, fb, ft, tb, tt = build_q8(nx, shell_class)
    elif kind == "T6":
        model, fb, ft, tb, tt = build_t6(nx, shell_class)
    else:
        raise ValueError(kind)

    apply_load_case(model, fb, ft, tb, tt, lc)
    model.build()
    model.solve_static()

    dof = REF_DOF[lc]
    ub = model.get_displacement(tb)[dof]
    ut = model.get_displacement(tt)[dof]

    if USE_ABS.get(lc, False):
        value = 0.5 * (abs(ub) + abs(ut))
    else:
        value = 0.5 * (ub + ut)

    return float(value)


# ---------------------------------------------------------------------
# Python convergence
# ---------------------------------------------------------------------
def run_python():
    all_results = {}
    failures = []

    for label, cls in Q8_ELEMENTS.items():
        print(f"\nPython Q8: {label}")
        vals = {lc: [] for lc in range(1, 7)}

        for nx in NX_LIST:
            for lc in range(1, 7):
                try:
                    val = solve_one(nx, "Q8", cls, lc)
                except Exception as exc:
                    val = np.nan
                    failures.append((label, "Q8", nx, lc, repr(exc)))
                    print(f"  FAIL nx={nx} LC{lc}: {exc}")
                vals[lc].append(val)

        all_results[label] = vals

    for label, cls in T6_ELEMENTS.items():
        print(f"\nPython T6: {label}")
        vals = {lc: [] for lc in range(1, 7)}

        for nx in NX_LIST:
            for lc in range(1, 7):
                try:
                    val = solve_one(nx, "T6", cls, lc)
                except Exception as exc:
                    val = np.nan
                    failures.append((label, "T6", nx, lc, repr(exc)))
                    print(f"  FAIL nx={nx} LC{lc}: {exc}")
                vals[lc].append(val)

        all_results[label] = vals

    return all_results, failures


# ---------------------------------------------------------------------
# Nastran deck generation
# ---------------------------------------------------------------------
def _nid(nx, i, j, typ="c"):
    n_corner = 2 * (nx + 1)
    n_hmid = 2 * nx

    if typ == "c":
        return j * (nx + 1) + i + 1

    if typ == "h":
        return n_corner + j * nx + i + 1

    if typ == "v":
        return n_corner + n_hmid + i + 1

    # Diagonal midside nodes, one per parent quad.
    if typ == "d":
        return n_corner + n_hmid + (nx + 1) + i + 1

    raise ValueError(typ)


def _grid_lines(nx):
    dx = LENGTH / nx
    lines = []

    for j in range(2):
        for i in range(nx + 1):
            nid = _nid(nx, i, j, "c")
            lines.append(
                f"GRID,{nid},,{i*dx:.8E},0.0,{j*HEIGHT:.8E}"
            )

    for j in range(2):
        for i in range(nx):
            nid = _nid(nx, i, j, "h")
            lines.append(
                f"GRID,{nid},,{(i+0.5)*dx:.8E},0.0,{j*HEIGHT:.8E}"
            )

    for i in range(nx + 1):
        nid = _nid(nx, i, 0, "v")
        lines.append(
            f"GRID,{nid},,{i*dx:.8E},0.0,{0.5*HEIGHT:.8E}"
        )

    return lines


def _quad_elements(nx):
    lines = []
    eid = 1

    for i in range(nx):
        n1 = _nid(nx, i,   0, "c")
        n2 = _nid(nx, i+1, 0, "c")
        n3 = _nid(nx, i+1, 1, "c")
        n4 = _nid(nx, i,   1, "c")
        n5 = _nid(nx, i,   0, "h")
        n6 = _nid(nx, i+1, 0, "v")
        n7 = _nid(nx, i,   1, "h")
        n8 = _nid(nx, i,   0, "v")

        lines.append(f"CQUAD8,{eid},1,{n1},{n2},{n3},{n4},{n5},{n6}")
        lines.append(f"+,{n7},{n8}")
        eid += 1

    return lines


def _tri_elements(nx):
    lines = []
    eid = 1

    for i in range(nx):
        sw = _nid(nx, i,   0, "c")
        se = _nid(nx, i+1, 0, "c")
        ne = _nid(nx, i+1, 1, "c")
        nw = _nid(nx, i,   1, "c")

        m_sw_se = _nid(nx, i,   0, "h")
        m_se_ne = _nid(nx, i+1, 0, "v")
        m_sw_ne = _nid(nx, i,   0, "d")
        m_ne_nw = _nid(nx, i,   1, "h")
        m_nw_sw = _nid(nx, i,   0, "v")

        # SW-SE-NE
        lines.append(
            f"CTRIA6,{eid},1,{sw},{se},{ne},"
            f"{m_sw_se},{m_se_ne},{m_sw_ne}"
        )
        eid += 1

        # SW-NE-NW
        lines.append(
            f"CTRIA6,{eid},1,{sw},{ne},{nw},"
            f"{m_sw_ne},{m_ne_nw},{m_nw_sw}"
        )
        eid += 1

    return lines


def _load_lines(nx):
    bottom = _nid(nx, 0, 0, "c")
    top = _nid(nx, 0, 1, "c")
    tip_bottom = _nid(nx, nx, 0, "c")
    tip_top = _nid(nx, nx, 1, "c")

    lines = [
        # LC1 axial
        f"FORCE,1,{tip_bottom},,5.0E-1,1.,0.,0.",
        f"FORCE,1,{tip_top},,5.0E-1,1.,0.,0.",

        # LC2 in-plane shear + bending
        f"FORCE,2,{top},,-5.0E-1,0.,0.,1.",
        f"FORCE,2,{tip_bottom},,5.0E-1,0.,0.,1.",
        f"FORCE,2,{tip_top},,5.0E-1,0.,0.,1.",

        # LC3 out-of-plane shear + bending
        f"FORCE,3,{tip_bottom},,5.0E-1,0.,1.,0.",
        f"FORCE,3,{tip_top},,5.0E-1,0.,1.,0.",

        # LC4 twist
        f"FORCE,4,{tip_bottom},,-5.0E0,0.,1.,0.",
        f"FORCE,4,{tip_top},,5.0E0,0.,1.,0.",

        # LC5 in-plane moment
        f"FORCE,5,{tip_bottom},,-5.0E0,1.,0.,0.",
        f"FORCE,5,{tip_top},,5.0E0,1.,0.,0.",

        # LC6 out-of-plane moment
        f"MOMENT,6,{tip_bottom},,5.0E-1,0.,0.,1.",
        f"MOMENT,6,{tip_top},,5.0E-1,0.,0.,1.",
    ]
    return lines


def write_solver_deck(nx, element_type, out_dir, solver="NASTRAN",
                      param_name=None, selector=None, file_tag=None):
    """
    element_type = 'cquad8' or 'ctria6'
    """
    if element_type not in ("cquad8", "ctria6"):
        raise ValueError(element_type)

    suffix = file_tag or ("cquad8" if element_type == "cquad8" else "ctria6")
    path = out_dir / f"prob_2_002_nx{nx:02d}_{suffix}.dat"

    lines = [
        "SOL 101",
        "CEND",
        f"TITLE = PROBLEM 2-002 {solver} {suffix.upper()} NX={nx}",
        "ECHO = NONE",
    ]

    for lc in range(1, 7):
        lines += [
            f"SUBCASE {lc}",
            f"  LABEL = LC{lc}",
            "  SPC = 1",
            f"  LOAD = {lc}",
            "  DISPLACEMENT = ALL",
        ]

    lines += [
        "BEGIN BULK",
        f"MAT1,1,{E:.8E},,{NU:.8E}",
        f"PSHELL,1,1,{THICK:.8E},1",
        "PARAM,POST,-1",
    ]
    if solver.upper() == "MYSTRAN":
        lines.append("PARAM,COUPMASS,0")
    if param_name and selector:
        lines.append(f"PARAM,{param_name},{selector}")

    lines += _grid_lines(nx)

    if element_type == "cquad8":
        lines += _quad_elements(nx)
    else:
        # Add diagonal midside GRID cards after the common grids.
        dx = LENGTH / nx
        for i in range(nx):
            nid = _nid(nx, i, 0, "d")
            lines.append(
                f"GRID,{nid},,{(i+0.5)*dx:.8E},0.0,{0.5*HEIGHT:.8E}"
            )
        lines += _tri_elements(nx)

    bottom = _nid(nx, 0, 0, "c")
    top = _nid(nx, 0, 1, "c")

    # bottom = UX UY UZ RZ
    # top    = UX UY RZ
    lines += [
        f"SPC1,1,1236,{bottom}",
        f"SPC1,1,126,{top}",
    ]

    lines += _load_lines(nx)
    lines += ["ENDDATA", ""]

    path.write_text("\n".join(lines), encoding="ascii")
    return path


def write_nastran_deck(nx, element_type):
    return write_solver_deck(nx, element_type, NASTRAN_DIR, solver="NASTRAN")


def write_mystran_deck(nx, element_type, param_name, selector, file_tag):
    return write_solver_deck(
        nx,
        element_type,
        MYSTRAN_DIR,
        solver="MYSTRAN",
        param_name=param_name,
        selector=selector,
        file_tag=file_tag,
    )


def write_all_nastran_decks():
    print(f"\nWriting Nastran decks to: {NASTRAN_DIR}")
    for nx in NX_LIST:
        for etype in ("cquad8", "ctria6"):
            path = write_nastran_deck(nx, etype)
            print(f"  {path.name}")


def write_all_mystran_decks():
    print(f"\nWriting MYSTRAN decks to: {MYSTRAN_DIR}")
    for nx in NX_LIST:
        for etype, param_name, selector, _label, file_tag in MYSTRAN_FORMULATIONS:
            path = write_mystran_deck(nx, etype, param_name, selector, file_tag)
            print(f"  {path.name}")


# ---------------------------------------------------------------------
# Nastran F06 parsing
# ---------------------------------------------------------------------
def _find_f06(deck):
    for suffix in (".f06", ".F06"):
        p = deck.with_suffix(suffix)
        if p.exists():
            return p
    return None


def _parse_one_f06_series(deck, nx, parse_f06_displacements):
    f06 = _find_f06(deck)
    if f06 is None:
        return None

    try:
        subcases = parse_f06_displacements(str(f06))
    except Exception as exc:
        print(f"  Could not parse {f06.name}: {exc}")
        return None

    tip_bottom = _nid(nx, nx, 0, "c")
    tip_top = _nid(nx, nx, 1, "c")

    vals = {}
    for lc in range(1, 7):
        nodes = subcases.get(lc, {})

        if tip_bottom not in nodes or tip_top not in nodes:
            vals[lc] = np.nan
            print(
                f"  {f06.name}: LC{lc} missing tip node "
                f"{tip_bottom}/{tip_top}"
            )
            continue

        dof = REF_DOF[lc]

        try:
            ub = float(nodes[tip_bottom][dof])
            ut = float(nodes[tip_top][dof])
        except (KeyError, IndexError, TypeError, ValueError):
            vals[lc] = np.nan
            print(
                f"  {f06.name}: LC{lc} missing DOF {dof} "
                f"at tip nodes {tip_bottom}/{tip_top}"
            )
            continue

        if USE_ABS.get(lc, False):
            vals[lc] = 0.5 * (abs(ub) + abs(ut))
        else:
            vals[lc] = 0.5 * (ub + ut)

    print(f"  Read {f06.name}")
    return vals


def parse_nastran_results():
    """
    Read existing Nastran F06 files using the established
    solver_f06_convergence parser.

    Important:
      - LC1/2/3/6 use signed average of the two tip nodes.
      - LC4/5 use average of absolute tip values.
      - Missing DOFs are reported as NaN, never silently converted to zero.
    """
    try:
        from solver_f06_convergence import parse_f06_displacements
    except ImportError:
        print("solver_f06_convergence.py not found; Nastran results skipped.")
        return {}, {}

    cquad8 = {}
    ctria6 = {}

    for nx in NX_LIST:
        for etype, target in (("cquad8", cquad8), ("ctria6", ctria6)):
            deck = NASTRAN_DIR / f"prob_2_002_nx{nx:02d}_{etype}.dat"
            vals = _parse_one_f06_series(deck, nx, parse_f06_displacements)
            if vals is not None:
                target[nx] = vals

    return cquad8, ctria6


def parse_mystran_results():
    try:
        from solver_f06_convergence import parse_f06_displacements
    except ImportError:
        print("solver_f06_convergence.py not found; MYSTRAN results skipped.")
        return {}

    results = {}
    for nx in NX_LIST:
        for _etype, _param_name, _selector, label, file_tag in MYSTRAN_FORMULATIONS:
            deck = MYSTRAN_DIR / f"prob_2_002_nx{nx:02d}_{file_tag}.dat"
            vals = _parse_one_f06_series(deck, nx, parse_f06_displacements)
            if vals is not None:
                results.setdefault(label, {})[nx] = vals

    return results


# ---------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------
def load_python_results_from_csv(path=CSV_PATH):
    path = HERE / path
    if not path.exists():
        return None

    results = {}
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("source") != "Python":
                continue
            label = row["element"]
            nx = int(row["nx"])
            if nx not in NX_LIST:
                continue
            data = results.setdefault(label, {lc: [np.nan] * len(NX_LIST) for lc in range(1, 7)})
            idx = NX_LIST.index(nx)
            for lc, col in [
                (1, "LC1_UX"),
                (2, "LC2_UZ"),
                (3, "LC3_UY"),
                (4, "LC4_UY_abs"),
                (5, "LC5_UX_abs"),
                (6, "LC6_RZ"),
            ]:
                data[lc][idx] = float(row[col])

    if not results:
        return None
    return results


def write_csv(python_results, external_results, path=CSV_PATH):

    path = HERE / path

    labels = list(python_results.keys())
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "source", "element", "nx",
            "LC1_UX", "LC2_UZ", "LC3_UY",
            "LC4_UY_abs", "LC5_UX_abs", "LC6_RZ",
        ])

        for label in labels:
            data = python_results[label]
            for k, nx in enumerate(NX_LIST):
                w.writerow([
                    "Python", label, nx,
                    *[data[lc][k] for lc in range(1, 7)]
                ])

        for label, series in external_results.items():
            source = "NASTRAN" if label.startswith("NASTRAN") else "MYSTRAN"
            for nx in sorted(series):
                w.writerow([
                    source, label, nx,
                    *[series[nx][lc] for lc in range(1, 7)]
                ])

    print(f"CSV written: {path}")


# ---------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------
def plot_results(python_results, external_results=None, include_python=True,
                 include_external=True, path=PLOT_ALL_PATH, title=None):

    path = HERE / path
    external_results = external_results or {}

    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    axes = axes.ravel()

    titles = [
        "LC1 — Axial extension (UX)",
        "LC2 — In-plane shear + bending (UZ)",
        "LC3 — Out-of-plane shear + bending (UY)",
        "LC4 — Twist (|UY|)",
        "LC5 — In-plane moment (|UX|)",
        "LC6 — Out-of-plane moment (RZ)",
    ]

    py_markers = ["o", "s", "^", "D", "v", "p", "*", "h", "P", "X"]
    external_markers = ["X", "P", "D", "s", "^", "v", "<", ">", "h", "8", "p", "*"]

    for lc, ax in enumerate(axes, start=1):
        if include_python:
            for idx, (label, data) in enumerate(python_results.items()):
                ax.plot(
                    NX_LIST, data[lc],
                    marker=py_markers[idx % len(py_markers)],
                    linestyle="-",
                    linewidth=1.4,
                    markersize=4.5,
                    label=label,
                )

        if include_external:
            for idx, (label, series) in enumerate(external_results.items()):
                if not series:
                    continue
                xs = sorted(series)
                ys = [series[x][lc] for x in xs]
                linestyle = "--" if label.startswith("NASTRAN") else "-."
                linewidth = 2.2 if label.startswith("NASTRAN") else 1.6
                ax.plot(
                    xs, ys,
                    marker=external_markers[idx % len(external_markers)],
                    linestyle=linestyle,
                    linewidth=linewidth,
                    markersize=6.0,
                    label=label,
                )

        ax.axhline(
            REF[lc][1],
            linestyle=":",
            linewidth=1.5,
            label="Independent",
        )

        ax.set_title(titles[lc - 1])
        ax.set_xlabel("Parent mesh Nx")
        ax.set_ylabel(REF[lc][0])
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=7, loc="best")

    fig.suptitle(
        title or "Problem 2-002 — Q8/T6 Python, MYSTRAN, and NASTRAN",
        fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)

    print(f"Plot written: {path}")


# ---------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------
def print_summary(python_results, external_results, failures):
    nx = NX_LIST[-1]

    print("\n" + "=" * 110)
    print(f"Problem 2-002 — finest Python mesh Nx={nx}")
    print("=" * 110)

    print("\nPython elements:")
    for label, data in python_results.items():
        vals = [data[lc][-1] for lc in range(1, 7)]
        print(
            f"{label:<16}"
            + "".join(f"  {v: .6e}" for v in vals)
        )

    if external_results:
        print("\nExternal solver elements:")
        for label, series in external_results.items():
            if nx in series:
                print(
                    f"{label:<28}"
                    + "".join(f"  {series[nx][lc]: .6e}" for lc in range(1, 7))
                )

    print("\nIndependent reference:")
    print(
        "REFERENCE        "
        + "".join(f"  {REF[lc][1]: .6e}" for lc in range(1, 7))
    )

    if failures:
        print("\nPython failures:")
        for failure in failures:
            print(" ", failure)

    n_nastran_q8 = len(external_results.get("NASTRAN CQUAD8", {}))
    n_nastran_t6 = len(external_results.get("NASTRAN CTRIA6", {}))
    n_mystran = sum(
        len(series) for label, series in external_results.items()
        if label.startswith("MYSTRAN")
    )
    print(
        f"\nExisting Nastran results found: "
        f"CQUAD8={n_nastran_q8}/{len(NX_LIST)}, "
        f"CTRIA6={n_nastran_t6}/{len(NX_LIST)}"
    )
    print(f"Existing MYSTRAN result rows found: {n_mystran}")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    print("Problem 2-002 Q8/T6 Nastran benchmark")
    print(f"Working directory: {HERE}")
    print(f"Nastran decks/results: {NASTRAN_DIR}")
    print(f"MYSTRAN decks/results: {MYSTRAN_DIR}")

    # Every run writes deterministic decks. It does not launch external solvers.
    # After external solving, rerunning this script automatically sees the F06 files.
    write_all_nastran_decks()
    write_all_mystran_decks()

    print("\nLoading Python element results...")
    python_results = load_python_results_from_csv()
    if python_results is None:
        print("  No Python CSV cache found; running Python elements.")
        python_results, failures = run_python()
    else:
        failures = []
        print(f"  Loaded Python cache from {HERE / CSV_PATH}")

    print("\nChecking existing Nastran F06 files...")
    nastran_q8, nastran_t6 = parse_nastran_results()
    external_results = {
        "NASTRAN CQUAD8": nastran_q8,
        "NASTRAN CTRIA6": nastran_t6,
    }

    print("\nChecking existing MYSTRAN F06 files...")
    external_results.update(parse_mystran_results())

    write_csv(python_results, external_results)
    plot_results(
        python_results,
        external_results,
        include_python=True,
        include_external=True,
        path=PLOT_ALL_PATH,
        title="Problem 2-002 — Q8/T6 Python + MYSTRAN + NASTRAN",
    )
    plot_results(
        {},
        external_results,
        include_python=False,
        include_external=True,
        path=PLOT_MYSTRAN_VS_NASTRAN_PATH,
        title="Problem 2-002 — MYSTRAN CQUAD8/CTRIA6 vs NASTRAN",
    )
    plot_results(
        python_results,
        {
            "NASTRAN CQUAD8": nastran_q8,
            "NASTRAN CTRIA6": nastran_t6,
        },
        include_python=True,
        include_external=True,
        path=PLOT_PYTHON_VS_NASTRAN_PATH,
        title="Problem 2-002 — Python Q8/T6 vs NASTRAN",
    )
    print_summary(
        python_results, external_results, failures
    )

    print("\nDone.")
    print("Next step: solve any missing .dat files in working_nastran/working_mystran.")
    print("Then run this Python script again to refresh the combined CSV and plots.")


if __name__ == "__main__":
    main()

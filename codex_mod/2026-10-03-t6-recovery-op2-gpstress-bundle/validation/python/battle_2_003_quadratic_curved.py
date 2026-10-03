#!/usr/bin/env python3
"""
battle_2_003_quadratic_curved.py
Problem 2-003 - curved beam quadratic benchmark with Python, MYSTRAN,
and NASTRAN deck generation.

Geometry is a quarter annulus in the XY plane:
  ri = 4.12, ro = 4.32, angle = 90 deg, thickness = 0.1.

Load cases:
  LC1 = in-plane shear (UY)
  LC2 = out-of-plane response (UZ)

Support method is aligned to the reference decks: clamp the two nodes on
the fixed radial edge with SPC 1234, not all 6 DOF.
"""

from pathlib import Path
import csv
import sys

import matplotlib.pyplot as plt
import numpy as np

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

from ANS8_BDG6_ShellElement import ANS8_BDG6_ShellElement
from HBQ8_ShellElement import HBQ8_ShellElement
from MITC8_ShellElement_v1p2 import MITC8_ShellElement_v1p2
from MITC8D_ShellElement import MITC8D_ShellElement
from MacNealQ8_1992_native_v2 import MacNealQ8_1992_native_v2
from MacNealQ8_1992_native_v3_drill import MacNealQ8_1992_native_v3_drill
from KikuchiMacNeal_Q8_ShellElement_v1 import KikuchiMacNeal_Q8_ShellElement_v1 
from KikuchiMacNeal_ANS8_v1 import KikuchiMacNeal_ANS8_v1
from KikuchiMacNeal_HBQ8_v1 import KikuchiMacNeal_HBQ8_v1

from Simo1993_Q8_ShellElement_v1p8_standalone import (
    Simo1993_Q8_ShellElement_v1p8,
)

from element_adapters import MITC8v12, MITC8D, SimoT6v18

RI = 4.12
RO = 4.32
THICK = 0.1
E = 1.0e7
NU = 0.25
ANGLE = np.pi / 2.0

REF = {
    1: ("UY", 8.8654708e-02),
    2: ("UZ", 5.0040000e-01),
}

REF_DOF = {1: 1, 2: 2}
NX_LIST = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64]

NASTRAN_DIR = HERE / "working_nastran"
MYSTRAN_DIR = HERE / "working_mystran"
NASTRAN_DIR.mkdir(parents=True, exist_ok=True)
MYSTRAN_DIR.mkdir(parents=True, exist_ok=True)

PLOT_ALL_PATH = "prob_2_003_q8_t6_nastran.png"
PLOT_MYSTRAN_VS_NASTRAN_PATH = "prob_2_003_q8_t6_mystran_Vs_nastran.png"
PLOT_PYTHON_VS_NASTRAN_PATH = "prob_2_003_q8_t6_python_Vs_nastran.png"
CSV_PATH = "prob_2_003_q8_t6_nastran.csv"

Q8_ELEMENTS = {
    "SimoQ8": Simo1993_Q8_ShellElement_v1p8,
    "SimoQ8_Kikuchi": KikuchiMacNeal_Q8_ShellElement_v1 ,
    "ANS8_BDG6": ANS8_BDG6_ShellElement,
    "ANS8_KIKU":KikuchiMacNeal_ANS8_v1,
    "HBQ8": HBQ8_ShellElement,
    "HBQ8_KIKU": KikuchiMacNeal_HBQ8_v1,

    "MITC8": MITC8_ShellElement_v1p2,
    "MITC8v12": MITC8v12,
    "MITC8D": MITC8D,
    "MNEALQ8": MacNealQ8_1992_native_v2,
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


def build_nodes(n_arc, model, with_diagonal=False):
    theta = np.linspace(0.0, ANGLE, n_arc + 1)
    r_vals = np.array([RI, RO])
    node_map = {}
    dmid = {}
    nid = 1

    for i in range(n_arc + 1):
        for j in range(2):
            r = r_vals[j]
            th = theta[i]
            model.add_node(nid, r * np.cos(th), r * np.sin(th), 0.0)
            node_map[(i, j, "c")] = nid
            nid += 1

    for i in range(n_arc):
        th_mid = 0.5 * (theta[i] + theta[i + 1])
        for j in range(2):
            r = r_vals[j]
            model.add_node(nid, r * np.cos(th_mid), r * np.sin(th_mid), 0.0)
            node_map[(i, j, "h")] = nid
            nid += 1

    for i in range(n_arc + 1):
        r_mid = 0.5 * (RI + RO)
        th = theta[i]
        model.add_node(nid, r_mid * np.cos(th), r_mid * np.sin(th), 0.0)
        node_map[(i, 0, "v")] = nid
        nid += 1

    if with_diagonal:
        th_mids = 0.5 * (theta[:-1] + theta[1:])
        for i in range(n_arc):
            r_mid = 0.5 * (RI + RO)
            th_mid = th_mids[i]
            model.add_node(nid, r_mid * np.cos(th_mid), r_mid * np.sin(th_mid), 0.0)
            dmid[(i, 0)] = nid
            nid += 1

    return node_map, dmid


def build_q8(n_arc, elem_cls):
    model = Model(ndof_per_node=6, autospc=False)
    node_map, _ = build_nodes(n_arc, model, with_diagonal=False)

    for i in range(n_arc):
        ids = [
            node_map[(i, 0, "c")], node_map[(i + 1, 0, "c")],
            node_map[(i + 1, 1, "c")], node_map[(i, 1, "c")],
            node_map[(i, 0, "h")], node_map[(i + 1, 0, "v")],
            node_map[(i, 1, "h")], node_map[(i, 0, "v")],
        ]
        model.add_element(elem_cls(len(model.elements) + 1, [model.nodes[n] for n in ids], E, NU, THICK))

    return (
        model,
        node_map[(0, 0, "c")],
        node_map[(0, 1, "c")],
        node_map[(n_arc, 0, "c")],
        node_map[(n_arc, 1, "c")],
    )


def build_t6(n_arc, elem_cls, alternate_diag=True):
    model = Model(ndof_per_node=6, autospc=False)
    node_map, dmid = build_nodes(n_arc, model, with_diagonal=True)

    for i in range(n_arc):
        sw = node_map[(i, 0, "c")]
        se = node_map[(i + 1, 0, "c")]
        ne = node_map[(i + 1, 1, "c")]
        nw = node_map[(i, 1, "c")]
        m_sw_se = node_map[(i, 0, "h")]
        m_se_ne = node_map[(i + 1, 0, "v")]
        m_ne_nw = node_map[(i, 1, "h")]
        m_nw_sw = node_map[(i, 0, "v")]
        m_diag = dmid[(i, 0)]

        if (not alternate_diag) or (i % 2 == 0):
            ids1 = [sw, se, ne, m_sw_se, m_se_ne, m_diag]
            ids2 = [sw, ne, nw, m_diag, m_ne_nw, m_nw_sw]
        else:
            ids1 = [sw, se, nw, m_sw_se, m_diag, m_nw_sw]
            ids2 = [se, ne, nw, m_se_ne, m_ne_nw, m_diag]

        model.add_element(elem_cls(len(model.elements) + 1, [model.nodes[n] for n in ids1], E, NU, THICK))
        model.add_element(elem_cls(len(model.elements) + 1, [model.nodes[n] for n in ids2], E, NU, THICK))

    return (
        model,
        node_map[(0, 0, "c")],
        node_map[(0, 1, "c")],
        node_map[(n_arc, 0, "c")],
        node_map[(n_arc, 1, "c")],
    )


def apply_bc_load(model, fixed_inner, fixed_outer, tip_inner, tip_outer, lc):
    for nid in (fixed_inner, fixed_outer):
        for dof in (0, 1, 2, 3):
            model.add_bc(nid, dof, 0.0)
    model.autospc = True
    dof = 1 if lc == 1 else 2
    model.add_load(tip_inner, dof, 0.5)
    model.add_load(tip_outer, dof, 0.5)


def solve_one(label, nx, lc, alternate_diag=True):
    if label in Q8_ELEMENTS:
        model, fi, fo, ti, to = build_q8(nx, Q8_ELEMENTS[label])
    elif label in T6_ELEMENTS:
        model, fi, fo, ti, to = build_t6(nx, T6_ELEMENTS[label], alternate_diag)
    else:
        raise ValueError(label)

    apply_bc_load(model, fi, fo, ti, to, lc)
    model.build()
    model.solve_static()
    dof = REF_DOF[lc]
    ub = model.get_displacement(ti)[dof]
    ut = model.get_displacement(to)[dof]
    return float(0.5 * (ub + ut))


def run_python(alternate_diag=True):
    results = {}
    failures = []

    for label, _cls in Q8_ELEMENTS.items():
        vals = {1: [], 2: []}
        for nx in NX_LIST:
            for lc in (1, 2):
                try:
                    vals[lc].append(solve_one(label, nx, lc, alternate_diag))
                except Exception as exc:
                    vals[lc].append(np.nan)
                    failures.append((label, nx, lc, repr(exc)))
        results[label] = vals

    for label, _cls in T6_ELEMENTS.items():
        vals = {1: [], 2: []}
        for nx in NX_LIST:
            for lc in (1, 2):
                try:
                    vals[lc].append(solve_one(label, nx, lc, alternate_diag))
                except Exception as exc:
                    vals[lc].append(np.nan)
                    failures.append((label, nx, lc, repr(exc)))
        results[label] = vals

    return results, failures


def print_mesh_method_check():
    print("\nMesh/method check at Nx=24")
    for label in T6_ELEMENTS:
        try:
            fixed = solve_one(label, 24, 2, alternate_diag=False)
            alt = solve_one(label, 24, 2, alternate_diag=True)
            print(f"  {label:<10s} LC2 fixed={fixed: .6e} alt={alt: .6e}")
        except Exception as exc:
            print(f"  {label:<10s} ERROR ({exc})")


def _nid(nx, i, j, typ="c"):
    n_corner = 2 * (nx + 1)
    n_hmid = 2 * nx
    n_vmid = nx + 1
    if typ == "c":
        return 1 + 2 * i + j
    if typ == "h":
        return n_corner + 2 * i + j + 1
    if typ == "v":
        return n_corner + n_hmid + i + 1
    if typ == "d":
        return n_corner + n_hmid + n_vmid + i + 1
    raise ValueError(typ)


def _grid_lines(nx, diagonal=False):
    theta = np.linspace(0.0, ANGLE, nx + 1)
    lines = []

    for i in range(nx + 1):
        for j, r in enumerate((RI, RO)):
            nid = _nid(nx, i, j, "c")
            lines.append(f"GRID,{nid},,{r*np.cos(theta[i]):.8E},{r*np.sin(theta[i]):.8E},0.0")

    for i in range(nx):
        th = 0.5 * (theta[i] + theta[i + 1])
        for j, r in enumerate((RI, RO)):
            nid = _nid(nx, i, j, "h")
            lines.append(f"GRID,{nid},,{r*np.cos(th):.8E},{r*np.sin(th):.8E},0.0")

    r_mid = 0.5 * (RI + RO)
    for i in range(nx + 1):
        nid = _nid(nx, i, 0, "v")
        lines.append(f"GRID,{nid},,{r_mid*np.cos(theta[i]):.8E},{r_mid*np.sin(theta[i]):.8E},0.0")

    if diagonal:
        for i in range(nx):
            th = 0.5 * (theta[i] + theta[i + 1])
            nid = _nid(nx, i, 0, "d")
            lines.append(f"GRID,{nid},,{r_mid*np.cos(th):.8E},{r_mid*np.sin(th):.8E},0.0")

    return lines


def _quad_elements(nx):
    lines = []
    eid = 1
    for i in range(nx):
        n1 = _nid(nx, i, 0, "c")
        n2 = _nid(nx, i + 1, 0, "c")
        n3 = _nid(nx, i + 1, 1, "c")
        n4 = _nid(nx, i, 1, "c")
        n5 = _nid(nx, i, 0, "h")
        n6 = _nid(nx, i + 1, 0, "v")
        n7 = _nid(nx, i, 1, "h")
        n8 = _nid(nx, i, 0, "v")
        lines.append(f"CQUAD8,{eid},1,{n1},{n2},{n3},{n4},{n5},{n6}")
        lines.append(f"+,{n7},{n8}")
        eid += 1
    return lines


def _tri_elements(nx):
    lines = []
    eid = 1
    for i in range(nx):
        sw = _nid(nx, i, 0, "c")
        se = _nid(nx, i + 1, 0, "c")
        ne = _nid(nx, i + 1, 1, "c")
        nw = _nid(nx, i, 1, "c")
        m_sw_se = _nid(nx, i, 0, "h")
        m_se_ne = _nid(nx, i + 1, 0, "v")
        m_sw_ne = _nid(nx, i, 0, "d")
        m_ne_nw = _nid(nx, i, 1, "h")
        m_nw_sw = _nid(nx, i, 0, "v")
        lines.append(f"CTRIA6,{eid},1,{sw},{se},{ne},{m_sw_se},{m_se_ne},{m_sw_ne}")
        eid += 1
        lines.append(f"CTRIA6,{eid},1,{sw},{ne},{nw},{m_sw_ne},{m_ne_nw},{m_nw_sw}")
        eid += 1
    return lines


def _load_lines(nx):
    tip_inner = _nid(nx, nx, 0, "c")
    tip_outer = _nid(nx, nx, 1, "c")
    top = _nid(nx, 0, 1, "c")
    return [
        f"FORCE,1,{tip_inner},,5.0E-1,0.,1.,0.",
        f"FORCE,1,{tip_outer},,5.0E-1,0.,1.,0.",
        f"FORCE,2,{top},,-5.0E-1,0.,0.,1.",
        f"FORCE,2,{tip_inner},,5.0E-1,0.,0.,1.",
        f"FORCE,2,{tip_outer},,5.0E-1,0.,0.,1.",
    ]


def write_solver_deck(nx, element_type, out_dir, solver="NASTRAN", param_name=None, selector=None, file_tag=None):
    suffix = file_tag or element_type
    path = out_dir / f"prob_2_003_nx{nx:02d}_{suffix}.dat"
    lines = [
        "SOL 101",
        "CEND",
        f"TITLE = PROBLEM 2-003 {solver} {suffix.upper()} NX={nx}",
        "ECHO = NONE",
        "SUBCASE 1",
        "  LABEL = LC1",
        "  SPC = 1",
        "  LOAD = 1",
        "  DISPLACEMENT = ALL",
        "SUBCASE 2",
        "  LABEL = LC2",
        "  SPC = 1",
        "  LOAD = 2",
        "  DISPLACEMENT = ALL",
        "BEGIN BULK",
        f"MAT1,1,{E:.8E},,{NU:.8E}",
        f"PSHELL,1,1,{THICK:.8E},1",
        "PARAM,POST,-1",
    ]
    if solver.upper() == "MYSTRAN":
        lines.append("PARAM,COUPMASS,0")
    if param_name and selector:
        lines.append(f"PARAM,{param_name},{selector}")
    lines += _grid_lines(nx, diagonal=(element_type == "ctria6"))
    if element_type == "cquad8":
        lines += _quad_elements(nx)
    else:
        lines += _tri_elements(nx)
    lines += [
        f"SPC1,1,1234,{_nid(nx, 0, 0, 'c')}",
        f"SPC1,1,1234,{_nid(nx, 0, 1, 'c')}",
    ]
    lines += _load_lines(nx)
    lines += ["ENDDATA", ""]
    path.write_text("\n".join(lines), encoding="ascii")
    return path


def write_nastran_deck(nx, element_type):
    return write_solver_deck(nx, element_type, NASTRAN_DIR, solver="NASTRAN")


def write_mystran_deck(nx, element_type, param_name, selector, file_tag):
    return write_solver_deck(nx, element_type, MYSTRAN_DIR, solver="MYSTRAN", param_name=param_name, selector=selector, file_tag=file_tag)


def write_all_nastran_decks():
    print(f"\nWriting Nastran decks to: {NASTRAN_DIR}")
    for nx in NX_LIST:
        for etype in ("cquad8", "ctria6"):
            print(f"  {write_nastran_deck(nx, etype).name}")


def write_all_mystran_decks():
    print(f"\nWriting MYSTRAN decks to: {MYSTRAN_DIR}")
    for nx in NX_LIST:
        for etype, param_name, selector, _label, file_tag in MYSTRAN_FORMULATIONS:
            print(f"  {write_mystran_deck(nx, etype, param_name, selector, file_tag).name}")


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

    tip_inner = _nid(nx, nx, 0, "c")
    tip_outer = _nid(nx, nx, 1, "c")
    vals = {}
    for lc in (1, 2):
        nodes = subcases.get(lc, {})
        if tip_inner not in nodes or tip_outer not in nodes:
            vals[lc] = np.nan
            continue
        dof = REF_DOF[lc]
        try:
            ub = float(nodes[tip_inner][dof])
            ut = float(nodes[tip_outer][dof])
        except Exception:
            vals[lc] = np.nan
            continue
        vals[lc] = 0.5 * (ub + ut)
    print(f"  Read {f06.name}")
    return vals


def parse_nastran_results():
    try:
        from solver_f06_convergence import parse_f06_displacements
    except ImportError:
        print("solver_f06_convergence.py not found; Nastran results skipped.")
        return {}, {}

    cquad8 = {}
    ctria6 = {}
    for nx in NX_LIST:
        for etype, target in (("cquad8", cquad8), ("ctria6", ctria6)):
            deck = NASTRAN_DIR / f"prob_2_003_nx{nx:02d}_{etype}.dat"
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
            deck = MYSTRAN_DIR / f"prob_2_003_nx{nx:02d}_{file_tag}.dat"
            vals = _parse_one_f06_series(deck, nx, parse_f06_displacements)
            if vals is not None:
                results.setdefault(label, {})[nx] = vals
    return results


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
            data = results.setdefault(label, {lc: [np.nan] * len(NX_LIST) for lc in (1, 2)})
            idx = NX_LIST.index(nx)
            data[1][idx] = float(row["LC1_UY"])
            data[2][idx] = float(row["LC2_UZ"])
    return results or None


def write_csv(python_results, external_results, path=CSV_PATH):
    path = HERE / path
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["source", "element", "nx", "LC1_UY", "LC2_UZ"])
        for label, data in python_results.items():
            for k, nx in enumerate(NX_LIST):
                w.writerow(["Python", label, nx, data[1][k], data[2][k]])
        for label, series in external_results.items():
            source = "NASTRAN" if label.startswith("NASTRAN") else "MYSTRAN"
            for nx in sorted(series):
                w.writerow([source, label, nx, series[nx][1], series[nx][2]])
    print(f"CSV written: {path}")


def plot_results(python_results, external_results=None, include_python=True, include_external=True, path=PLOT_ALL_PATH, title=None):
    path = HERE / path
    external_results = external_results or {}
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2))
    titles = ["LC1 - In-plane shear (UY)", "LC2 - Out-of-plane response (UZ)"]
    py_markers = ["o", "s", "^", "D", "v", "p", "*", "h", "P", "X"]
    ext_markers = ["X", "P", "D", "s", "^", "v", "<", ">", "h", "8", "p", "*"]

    for idx, lc in enumerate((1, 2)):
        ax = axes[idx]
        if include_python:
            for j, (label, data) in enumerate(python_results.items()):
                ax.plot(NX_LIST, data[lc], marker=py_markers[j % len(py_markers)], linewidth=1.4, markersize=4.5, label=label)
        if include_external:
            for j, (label, series) in enumerate(external_results.items()):
                if not series:
                    continue
                xs = sorted(series)
                ys = [series[x][lc] for x in xs]
                ax.plot(xs, ys, marker=ext_markers[j % len(ext_markers)], linestyle="--" if label.startswith("NASTRAN") else "-.", linewidth=1.8, markersize=5.5, label=label)
        ax.axhline(REF[lc][1], linestyle=":", linewidth=1.4, label="Reference")
        ax.set_title(titles[idx])
        ax.set_xlabel("Parent mesh Nx")
        ax.set_ylabel(REF[lc][0])
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=7, loc="best")

    fig.suptitle(title or "Problem 2-003 - Q8/T6 Python, MYSTRAN, and NASTRAN", fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)
    print(f"Plot written: {path}")


def print_summary(python_results, external_results, failures):
    nx = NX_LIST[-1]
    print("\n" + "=" * 110)
    print(f"Problem 2-003 - finest Python mesh Nx={nx}")
    print("=" * 110)
    print("\nPython elements:")
    for label, data in python_results.items():
        print(f"{label:<16}" + "".join(f"  {data[lc][-1]: .6e}" for lc in (1, 2)))
    if external_results:
        print("\nExternal solver elements:")
        for label, series in external_results.items():
            if nx in series:
                print(f"{label:<28}" + "".join(f"  {series[nx][lc]: .6e}" for lc in (1, 2)))
    print("\nReference:")
    print("REFERENCE        " + "".join(f"  {REF[lc][1]: .6e}" for lc in (1, 2)))
    if failures:
        print("\nPython failures:")
        for failure in failures:
            print(" ", failure)


def main():
    print("Problem 2-003 quadratic curved benchmark")
    print(f"Working directory: {HERE}")
    print(f"Nastran decks/results: {NASTRAN_DIR}")
    print(f"MYSTRAN decks/results: {MYSTRAN_DIR}")

    write_all_nastran_decks()
    write_all_mystran_decks()

    print("\nLoading Python element results...")
    python_results = load_python_results_from_csv()
    if python_results is None:
        print("  No Python CSV cache found; running Python elements.")
        python_results, failures = run_python(alternate_diag=True)
    else:
        failures = []
        print(f"  Loaded Python cache from {HERE / CSV_PATH}")

    print_mesh_method_check()

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
        title="Problem 2-003 - Q8/T6 Python + MYSTRAN + NASTRAN",
    )
    plot_results(
        {},
        external_results,
        include_python=False,
        include_external=True,
        path=PLOT_MYSTRAN_VS_NASTRAN_PATH,
        title="Problem 2-003 - MYSTRAN CQUAD8/CTRIA6 vs NASTRAN",
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
        title="Problem 2-003 - Python Q8/T6 vs NASTRAN",
    )
    print_summary(python_results, external_results, failures)
    print("\nDone.")
    print("Next step: solve any missing .dat files in working_nastran/working_mystran.")
    print("Then rerun this script to refresh the combined CSV and plots.")


if __name__ == "__main__":
    main()

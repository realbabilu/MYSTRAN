from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\mystran2")
MITC3P_DIR = ROOT / "mitc3+"
CLAUDEAI_DIR = ROOT / "pythonfem" / "claudeai"
if str(MITC3P_DIR) not in sys.path:
    sys.path.insert(0, str(MITC3P_DIR))
if str(CLAUDEAI_DIR) not in sys.path:
    sys.path.insert(0, str(CLAUDEAI_DIR))

from core import Node  # type: ignore
import mitc3plus as mitc3plus_mod  # type: ignore


FLOAT_RE = re.compile(r"[-+]?\d+\.\d+E[-+]\d+")


def parse_err_named_matrix(path: Path, name: str, nrow: int, ncol: int) -> np.ndarray:
    lines = path.read_text().splitlines()
    begin = f"MITC3P_MATRIX_BEGIN {name}"
    end = f"MITC3P_MATRIX_END {name}"
    start = None
    for i, line in enumerate(lines):
        if begin in " ".join(line.split()):
            start = i + 1
            break
    if start is None:
        raise RuntimeError(f"Could not find {begin} in {path}")
    rows: list[list[float]] = []
    for i in range(start, len(lines)):
        norm = " ".join(lines[i].split())
        if end in norm:
            break
        if norm.startswith("ROW"):
            vals = [float(x) for x in FLOAT_RE.findall(lines[i])]
            rows.append(vals)
    if len(rows) != nrow:
        raise RuntimeError(f"Expected {nrow} rows in {name}, got {len(rows)}")
    return np.array(rows, dtype=float)


def rel_diff(a: np.ndarray, b: np.ndarray) -> float:
    denom = np.linalg.norm(b)
    if denom == 0.0:
        return float("inf")
    return np.linalg.norm(a - b) / denom


def build_python_component_full(nodes: list[Node], e: float, nu: float, h: float, drilling_penalty: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    _, _, _, x_local = mitc3plus_mod._local_frame(nodes)
    xy = x_local[:, :2]
    j = mitc3plus_mod._dh().T @ xy
    det_j = np.linalg.det(j)
    j_inv = np.linalg.inv(j)
    cb = mitc3plus_mod._C_bending(e, nu, h)
    cs = mitc3plus_mod._C_shear(e, nu, h)
    bs_factory = mitc3plus_mod._Bs_mitc3plus(xy, h)
    r_gp, s_gp, w_gp = mitc3plus_mod._gauss_triangle_7pt()

    idx_b = [3, 4, 9, 10, 15, 16, 18, 19]
    idx_s = [2, 3, 4, 8, 9, 10, 14, 15, 16, 18, 19]

    kb = np.zeros((20, 20))
    ks = np.zeros((20, 20))
    kd = np.zeros((20, 20))

    for i_gp in range(len(r_gp)):
        r = r_gp[i_gp]
        s = s_gp[i_gp]
        fac = w_gp[i_gp] * det_j

        bb = mitc3plus_mod._Bb_at(r, s, j_inv)
        kb_gp = bb.T @ cb @ bb * fac
        for ii, gi in enumerate(idx_b):
            for jj, gj in enumerate(idx_b):
                kb[gi, gj] += kb_gp[ii, jj]

        bs = bs_factory(r, s)
        cs_cov = j_inv.T @ cs @ j_inv
        ks_gp = bs.T @ cs_cov @ bs * fac
        for ii, gi in enumerate(idx_s):
            for jj, gj in enumerate(idx_s):
                ks[gi, gj] += ks_gp[ii, jj]

    for idx_rz in [5, 11, 17]:
        kd[idx_rz, idx_rz] += drilling_penalty
    return kb, ks, kd


def main() -> None:
    err_path = ROOT / "MYSTRANSolver-18.0.0" / "run_debug" / "dkmq" / "shell_static_mct" / "ctria3_mitc3p_one_elem_shell_debug.ERR"
    kaa_err = parse_err_named_matrix(err_path, "KAA", 18, 18)
    kab_err = parse_err_named_matrix(err_path, "KAB", 18, 2)
    kbb_err = parse_err_named_matrix(err_path, "KBB", 2, 2)
    kcond_err = parse_err_named_matrix(err_path, "KCOND", 18, 18)

    nodes = [
        Node(1, 0.0, 0.0, 0.0),
        Node(2, 0.5, 0.0, 0.0),
        Node(9, 0.5, 0.3, 0.0),
    ]
    for i, node in enumerate(nodes):
        node.dofs = list(range(6 * i, 6 * i + 6))

    elem = mitc3plus_mod.MITC3Plus(1, nodes, 1.07e7, 0.30, 0.1)
    kb20, ks20, kd20 = build_python_component_full(nodes, 1.07e7, 0.30, 0.1, elem.drilling_penalty)
    kfull = kb20 + ks20 + kd20

    def split(mat20: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        kaa = mat20[:18, :18]
        kab = mat20[:18, 18:]
        kba = mat20[18:, :18]
        kbb = mat20[18:, 18:]
        kcond = kaa - kab @ np.linalg.pinv(kbb) @ kba
        return kaa, kab, kbb, kcond

    kaa_b, kab_b, kbb_b, kcond_b = split(kb20)
    kaa_s, kab_s, kbb_s, kcond_s = split(ks20)
    kaa_d, kab_d, kbb_d, kcond_d = split(kd20)
    kaa_p, kab_p, kbb_p, kcond_p = split(kfull)

    print("MITC3+ plate block decomposition vs MYSTRAN ERR")
    print(f"  err_vs_python_plate_kaa   = {rel_diff(kaa_err, kaa_p):.9E}")
    print(f"  err_vs_python_plate_kab   = {rel_diff(kab_err, kab_p):.9E}")
    print(f"  err_vs_python_plate_kbb   = {rel_diff(kbb_err, kbb_p):.9E}")
    print(f"  err_vs_python_plate_kcond = {rel_diff(kcond_err, kcond_p):.9E}")
    print(f"  err_vs_bending_kaa        = {rel_diff(kaa_err, kaa_b):.9E}")
    print(f"  err_vs_shear_kaa          = {rel_diff(kaa_err, kaa_s):.9E}")
    print(f"  err_vs_drill_kaa          = {rel_diff(kaa_err, kaa_d):.9E}")
    print(f"  err_vs_bending_kab        = {rel_diff(kab_err, kab_b):.9E}")
    print(f"  err_vs_shear_kab          = {rel_diff(kab_err, kab_s):.9E}")
    print(f"  err_vs_bending_kbb        = {rel_diff(kbb_err, kbb_b):.9E}")
    print(f"  err_vs_shear_kbb          = {rel_diff(kbb_err, kbb_s):.9E}")
    print(f"  ||kaa_err||               = {np.linalg.norm(kaa_err):.9E}")
    print(f"  ||kaa_b||                 = {np.linalg.norm(kaa_b):.9E}")
    print(f"  ||kaa_s||                 = {np.linalg.norm(kaa_s):.9E}")
    print(f"  ||kab_err||               = {np.linalg.norm(kab_err):.9E}")
    print(f"  ||kab_b||                 = {np.linalg.norm(kab_b):.9E}")
    print(f"  ||kab_s||                 = {np.linalg.norm(kab_s):.9E}")
    print(f"  ||kbb_err||               = {np.linalg.norm(kbb_err):.9E}")
    print(f"  ||kbb_b||                 = {np.linalg.norm(kbb_b):.9E}")
    print(f"  ||kbb_s||                 = {np.linalg.norm(kbb_s):.9E}")


if __name__ == "__main__":
    main()

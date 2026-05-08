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
from mitc3plus import MITC3Plus  # type: ignore


FLOAT_RE = re.compile(r"[-+]?\d+\.\d+E[-+]\d+")


def parse_bug_full_ke(path: Path, ndof: int) -> np.ndarray:
    lines = path.read_text().splitlines()
    start = None
    for i, line in enumerate(lines):
        if "KE element stiffness matrix in local element coordinate system" in line:
            start = i
            break
    if start is None:
        raise RuntimeError(f"Could not find KE block in {path}")

    rows: list[list[float]] = []
    i = start + 1
    while i < len(lines) and len(rows) < ndof:
        line = lines[i]
        if line.strip().startswith("Row"):
            row_vals: list[float] = []
            i += 1
            while i < len(lines):
                nums = [float(x) for x in FLOAT_RE.findall(lines[i])]
                if not nums:
                    break
                row_vals.extend(nums)
                i += 1
                if len(row_vals) >= ndof:
                    break
            if len(row_vals) != ndof:
                raise RuntimeError(f"Expected {ndof} values in {path}, got {len(row_vals)}")
            rows.append(row_vals)
        else:
            i += 1
    if len(rows) != ndof:
        raise RuntimeError(f"Expected {ndof} rows in {path}, got {len(rows)}")
    return np.array(rows, dtype=float)


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
            if len(vals) != ncol:
                raise RuntimeError(f"Expected {ncol} values in {name} row, got {len(vals)}")
            rows.append(vals)
    if len(rows) != nrow:
        raise RuntimeError(f"Expected {nrow} rows in {name}, got {len(rows)}")
    return np.array(rows, dtype=float)


def rel_diff(a: np.ndarray, b: np.ndarray) -> float:
    denom = np.linalg.norm(b)
    if denom == 0.0:
        return float("inf")
    return np.linalg.norm(a - b) / denom


def tri_shell_blocks() -> tuple[list[int], list[int], list[int]]:
    mem_idx: list[int] = []
    plate_idx: list[int] = []
    rz_idx: list[int] = []
    for n in range(3):
        base = 6 * n
        mem_idx.extend([base + 0, base + 1])
        plate_idx.extend([base + 2, base + 3, base + 4])
        rz_idx.append(base + 5)
    return mem_idx, plate_idx, rz_idx


def build_python_plate_full(nodes: list[Node], e: float, nu: float, h: float, drilling_penalty: float) -> np.ndarray:
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
    k_full = np.zeros((20, 20))

    for i_gp in range(len(r_gp)):
        r = r_gp[i_gp]
        s = s_gp[i_gp]
        fac = w_gp[i_gp] * det_j

        bb = mitc3plus_mod._Bb_at(r, s, j_inv)
        kb_gp = bb.T @ cb @ bb * fac
        for ii, gi in enumerate(idx_b):
            for jj, gj in enumerate(idx_b):
                k_full[gi, gj] += kb_gp[ii, jj]

        bs = bs_factory(r, s)
        cs_cov = j_inv.T @ cs @ j_inv
        ks_gp = bs.T @ cs_cov @ bs * fac
        for ii, gi in enumerate(idx_s):
            for jj, gj in enumerate(idx_s):
                k_full[gi, gj] += ks_gp[ii, jj]

    for idx_rz in [5, 11, 17]:
        k_full[idx_rz, idx_rz] += drilling_penalty
    return k_full


def main() -> None:
    bug_path = (
        ROOT
        / "MYSTRANSolver-18.0.0"
        / "run_debug"
        / "dkmq"
        / "shell_static_mct"
        / "ctria3_mitc3p_one_elem_shell_debug.BUG"
    )
    err_path = bug_path.with_suffix(".ERR")
    legacy_bug_path = (
        ROOT
        / "MYSTRANSolver-18.0.0"
        / "run_debug"
        / "dkmq"
        / "shell_static_mct"
        / "ctria3_one_elem_shell_debug.BUG"
    )
    k_mystran = parse_bug_full_ke(bug_path, 18)
    k_legacy = parse_bug_full_ke(legacy_bug_path, 18)
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

    elem = MITC3Plus(1, nodes, 1.07e7, 0.30, 0.1)
    k_python = elem.k_local()
    k_full_python = mitc3plus_mod._build_K_full(nodes, 1.07e7, 0.30, 0.1, elem.drilling_penalty)
    k_full_plate_py = build_python_plate_full(nodes, 1.07e7, 0.30, 0.1, elem.drilling_penalty)
    kaa_py = k_full_python[:18, :18]
    kab_py = k_full_python[:18, 18:]
    kbb_py = k_full_python[18:, 18:]
    kcond_py = mitc3plus_mod._static_condense(k_full_python)
    kaa_plate_py = k_full_plate_py[:18, :18]
    kab_plate_py = k_full_plate_py[:18, 18:]
    kbb_plate_py = k_full_plate_py[18:, 18:]
    kcond_plate_py = mitc3plus_mod._static_condense(k_full_plate_py)
    mem_idx, plate_idx, rz_idx = tri_shell_blocks()

    print("MITC3+ one-element shell stiffness compare")
    print(f"  mystran_bug_path = {bug_path}")
    print(f"  mystran_err_path = {err_path}")
    print(f"  bug_vs_python_max_abs   = {np.max(np.abs(k_mystran - k_python)):.9E}")
    print(f"  bug_vs_python_rel_fro   = {rel_diff(k_mystran, k_python):.9E}")
    print(f"  bug_vs_python_mem_rel   = {rel_diff(k_mystran[np.ix_(mem_idx, mem_idx)], k_python[np.ix_(mem_idx, mem_idx)]):.9E}")
    print(f"  bug_vs_python_plate_rel = {rel_diff(k_mystran[np.ix_(plate_idx, plate_idx)], k_python[np.ix_(plate_idx, plate_idx)]):.9E}")
    print(f"  bug_vs_python_rz_rel    = {rel_diff(k_mystran[np.ix_(rz_idx, rz_idx)], k_python[np.ix_(rz_idx, rz_idx)]):.9E}")
    print(f"  bug_vs_legacy_full      = {rel_diff(k_mystran, k_legacy):.9E}")
    print(f"  bug_vs_legacy_plate     = {rel_diff(k_mystran[np.ix_(plate_idx, plate_idx)], k_legacy[np.ix_(plate_idx, plate_idx)]):.9E}")
    print(f"  err_kaa_rel_diff        = {rel_diff(kaa_err, kaa_py):.9E}")
    print(f"  err_kab_rel_diff        = {rel_diff(kab_err, kab_py):.9E}")
    print(f"  err_kbb_rel_diff        = {rel_diff(kbb_err, kbb_py):.9E}")
    print(f"  err_kcond_rel_diff      = {rel_diff(kcond_err, kcond_py):.9E}")
    print(f"  err_kcond_mem_rel       = {rel_diff(kcond_err[np.ix_(mem_idx, mem_idx)], kcond_py[np.ix_(mem_idx, mem_idx)]):.9E}")
    print(f"  err_kcond_plate_rel     = {rel_diff(kcond_err[np.ix_(plate_idx, plate_idx)], kcond_py[np.ix_(plate_idx, plate_idx)]):.9E}")
    print(f"  err_kcond_rz_rel        = {rel_diff(kcond_err[np.ix_(rz_idx, rz_idx)], kcond_py[np.ix_(rz_idx, rz_idx)]):.9E}")
    print(f"  err_kaa_plateonly_rel   = {rel_diff(kaa_err, kaa_plate_py):.9E}")
    print(f"  err_kab_plateonly_rel   = {rel_diff(kab_err, kab_plate_py):.9E}")
    print(f"  err_kbb_plateonly_rel   = {rel_diff(kbb_err, kbb_plate_py):.9E}")
    print(f"  err_kcond_plateonly_rel = {rel_diff(kcond_err, kcond_plate_py):.9E}")
    print(f"  err_kcond_plateonly_blk = {rel_diff(kcond_err[np.ix_(plate_idx, plate_idx)], kcond_plate_py[np.ix_(plate_idx, plate_idx)]):.9E}")
    print(f"  err_kcond_plateonly_rz  = {rel_diff(kcond_err[np.ix_(rz_idx, rz_idx)], kcond_plate_py[np.ix_(rz_idx, rz_idx)]):.9E}")
    print(f"  mystran_bug_sym_err     = {np.max(np.abs(k_mystran - k_mystran.T)):.9E}")
    print(f"  mystran_err_kcond_sym   = {np.max(np.abs(kcond_err - kcond_err.T)):.9E}")
    print(f"  python_sym_err   = {np.max(np.abs(k_python - k_python.T)):.9E}")


if __name__ == "__main__":
    main()

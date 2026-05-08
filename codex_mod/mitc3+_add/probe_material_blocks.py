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
            rows.append([float(x) for x in FLOAT_RE.findall(lines[i])])
    return np.array(rows, dtype=float)


def parse_scalar(path: Path, key: str) -> float:
    for line in path.read_text().splitlines():
        if key in line:
            vals = [float(x) for x in FLOAT_RE.findall(line)]
            if vals:
                return vals[0]
    raise RuntimeError(f"Could not find scalar {key}")


def rel_diff(a: np.ndarray, b: np.ndarray) -> float:
    denom = np.linalg.norm(b)
    if denom == 0.0:
        return float("inf")
    return np.linalg.norm(a - b) / denom


def main() -> None:
    err_path = ROOT / "MYSTRANSolver-18.0.0" / "run_debug" / "dkmq" / "shell_static_mct" / "ctria3_mitc3p_one_elem_shell_debug.ERR"
    shell_d_my = parse_err_named_matrix(err_path, "SHELL_D", 3, 3)
    shell_t_my = parse_err_named_matrix(err_path, "SHELL_T", 2, 2)
    cov_s_my = parse_err_named_matrix(err_path, "COV_S", 2, 2)
    drill_pen_my = parse_scalar(err_path, "MITC3P_DRILL_PEN")

    nodes = [
        Node(1, 0.0, 0.0, 0.0),
        Node(2, 0.5, 0.0, 0.0),
        Node(9, 0.5, 0.3, 0.0),
    ]
    elem = mitc3plus_mod.MITC3Plus(1, nodes, 1.07e7, 0.30, 0.1)
    _, _, _, x_local = mitc3plus_mod._local_frame(nodes)
    xy = x_local[:, :2]
    j = mitc3plus_mod._dh().T @ xy
    j_inv = np.linalg.inv(j)

    cb_py = mitc3plus_mod._C_bending(1.07e7, 0.30, 0.1)
    cs_py = mitc3plus_mod._C_shear(1.07e7, 0.30, 0.1)
    cov_s_py = j_inv.T @ cs_py @ j_inv
    drill_pen_py = elem.drilling_penalty

    print("MITC3+ material block probe")
    print(f"  shell_d_rel_diff   = {rel_diff(shell_d_my, cb_py):.9E}")
    print(f"  shell_t_rel_diff   = {rel_diff(shell_t_my, cs_py):.9E}")
    print(f"  cov_s_rel_diff     = {rel_diff(cov_s_my, cov_s_py):.9E}")
    print(f"  drill_pen_my       = {drill_pen_my:.9E}")
    print(f"  drill_pen_py       = {drill_pen_py:.9E}")
    print(f"  drill_pen_rel_diff = {abs(drill_pen_my - drill_pen_py) / abs(drill_pen_py):.9E}")
    print("  shell_d_my =")
    print(shell_d_my)
    print("  cb_py =")
    print(cb_py)
    print("  shell_t_my =")
    print(shell_t_my)
    print("  cs_py =")
    print(cs_py)
    print("  cov_s_my =")
    print(cov_s_my)
    print("  cov_s_py =")
    print(cov_s_py)


if __name__ == "__main__":
    main()

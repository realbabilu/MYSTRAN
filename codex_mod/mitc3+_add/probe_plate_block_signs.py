from __future__ import annotations

import itertools
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
            rows.append(row_vals[:ndof])
        else:
            i += 1
    return np.array(rows, dtype=float)


def rel_diff(a: np.ndarray, b: np.ndarray) -> float:
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def plate_idx() -> list[int]:
    out: list[int] = []
    for n in range(3):
        base = 6 * n
        out.extend([base + 2, base + 3, base + 4])
    return out


def apply_sign_case(k: np.ndarray, flips: tuple[int, int, int]) -> np.ndarray:
    idx = plate_idx()
    s = np.ones(len(idx))
    for n, flip in enumerate(flips):
        if flip & 1:
            s[3 * n + 1] *= -1.0
        if flip & 2:
            s[3 * n + 2] *= -1.0
    d = np.diag(s)
    kp = k[np.ix_(idx, idx)]
    return d @ kp @ d


def main() -> None:
    bug_path = (
        ROOT
        / "MYSTRANSolver-18.0.0"
        / "run_debug"
        / "dkmq"
        / "shell_static_mct"
        / "ctria3_mitc3p_one_elem_shell_debug.BUG"
    )
    k_mystran = parse_bug_full_ke(bug_path, 18)
    nodes = [
        Node(1, 0.0, 0.0, 0.0),
        Node(2, 0.5, 0.0, 0.0),
        Node(9, 0.5, 0.3, 0.0),
    ]
    for i, node in enumerate(nodes):
        node.dofs = list(range(6 * i, 6 * i + 6))
    k_python = MITC3Plus(1, nodes, 1.07e7, 0.30, 0.1).k_local()

    kp = k_python[np.ix_(plate_idx(), plate_idx())]
    print("MITC3+ plate-block sign probe")
    best = None
    for flips in itertools.product(range(4), repeat=3):
        km = apply_sign_case(k_mystran, flips)
        rd = rel_diff(km, kp)
        print(f"  flips={flips} rel_diff={rd:.9E}")
        if best is None or rd < best[1]:
            best = (flips, rd)
    print(f"best = flips={best[0]} rel_diff={best[1]:.9E}")


if __name__ == "__main__":
    main()

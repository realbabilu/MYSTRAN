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
from mitc3plus import MITC3Plus  # type: ignore


FLOAT_RE = re.compile(r"[-+]?\d+\.\d+E[-+]\d+")


def parse_bug_full_ke(path: Path, ndof: int) -> np.ndarray:
    lines = path.read_text().splitlines()
    start = None
    for i, line in enumerate(lines):
        if "KE element stiffness matrix in local element coordinate system" in line:
            start = i
            break
    rows = []
    i = start + 1
    while i < len(lines) and len(rows) < ndof:
        if lines[i].strip().startswith("Row"):
            vals = []
            i += 1
            while i < len(lines):
                nums = [float(x) for x in FLOAT_RE.findall(lines[i])]
                if not nums:
                    break
                vals.extend(nums)
                i += 1
                if len(vals) >= ndof:
                    break
            rows.append(vals[:ndof])
        else:
            i += 1
    return np.array(rows, dtype=float)


def plate_idx() -> list[int]:
    idx = []
    for n in range(3):
        base = 6 * n
        idx.extend([base + 2, base + 3, base + 4])
    return idx


def main() -> None:
    bug_path = ROOT / "MYSTRANSolver-18.0.0" / "run_debug" / "dkmq" / "shell_static_mct" / "ctria3_mitc3p_one_elem_shell_debug.BUG"
    k_my = parse_bug_full_ke(bug_path, 18)
    nodes = [Node(1, 0.0, 0.0, 0.0), Node(2, 0.5, 0.0, 0.0), Node(9, 0.5, 0.3, 0.0)]
    for i, node in enumerate(nodes):
        node.dofs = list(range(6 * i, 6 * i + 6))
    k_py = MITC3Plus(1, nodes, 1.07e7, 0.30, 0.1).k_local()

    idx = plate_idx()
    myp = k_my[np.ix_(idx, idx)]
    pyp = k_py[np.ix_(idx, idx)]
    print("MITC3+ plate-block eig/diag probe")
    print("  my_diag =", " ".join(f"{x:.6E}" for x in np.diag(myp)))
    print("  py_diag =", " ".join(f"{x:.6E}" for x in np.diag(pyp)))
    print("  my_eigs =", " ".join(f"{x:.6E}" for x in np.linalg.eigvalsh(myp)))
    print("  py_eigs =", " ".join(f"{x:.6E}" for x in np.linalg.eigvalsh(pyp)))


if __name__ == "__main__":
    main()

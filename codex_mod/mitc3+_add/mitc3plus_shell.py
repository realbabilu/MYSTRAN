"""
MITC3+ shell scaffold for the local ``core.py`` FEM engine.

This wrapper embeds the plate-side MITC3+ scaffold into shell-style
6 DOF/node slots:
    [ux, uy, uz, rx, ry, rz]

Active bending slots:
    [uz, rx, ry]

Current status:
- Source-backed shell kinematics and membrane coupling are not ported yet.
- The bending block comes from ``MITC3PlusBendingOnly``.
- This is suitable for smoke tests and gradual development beside ``core.py``.
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(r"D:\mystran2")
CLAUDEAI_DIR = ROOT / "pythonfem" / "claudeai"
MITC3PLUS_DIR = ROOT / "mitc3+"

for path in (CLAUDEAI_DIR, MITC3PLUS_DIR):
    path_str = str(path)
    if path_str not in sys.path:
        sys.path.insert(0, path_str)

from core import Element, Node  # noqa: E402
from mitc3plus_bending_only import MITC3PlusBendingOnly  # noqa: E402


def _tria_local_frame(coords3d: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    c = np.asarray(coords3d, dtype=float)
    origin = c[0]
    v1 = c[1] - c[0]
    e1 = v1 / np.linalg.norm(v1)
    v2 = c[2] - c[0]
    e3 = np.cross(v1, v2)
    e3 /= np.linalg.norm(e3)
    e2 = np.cross(e3, e1)
    e2 /= np.linalg.norm(e2)
    return e1, e2, e3, origin


def _project_to_local(coords3d: np.ndarray, e1: np.ndarray, e2: np.ndarray, origin: np.ndarray) -> np.ndarray:
    coords3d = np.asarray(coords3d, dtype=float)
    out = np.zeros((coords3d.shape[0], 2), dtype=float)
    for i, xyz in enumerate(coords3d):
        r = xyz - origin
        out[i, 0] = float(np.dot(r, e1))
        out[i, 1] = float(np.dot(r, e2))
    return out


def _shell_T_global_to_local(e1: np.ndarray, e2: np.ndarray, e3: np.ndarray, n_nodes: int) -> np.ndarray:
    lam = np.column_stack([e1, e2, e3])
    blk = np.zeros((6, 6), dtype=float)
    blk[0:3, 0:3] = lam.T
    blk[3:6, 3:6] = lam.T
    tmat = np.zeros((6 * n_nodes, 6 * n_nodes), dtype=float)
    for i in range(n_nodes):
        s = 6 * i
        tmat[s : s + 6, s : s + 6] = blk
    return tmat


def _embed_plate_block(k_plate: np.ndarray, n_nodes: int) -> np.ndarray:
    kout = np.zeros((6 * n_nodes, 6 * n_nodes), dtype=float)
    shell_slots = [6 * i + k for i in range(n_nodes) for k in (2, 3, 4)]
    for i, si in enumerate(shell_slots):
        for j, sj in enumerate(shell_slots):
            kout[si, sj] = k_plate[i, j]
    return kout


class MITC3PlusShell(Element):
    def __init__(
        self,
        eid: int,
        n1: Node,
        n2: Node,
        n3: Node,
        E: float,
        nu: float,
        t: float,
        rho: float = 0.0,
        kappa: float = 5.0 / 6.0,
        use_fallback_mitc3: bool = True,
    ):
        self.eid = eid
        self.nodes = [n1, n2, n3]
        self.E = float(E)
        self.nu = float(nu)
        self.t = float(t)
        self.rho = float(rho)
        self.kappa = float(kappa)
        self.use_fallback_mitc3 = bool(use_fallback_mitc3)

    @property
    def ndof_per_node(self) -> int:
        return 6

    def _local_plate_element(self) -> MITC3PlusBendingOnly:
        coords3d = np.array([n.coords for n in self.nodes], dtype=float)
        e1, e2, _e3, origin = _tria_local_frame(coords3d)
        coords2d = _project_to_local(coords3d, e1, e2, origin)
        local_nodes = [Node(i + 1, float(x), float(y), 0.0) for i, (x, y) in enumerate(coords2d)]
        return MITC3PlusBendingOnly(
            self.eid,
            local_nodes[0],
            local_nodes[1],
            local_nodes[2],
            self.E,
            self.nu,
            self.t,
            rho=self.rho,
            kappa=self.kappa,
            use_fallback_mitc3=self.use_fallback_mitc3,
        )

    def k_local(self) -> np.ndarray:
        k_plate = self._local_plate_element().k_local()
        coords3d = np.array([n.coords for n in self.nodes], dtype=float)
        e1, e2, e3, _origin = _tria_local_frame(coords3d)
        k_shell_local = _embed_plate_block(k_plate, 3)
        t_shell = _shell_T_global_to_local(e1, e2, e3, 3)
        return t_shell.T @ k_shell_local @ t_shell

    def T_matrix(self) -> np.ndarray:
        return np.eye(18, dtype=float)


__all__ = ["MITC3PlusShell"]

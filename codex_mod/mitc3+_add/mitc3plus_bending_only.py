"""
MITC3+ bending-only scaffold for the local ``core.py`` FEM engine.

Status:
- This file is a Phase-1 scaffold for the MITC3+ element of
  Lee, Lee, and Bathe (2014).
- The class layout, local DOF plan, and condensation hooks are prepared.
- The current stiffness used in ``k_local`` still falls back to the
  existing MITC3 baseline until the cubic-bubble rotation enrichment and
  the new assumed transverse shear strain field are ported in full.

Intended final element:
- 3 node plate-bending element
- external nodal DOFs: [w, rx, ry] at each node -> 9 external DOFs
- internal bubble rotational DOFs: 2
- internal bubble DOFs statically condensed on the element level
"""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


ROOT = Path(r"D:\mystran2")
CLAUDEAI_DIR = ROOT / "pythonfem" / "claudeai"

for path in (CLAUDEAI_DIR,):
    path_str = str(path)
    if path_str not in sys.path:
        sys.path.insert(0, path_str)

from core import Element, Node  # noqa: E402


def bubble_shape_cubic(l1: float, l2: float, l3: float) -> float:
    """Cubic bubble function used as the natural enrichment seed."""
    return 27.0 * l1 * l2 * l3


class MITC3Baseline(Element):
    """Local self-contained MITC3 baseline used by the MITC3+ scaffold."""

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
    ):
        self.eid = eid
        self.nodes = [n1, n2, n3]
        self.E = float(E)
        self.nu = float(nu)
        self.t = float(t)
        self.rho = float(rho)
        self.kappa = float(kappa)

    @property
    def ndof_per_node(self) -> int:
        return 3

    @property
    def D_b(self) -> np.ndarray:
        d = self.E * self.t**3 / (12.0 * (1.0 - self.nu**2))
        return d * np.array(
            [
                [1.0, self.nu, 0.0],
                [self.nu, 1.0, 0.0],
                [0.0, 0.0, 0.5 * (1.0 - self.nu)],
            ],
            dtype=float,
        )

    @property
    def G(self) -> float:
        return self.E / (2.0 * (1.0 + self.nu))

    @property
    def D_s(self) -> np.ndarray:
        return self.kappa * self.G * self.t * np.eye(2, dtype=float)

    def _area(self) -> float:
        x1, y1 = self.nodes[0].x, self.nodes[0].y
        x2, y2 = self.nodes[1].x, self.nodes[1].y
        x3, y3 = self.nodes[2].x, self.nodes[2].y
        return 0.5 * abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))

    def dN_dxdy(self) -> np.ndarray:
        x1, y1 = self.nodes[0].x, self.nodes[0].y
        x2, y2 = self.nodes[1].x, self.nodes[1].y
        x3, y3 = self.nodes[2].x, self.nodes[2].y
        two_a = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
        dN = np.zeros((3, 2), dtype=float)
        dN[0, 0] = (y2 - y3) / two_a
        dN[0, 1] = (x3 - x2) / two_a
        dN[1, 0] = (y3 - y1) / two_a
        dN[1, 1] = (x1 - x3) / two_a
        dN[2, 0] = (y1 - y2) / two_a
        dN[2, 1] = (x2 - x1) / two_a
        return dN

    def B_bending(self) -> np.ndarray:
        dN = self.dN_dxdy()
        bb = np.zeros((3, 9), dtype=float)
        for i in range(3):
            bb[0, 3 * i + 1] = dN[i, 0]
            bb[1, 3 * i + 2] = dN[i, 1]
            bb[2, 3 * i + 1] = dN[i, 1]
            bb[2, 3 * i + 2] = dN[i, 0]
        return bb

    def shear_strains(self, u_e: np.ndarray) -> tuple[float, float]:
        tying = ((0.5, 0.5, 0.0), (0.0, 0.5, 0.5), (0.5, 0.0, 0.5))
        dN = self.dN_dxdy()
        gamma = []
        for l1, l2, l3 in tying:
            n = np.array([l1, l2, l3], dtype=float)
            gxz = 0.0
            gyz = 0.0
            for i in range(3):
                gxz += dN[i, 0] * u_e[3 * i] + n[i] * u_e[3 * i + 1]
                gyz += dN[i, 1] * u_e[3 * i] + n[i] * u_e[3 * i + 2]
            gamma.append((gxz, gyz))
        gxz = (gamma[0][0] + gamma[1][0] + gamma[2][0]) / 3.0
        gyz = (gamma[0][1] + gamma[1][1] + gamma[2][1]) / 3.0
        return float(gxz), float(gyz)

    def k_local(self) -> np.ndarray:
        area = self._area()
        k = area * (self.B_bending().T @ self.D_b @ self.B_bending())
        bs = np.zeros((2, 9), dtype=float)
        for i in range(9):
            u_base = np.zeros(9, dtype=float)
            u_base[i] = 1.0
            gxz, gyz = self.shear_strains(u_base)
            bs[0, i] = gxz
            bs[1, i] = gyz
        k += area * (bs.T @ self.D_s @ bs)
        return k


class MITC3PlusBendingOnly(Element):
    """
    Plate-only MITC3+ scaffold with 3 DOF/node.

    External DOF ordering:
    [w1, rx1, ry1, w2, rx2, ry2, w3, rx3, ry3]

    Internal bubble DOF ordering placeholder:
    [alpha_bubble, beta_bubble]
    """

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
        d_twist: float = 1.0e-4,
    ):
        self.eid = eid
        self.nodes = [n1, n2, n3]
        self.E = float(E)
        self.nu = float(nu)
        self.t = float(t)
        self.rho = float(rho)
        self.kappa = float(kappa)
        self.use_fallback_mitc3 = bool(use_fallback_mitc3)
        self.d_twist = float(d_twist)

    @property
    def ndof_per_node(self) -> int:
        return 3

    @property
    def n_external_dof(self) -> int:
        return 9

    @property
    def n_internal_dof(self) -> int:
        return 2

    def T_matrix(self) -> np.ndarray:
        return np.eye(self.n_external_dof, dtype=float)

    def global_dof_indices(self) -> list[int]:
        dofs: list[int] = []
        for node in self.nodes:
            dofs.extend(node.dofs[:3])
        return dofs

    def _baseline_mitc3(self) -> MITC3Baseline:
        return MITC3Baseline(
            self.eid,
            self.nodes[0],
            self.nodes[1],
            self.nodes[2],
            self.E,
            self.nu,
            self.t,
            rho=self.rho,
            kappa=self.kappa,
        )

    def bubble_stiffness_blocks(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        k_full = self._enriched_full_stiffness()
        kaa = k_full[: self.n_external_dof, : self.n_external_dof]
        kab = k_full[: self.n_external_dof, self.n_external_dof :]
        kbb = k_full[self.n_external_dof :, self.n_external_dof :]
        return kaa, kab, kbb

    @staticmethod
    def condense_internal_dofs(kaa: np.ndarray, kab: np.ndarray, kbb: np.ndarray) -> np.ndarray:
        if kbb.size == 0:
            return kaa
        if np.allclose(kbb, 0.0):
            return kaa
        return kaa - kab @ np.linalg.solve(kbb, kab.T)

    def k_local(self) -> np.ndarray:
        k_base = self._baseline_mitc3().k_local()
        kaa, kab, kbb = self.bubble_stiffness_blocks()
        try:
            k_enriched = self.condense_internal_dofs(kaa, kab, kbb)
            if not np.all(np.isfinite(k_enriched)):
                raise FloatingPointError("non-finite enriched stiffness")
            if np.linalg.matrix_rank(k_enriched, tol=1.0e-10) < self.n_external_dof - 3:
                raise np.linalg.LinAlgError("rank-deficient enriched stiffness")
            return k_enriched
        except Exception:
            if self.use_fallback_mitc3:
                return k_base
            raise

    def m_local(self) -> np.ndarray | None:
        return None

    @property
    def D_b(self) -> np.ndarray:
        d = self.E * self.t**3 / (12.0 * (1.0 - self.nu**2))
        return d * np.array(
            [
                [1.0, self.nu, 0.0],
                [self.nu, 1.0, 0.0],
                [0.0, 0.0, 0.5 * (1.0 - self.nu)],
            ],
            dtype=float,
        )

    @property
    def D_s(self) -> np.ndarray:
        g = self.E / (2.0 * (1.0 + self.nu))
        return self.kappa * g * self.t * np.eye(2, dtype=float)

    def _area(self) -> float:
        x1, y1 = self.nodes[0].x, self.nodes[0].y
        x2, y2 = self.nodes[1].x, self.nodes[1].y
        x3, y3 = self.nodes[2].x, self.nodes[2].y
        return 0.5 * abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))

    def _dL_dxdy(self) -> np.ndarray:
        return self._baseline_mitc3().dN_dxdy()

    def _jacobian_rs_to_xy(self) -> np.ndarray:
        x1, y1 = self.nodes[0].x, self.nodes[0].y
        x2, y2 = self.nodes[1].x, self.nodes[1].y
        x3, y3 = self.nodes[2].x, self.nodes[2].y
        return np.array(
            [
                [x2 - x1, x3 - x1],
                [y2 - y1, y3 - y1],
            ],
            dtype=float,
        )

    @staticmethod
    def _shape_h(r: float, s: float) -> np.ndarray:
        return np.array([1.0 - r - s, r, s], dtype=float)

    def _shape_f(self, r: float, s: float) -> tuple[np.ndarray, float]:
        h = self._shape_h(r, s)
        f4 = bubble_shape_cubic(h[0], h[1], h[2])
        f = h - f4 / 3.0
        return f, f4

    def _shape_f_derivatives(self, r: float, s: float) -> tuple[np.ndarray, np.ndarray]:
        dL = self._dL_dxdy()
        l1, l2, l3 = self._shape_h(r, s)
        dl1dx, dl1dy = dL[0]
        dl2dx, dl2dy = dL[1]
        dl3dx, dl3dy = dL[2]
        df4dx = 27.0 * (dl1dx * l2 * l3 + l1 * dl2dx * l3 + l1 * l2 * dl3dx)
        df4dy = 27.0 * (dl1dy * l2 * l3 + l1 * dl2dy * l3 + l1 * l2 * dl3dy)
        dfdx = dL[:, 0] - df4dx / 3.0
        dfdy = dL[:, 1] - df4dy / 3.0
        return dfdx, dfdy

    def _bending_B_full(self, r: float, s: float) -> np.ndarray:
        dfdx, dfdy = self._shape_f_derivatives(r, s)
        l1, l2, l3 = self._shape_h(r, s)
        dL = self._dL_dxdy()
        dl1dx, dl1dy = dL[0]
        dl2dx, dl2dy = dL[1]
        dl3dx, dl3dy = dL[2]
        df4dx = 27.0 * (dl1dx * l2 * l3 + l1 * dl2dx * l3 + l1 * l2 * dl3dx)
        df4dy = 27.0 * (dl1dy * l2 * l3 + l1 * dl2dy * l3 + l1 * l2 * dl3dy)

        b = np.zeros((3, 11), dtype=float)
        for i in range(3):
            b[0, 3 * i + 1] = dfdx[i]
            b[1, 3 * i + 2] = dfdy[i]
            b[2, 3 * i + 1] = dfdy[i]
            b[2, 3 * i + 2] = dfdx[i]
        b[0, 9] = df4dx
        b[1, 10] = df4dy
        b[2, 9] = df4dy
        b[2, 10] = df4dx
        return b

    def _raw_shear_at(self, q: np.ndarray, r: float, s: float) -> tuple[float, float]:
        f, f4 = self._shape_f(r, s)
        dN = self._dL_dxdy()
        dwdx = float(np.dot(dN[:, 0], q[[0, 3, 6]]))
        dwdy = float(np.dot(dN[:, 1], q[[0, 3, 6]]))
        rx = float(np.dot(f, q[[1, 4, 7]]) + f4 * q[9])
        ry = float(np.dot(f, q[[2, 5, 8]]) + f4 * q[10])
        gamma_xz = dwdx + rx
        gamma_yz = dwdy + ry
        return gamma_xz, gamma_yz

    def _raw_covariant_shear_at(self, q: np.ndarray, r: float, s: float) -> tuple[float, float]:
        gamma_xz, gamma_yz = self._raw_shear_at(q, r, s)
        gamma_phys = np.array([gamma_xz, gamma_yz], dtype=float)
        e_cov = self._jacobian_rs_to_xy().T @ gamma_phys
        return float(e_cov[0]), float(e_cov[1])

    def _shear_B_full(self, r: float, s: float) -> np.ndarray:
        d = self.d_twist
        points = {
            "A": (1.0 / 6.0, 2.0 / 3.0),
            "B": (2.0 / 3.0, 1.0 / 6.0),
            "C": (1.0 / 6.0, 1.0 / 6.0),
            "D": (1.0 / 3.0 + d, 1.0 / 3.0 - 2.0 * d),
            "E": (1.0 / 3.0 - 2.0 * d, 1.0 / 3.0 + d),
            "F": (1.0 / 3.0 + d, 1.0 / 3.0 + d),
        }

        def eval_from_q(q: np.ndarray) -> tuple[float, float]:
            eA_rt, eA_st = self._raw_covariant_shear_at(q, *points["A"])
            eB_rt, eB_st = self._raw_covariant_shear_at(q, *points["B"])
            eC_rt, eC_st = self._raw_covariant_shear_at(q, *points["C"])
            eD_rt, eD_st = self._raw_covariant_shear_at(q, *points["D"])
            eE_rt, eE_st = self._raw_covariant_shear_at(q, *points["E"])
            eF_rt, eF_st = self._raw_covariant_shear_at(q, *points["F"])

            c_hat = (eF_rt - eD_rt) - (eF_st - eE_st)
            econst_rt = (2.0 / 3.0) * (eB_rt - 0.5 * eB_st) - (1.0 / 3.0) * (eC_rt - eC_st)
            econst_st = (2.0 / 3.0) * (eA_st - 0.5 * eA_rt) + (1.0 / 3.0) * (eC_rt - eC_st)
            ehat_rt = econst_rt + (1.0 / 3.0) * c_hat * (3.0 * s - 1.0)
            ehat_st = econst_st + (1.0 / 3.0) * c_hat * (1.0 - 3.0 * r)

            # The paper defines the assumed field in covariant r/s components.
            # Convert back to physical x/y shear components for the plate
            # constitutive law used by this local bending-only model.
            gamma_phys = np.linalg.solve(self._jacobian_rs_to_xy().T, np.array([ehat_rt, ehat_st], dtype=float))
            return float(gamma_phys[0]), float(gamma_phys[1])

        bs = np.zeros((2, 11), dtype=float)
        for i in range(11):
            q = np.zeros(11, dtype=float)
            q[i] = 1.0
            bs[0, i], bs[1, i] = eval_from_q(q)
        return bs

    def _enriched_full_stiffness(self) -> np.ndarray:
        area = self._area()
        gauss = (
            (1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0),
            (2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0),
            (1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0),
        )
        k = np.zeros((11, 11), dtype=float)
        jac = 2.0 * area
        for r, s, w in gauss:
            bb = self._bending_B_full(r, s)
            bs = self._shear_B_full(r, s)
            fac = jac * w
            k += fac * (bb.T @ self.D_b @ bb + bs.T @ self.D_s @ bs)
        return 0.5 * (k + k.T)


__all__ = ["MITC3PlusBendingOnly", "bubble_shape_cubic"]

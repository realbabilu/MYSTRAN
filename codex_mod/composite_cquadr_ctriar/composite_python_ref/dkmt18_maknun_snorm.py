"""
dkmt18_maknun_snorm.py — DKMT18 Maknun-style triangular shell element with SNORM/nodal normals
=================================================================

Experimental paper-oriented implementation for the local ``core.py`` / ``core_sparse.py``
FEM engine.  This is intentionally separate from ``dkmt18_allman.py``.

Target paper:
  I.J. Maknun, I. Katili, A. Ibrahimbegovic, A.M. Katili (2020),
  "A new triangular shell element for composites accounting for shear deformation",
  Composite Structures 243, 112214.

Scope of this file
------------------
- 3-node triangular shell, 6 DOF/node: [U,V,W, RX,RY,RZ] in global axes.
- Naghdi/Reissner-Mindlin shell-style covariant geometry.
- Membrane matrix Bm built from Co and local surface basis t1,t2.
- Curvature matrix Bb includes Maknun curvature coupling ``bc * u`` plus
  DKMT incomplete-quadratic rotational correction.
- Transverse shear matrix Bs uses the DKMT discrete shear influence factor phi_k.
- Fictitious stiffness for normal rotation theta_z / drilling-like modes.

Important notes
---------------
1. This is not the older flat-shell superposition:
       Allman membrane + DKMT plate + zero coupling.
   Geometry coupling enters through nodal normals and the ``bc`` matrix.

2. For a flat triangle with a constant normal, ``bc = 0`` and this element
   reduces close to a DKMT flat-shell form.

3. For curved shells, pass nodal normals (SNORM/director) via ``nodal_normals``.
   Without them, the element uses its flat element normal at all 3 nodes, so
   the curved-shell coupling vanishes.

4. Composite laminate A/B/D matrices are not implemented yet.  The current
   constitutive law is isotropic homogeneous, so material membrane-bending
   coupling B=0.  Geometric membrane-bending coupling is still present through
   the curvature matrix.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Optional

import numpy as np

from core import Element, Node


# 3-point rule on reference triangle: area(ref) = 1/2, weights sum to 1/2
_TRI3 = (
    (0.5, 0.0, 1.0 / 6.0),
    (0.0, 0.5, 1.0 / 6.0),
    (0.5, 0.5, 1.0 / 6.0),
)


def shape_N(xi: float, eta: float) -> np.ndarray:
    return np.array([1.0 - xi - eta, xi, eta], dtype=float)


def shape_dN() -> np.ndarray:
    # rows: d/dxi, d/deta; cols: nodes 1,2,3
    return np.array([[-1.0, 1.0, 0.0], [-1.0, 0.0, 1.0]], dtype=float)


def shape_P(xi: float, eta: float) -> np.ndarray:
    lam = 1.0 - xi - eta
    return np.array([4.0 * lam * xi, 4.0 * xi * eta, 4.0 * lam * eta], dtype=float)


def shape_dP(xi: float, eta: float) -> np.ndarray:
    lam = 1.0 - xi - eta
    return np.array(
        [
            [4.0 * (lam - xi), 4.0 * eta, -4.0 * eta],
            [-4.0 * xi, 4.0 * xi, 4.0 * (lam - eta)],
        ],
        dtype=float,
    )


def _unit(v: np.ndarray, fallback: np.ndarray | None = None) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n > 1.0e-14:
        return np.asarray(v, dtype=float) / n
    if fallback is None:
        fallback = np.array([0.0, 0.0, 1.0], dtype=float)
    return np.asarray(fallback, dtype=float).copy()


def _RN(nx: float, ny: float, nz: float) -> np.ndarray:
    """Cross-product matrix used in Maknun/DKMQ-style rotation coupling."""
    return np.array(
        [[0.0, nz, -ny], [-nz, 0.0, nx], [ny, -nx, 0.0]],
        dtype=float,
    )


@dataclass
class DKMT18MaknunOptions:
    fictitious_scale: float = 1.0e-3
    use_fictitious: bool = True
    use_geometric_coupling: bool = True
    use_dkmt_shear: bool = True


class DKMT18MaknunSNORM(Element):
    """Maknun-style DKMT18 triangular shell element with nodal normals/SNORM.

    This class is intentionally the SNORM-enabled variant: pass nodal_normals
    for curved shells. If nodal_normals is None, it falls back to element normal,
    but that should be considered DKMT18_FLAT behavior.
    

    Parameters
    ----------
    eid : int
        Element id.
    nodes : list[Node]
        Three corner nodes.
    E, nu, h : float
        Isotropic material constants and thickness.
    nodal_normals : dict[int, ndarray] or ndarray, optional
        Nodal directors/SNORM.  If omitted, element normal is used for all nodes.
        Pass analytic or geometric nodal normals for curved shells.
    """

    ndof_per_node = 6

    def __init__(
        self,
        eid: int,
        nodes: list[Node],
        E: float,
        nu: float,
        h: float,
        kappa: float = 5.0 / 6.0,
        nodal_normals: Optional[Dict[int, np.ndarray] | np.ndarray] = None,
        options: Optional[DKMT18MaknunOptions] = None,
    ):
        if len(nodes) != 3:
            raise ValueError("DKMT18Maknun requires exactly 3 nodes")
        self.eid = int(eid)
        self.nodes = list(nodes)
        self.E = float(E)
        self.nu = float(nu)
        self.h = float(h)
        self.kappa = float(kappa)
        self.options = options or DKMT18MaknunOptions()
        self._xyz = np.array([n.coords for n in nodes], dtype=float)
        self._side_pairs = [(0, 1), (1, 2), (2, 0)]  # sides 4,5,6
        self._normals = self._prepare_normals(nodal_normals)
        self._K18_cache: np.ndarray | None = None

    # ---- core.py interface -------------------------------------------------

    def k_local(self) -> np.ndarray:
        # The element is assembled directly in global DOF coordinates.
        return self.k_global()

    def T_matrix(self) -> np.ndarray:
        return np.eye(18)

    def k_global(self) -> np.ndarray:
        if self._K18_cache is None:
            self._K18_cache = self._compute_stiffness()
        return self._K18_cache

    def global_dof_indices(self) -> list[int]:
        dofs: list[int] = []
        for node in self.nodes:
            dofs.extend(node.dofs[:6])
        return dofs

    # ---- material ----------------------------------------------------------

    def _Hm(self) -> np.ndarray:
        c = self.E * self.h / (1.0 - self.nu**2)
        v = self.nu
        return c * np.array([[1.0, v, 0.0], [v, 1.0, 0.0], [0.0, 0.0, (1.0 - v) / 2.0]])

    def _Hb(self) -> np.ndarray:
        c = self.E * self.h**3 / (12.0 * (1.0 - self.nu**2))
        v = self.nu
        return c * np.array([[1.0, v, 0.0], [v, 1.0, 0.0], [0.0, 0.0, (1.0 - v) / 2.0]])

    def _Hs(self) -> np.ndarray:
        c = self.kappa * self.E * self.h / (2.0 * (1.0 + self.nu))
        return c * np.eye(2)

    # ---- geometry / normals ------------------------------------------------

    def _element_normal(self) -> np.ndarray:
        x = self._xyz
        return _unit(np.cross(x[1] - x[0], x[2] - x[0]))

    def _prepare_normals(self, nodal_normals) -> np.ndarray:
        en = self._element_normal()
        out = np.zeros((3, 3), dtype=float)
        if nodal_normals is None:
            out[:, :] = en
            return out
        if isinstance(nodal_normals, dict):
            for i, nd in enumerate(self.nodes):
                out[i] = _unit(np.asarray(nodal_normals.get(nd.nid, en), dtype=float), en)
        else:
            arr = np.asarray(nodal_normals, dtype=float)
            if arr.shape != (3, 3):
                raise ValueError("nodal_normals array must have shape (3,3)")
            for i in range(3):
                out[i] = _unit(arr[i], en)
        # Avoid accidental reversed SNORMs for this element orientation.
        for i in range(3):
            if float(out[i] @ en) < 0.0:
                out[i] *= -1.0
        return out

    def _geometry_at(self, xi: float, eta: float):
        dN = shape_dN()
        xyz = self._xyz
        normals = self._normals

        a1 = dN[0] @ xyz
        a2 = dN[1] @ xyz
        axb = np.cross(a1, a2)
        J = float(np.linalg.norm(axb))
        if J <= 1.0e-14:
            raise ValueError(f"DKMT18Maknun eid={self.eid}: degenerate triangle")
        n_vec = axb / J

        a11 = float(a1 @ a1)
        a12 = float(a1 @ a2)
        a22 = float(a2 @ a2)
        amat = np.array([[a11, a12], [a12, a22]], dtype=float)
        inva = np.linalg.inv(amat)
        a1c = inva[0, 0] * a1 + inva[0, 1] * a2
        a2c = inva[1, 0] * a1 + inva[1, 1] * a2

        k_hat = np.array([0.0, 0.0, 1.0])
        t1 = np.cross(n_vec, k_hat)
        if np.linalg.norm(t1) < 1.0e-10:
            t1 = np.cross(n_vec, np.array([0.0, 1.0, 0.0]))
        t1 = _unit(t1)
        t2 = _unit(np.cross(n_vec, t1))

        Co = np.array(
            [[a1c @ t1, a1c @ t2], [a2c @ t1, a2c @ t2]],
            dtype=float,
        )

        # Curvature coupling: derivative of interpolated nodal normals.
        n_xi = dN[0] @ normals
        n_eta = dN[1] @ normals
        bn11 = float(a1c @ n_xi)
        bn12 = float(a1c @ n_eta)
        bn21 = float(a2c @ n_xi)
        bn22 = float(a2c @ n_eta)
        bn_hat = np.array([[bn22, -bn12], [-bn21, bn11]], dtype=float)
        bc_mat = bn_hat @ Co if self.options.use_geometric_coupling else np.zeros((2, 2))

        return t1, t2, n_vec, J, Co, bc_mat, normals

    # ---- DKMT side equations ----------------------------------------------

    def _side_lengths_phi(self) -> tuple[np.ndarray, np.ndarray]:
        Ls = np.zeros(3, dtype=float)
        phi = np.zeros(3, dtype=float)
        for k, (i, j) in enumerate(self._side_pairs):
            Ls[k] = float(np.linalg.norm(self._xyz[j] - self._xyz[i]))
            if Ls[k] < 1.0e-14 or not self.options.use_dkmt_shear:
                phi[k] = 0.0
            else:
                phi[k] = 2.0 / (self.kappa * (1.0 - self.nu)) * (self.h / Ls[k]) ** 2
        return Ls, phi

    def _Au_matrix(self) -> np.ndarray:
        """3×18 matrix mapping nodal DOF to side rotational increments."""
        xyz = self._xyz
        normals = self._normals
        Au = np.zeros((3, 18), dtype=float)
        for k, (i, j) in enumerate(self._side_pairs):
            xji = xyz[j] - xyz[i]
            L = float(np.linalg.norm(xji))
            if L < 1.0e-14:
                continue
            tsk = xji / L
            nk = _unit(0.5 * (normals[i] + normals[j]), self._element_normal())
            RNi = _RN(*normals[i])
            RNj = _RN(*normals[j])
            ci = 6 * i
            cj = 6 * j
            # shell-normal transverse displacement difference along the side
            Au[k, ci : ci + 3] += -nk / L
            Au[k, cj : cj + 3] += nk / L
            # average tangential rotation along the side
            Au[k, ci + 3 : ci + 6] += 0.5 * (RNi.T @ tsk)
            Au[k, cj + 3 : cj + 6] += 0.5 * (RNj.T @ tsk)
        return Au

    def _Adelta_inv_Au(self) -> np.ndarray:
        _Ls, phi = self._side_lengths_phi()
        # sign convention follows the DKMT triangular side constraint.
        diag = -(2.0 / 3.0) * (1.0 + phi)
        return np.diag(1.0 / diag) @ self._Au_matrix()

    def _side_CS(self, t1: np.ndarray, t2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        C = np.zeros(3, dtype=float)
        S = np.zeros(3, dtype=float)
        for k, (i, j) in enumerate(self._side_pairs):
            tsk = _unit(self._xyz[j] - self._xyz[i], t1)
            C[k] = float(t1 @ tsk)
            S[k] = float(t2 @ tsk)
        return C, S

    @staticmethod
    def _Bsg_from_CS(C: np.ndarray, S: np.ndarray) -> np.ndarray:
        # Same topology as DKMT sides 4=(0-1), 5=(1-2), 6=(2-0).
        A1 = C[0] * S[2] - C[2] * S[0]
        A2 = C[1] * S[0] - C[0] * S[1]
        A3 = C[2] * S[1] - C[1] * S[2]
        eps = 1.0e-14
        if abs(A1) < eps or abs(A2) < eps or abs(A3) < eps:
            return np.zeros((2, 3), dtype=float)
        Bsg = np.array(
            [
                [(S[1] / A2 - S[2] / A1), (S[2] / A3 - S[0] / A2), (S[0] / A1 - S[1] / A3)],
                [(C[2] / A1 - C[1] / A2), (C[0] / A2 - C[2] / A3), (C[1] / A3 - C[0] / A1)],
            ],
            dtype=float,
        )
        return (2.0 / 3.0) * Bsg

    # ---- B matrices ---------------------------------------------------------

    def _Bm_at(self, xi, eta, t1, t2, J, Co) -> np.ndarray:
        dN = shape_dN()
        Nix = dN[0] * Co[0, 0] + dN[1] * Co[1, 0]
        Niy = dN[0] * Co[0, 1] + dN[1] * Co[1, 1]
        Bm = np.zeros((3, 18), dtype=float)
        for i in range(3):
            c = 6 * i
            for d in range(3):
                Bm[0, c + d] = t1[d] * Nix[i]
                Bm[1, c + d] = t2[d] * Niy[i]
                Bm[2, c + d] = t1[d] * Niy[i] + t2[d] * Nix[i]
        return Bm

    def _Bb_at(self, xi, eta, t1, t2, n_vec, Co, bc_mat, normals, Adelta_inv_Au) -> np.ndarray:
        dN = shape_dN()
        dP = shape_dP(xi, eta)
        Nix = dN[0] * Co[0, 0] + dN[1] * Co[1, 0]
        Niy = dN[0] * Co[0, 1] + dN[1] * Co[1, 1]
        Pkx = dP[0] * Co[0, 0] + dP[1] * Co[1, 0]
        Pky = dP[0] * Co[0, 1] + dP[1] * Co[1, 1]
        Nbc1 = dN[0] * bc_mat[0, 0] + dN[1] * bc_mat[1, 0]
        Nbc2 = dN[0] * bc_mat[0, 1] + dN[1] * bc_mat[1, 1]

        Bb_beta = np.zeros((3, 18), dtype=float)
        for i in range(3):
            c = 6 * i
            RNi = _RN(*normals[i])
            V1i = np.array([t1 @ RNi[:, 0], t1 @ RNi[:, 1], t1 @ RNi[:, 2]], dtype=float)
            V2i = np.array([t2 @ RNi[:, 0], t2 @ RNi[:, 1], t2 @ RNi[:, 2]], dtype=float)
            # curvature coupling from translations: bc * u
            for d in range(3):
                Bb_beta[0, c + d] = t1[d] * Nbc1[i]
                Bb_beta[1, c + d] = t2[d] * Nbc2[i]
                Bb_beta[2, c + d] = t1[d] * Nbc2[i] + t2[d] * Nbc1[i]
            # rotation contribution
            Bb_beta[0, c + 3 : c + 6] = V1i * Nix[i]
            Bb_beta[1, c + 3 : c + 6] = V2i * Niy[i]
            Bb_beta[2, c + 3 : c + 6] = V1i * Niy[i] + V2i * Nix[i]

        C, S = self._side_CS(t1, t2)
        BbD = np.zeros((3, 3), dtype=float)
        for k, (i, j) in enumerate(self._side_pairs):
            tsk = _unit(self._xyz[j] - self._xyz[i], t1)
            t1s = float(t1 @ tsk)
            t2s = float(t2 @ tsk)
            BbD[0, k] = t1s * Pkx[k]
            BbD[1, k] = t2s * Pky[k]
            BbD[2, k] = t1s * Pky[k] + t2s * Pkx[k]
        return Bb_beta + BbD @ Adelta_inv_Au

    def _Bs_at(self, t1, t2, Adelta_inv_Au) -> np.ndarray:
        if not self.options.use_dkmt_shear:
            return np.zeros((2, 18), dtype=float)
        C, S = self._side_CS(t1, t2)
        Bsg = self._Bsg_from_CS(C, S)
        _Ls, phi = self._side_lengths_phi()
        return Bsg @ np.diag(phi) @ Adelta_inv_Au

    # ---- fictitious drilling / normal rotation -----------------------------

    def _k_fictitious(self) -> np.ndarray:
        if not self.options.use_fictitious or self.options.fictitious_scale <= 0.0:
            return np.zeros((18, 18), dtype=float)
        scale = float(self.options.fictitious_scale)
        K = np.zeros((18, 18), dtype=float)
        normals = self._normals
        # Two parts: gradient of theta_n and small mass-like theta_n term.
        alpha_grad = scale * self.E * self.h**3 / 12.0
        alpha_mac = scale * self.E * self.h / (2.0 * (1.0 + self.nu))
        for xi, eta, w in _TRI3:
            N = shape_N(xi, eta)
            dN = shape_dN()
            t1, t2, n_vec, J, Co, bc_mat, _ = self._geometry_at(xi, eta)
            fac = w * J
            Nix = dN[0] * Co[0, 0] + dN[1] * Co[1, 0]
            Niy = dN[0] * Co[0, 1] + dN[1] * Co[1, 1]
            Gth = np.zeros((2, 18), dtype=float)
            Hth = np.zeros(18, dtype=float)
            for i in range(3):
                c = 6 * i
                ni = normals[i]
                Gth[0, c + 3 : c + 6] = Nix[i] * ni
                Gth[1, c + 3 : c + 6] = Niy[i] * ni
                Hth[c + 3 : c + 6] = N[i] * ni
            K += fac * (alpha_grad * (Gth.T @ Gth) + alpha_mac * np.outer(Hth, Hth))
        return K

    # ---- stiffness ---------------------------------------------------------

    def _compute_stiffness(self) -> np.ndarray:
        Hm, Hb, Hs = self._Hm(), self._Hb(), self._Hs()
        Adelta_inv_Au = self._Adelta_inv_Au()
        K = np.zeros((18, 18), dtype=float)
        for xi, eta, w in _TRI3:
            t1, t2, n_vec, J, Co, bc_mat, normals = self._geometry_at(xi, eta)
            Bm = self._Bm_at(xi, eta, t1, t2, J, Co)
            Bb = self._Bb_at(xi, eta, t1, t2, n_vec, Co, bc_mat, normals, Adelta_inv_Au)
            Bs = self._Bs_at(t1, t2, Adelta_inv_Au)
            fac = w * J
            K += fac * (Bm.T @ Hm @ Bm + Bb.T @ Hb @ Bb + Bs.T @ Hs @ Bs)
        K += self._k_fictitious()
        return 0.5 * (K + K.T)

    def area(self) -> float:
        return 0.5 * float(np.linalg.norm(np.cross(self._xyz[1] - self._xyz[0], self._xyz[2] - self._xyz[0])))

    def __repr__(self):
        return f"DKMT18Maknun(eid={self.eid}, nodes={[n.nid for n in self.nodes]}, h={self.h})"


# ---- reusable normal helpers ----------------------------------------------


def compute_geometric_nodal_normals(nodes: dict[int, Node], tri_conn: Iterable[tuple[int, int, int]]) -> dict[int, np.ndarray]:
    acc: dict[int, np.ndarray] = {nid: np.zeros(3, dtype=float) for nid in nodes}
    for tri in tri_conn:
        p = [nodes[nid].coords for nid in tri]
        area_vec = np.cross(p[1] - p[0], p[2] - p[0])
        for nid in tri:
            acc[nid] += area_vec
    out: dict[int, np.ndarray] = {}
    for nid, v in acc.items():
        out[nid] = _unit(v)
    return out


# Backward-friendly alias
DKMT18Maknun = DKMT18MaknunSNORM

__all__ = ["DKMT18MaknunSNORM", "DKMT18Maknun", "DKMT18MaknunOptions", "compute_geometric_nodal_normals"]

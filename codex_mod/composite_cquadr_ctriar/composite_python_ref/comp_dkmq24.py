"""
comp_dkmq24.py

Composite DKMQ24 quadrilateral shell element.

Basis:
  - dkmq24_katili.DKMQ24
  - laminate_utils.Laminate

This is the first composite-capable variant for CQUADR-like usage.

Key difference from isotropic DKMQ24:
  K = Bm.T A Bm
    + Bm.T B Bb + Bb.T B Bm
    + Bb.T D Bb
    + Bs.T As Bs
    + drilling/fictitious stiffness

So unsymmetric laminates with B != 0 are supported at static stiffness level.
"""

from __future__ import annotations

import numpy as np

from laminate_utils import Laminate

try:
    from dkmq24_katili import DKMQ24, GAUSS_POINTS
except Exception:
    from dkmq24_element_v3 import DKMQ24, GAUSS_POINTS


class CompDKMQ24(DKMQ24):
    """
    Composite DKMQ24 shell element.

    Parameters
    ----------
    eid, nodes:
        Same as DKMQ24.
    laminate:
        Laminate object with A/B/D/As matrices.
    kappa:
        Kept for API compatibility.  Shear stiffness comes from laminate.As.
    use_mindlin:
        If False, transverse shear block is skipped.
    """

    def __init__(self, eid: int, nodes: list, laminate: Laminate,
                 kappa: float = 5.0 / 6.0, use_mindlin: bool = True):
        self.laminate = laminate
        # E_ref/nu_ref are used only by inherited fictitious drilling scale.
        super().__init__(
            eid=eid,
            nodes=nodes,
            E=laminate.E_ref,
            nu=laminate.nu_ref,
            h=laminate.h,
            kappa=kappa,
            use_mindlin=use_mindlin,
        )

    def _Hm(self) -> np.ndarray:
        return np.asarray(self.laminate.A, dtype=float)

    def _Hb(self) -> np.ndarray:
        return np.asarray(self.laminate.D, dtype=float)

    def _Hs(self) -> np.ndarray:
        if not self.use_mindlin:
            return np.zeros((2, 2), dtype=float)
        return np.asarray(self.laminate.As, dtype=float)

    def _Hmb(self) -> np.ndarray:
        return np.asarray(self.laminate.B, dtype=float)

    def _compute_stiffness(self) -> np.ndarray:
        A = self._Hm()
        B = self._Hmb()
        D = self._Hb()
        As = self._Hs()

        Au = self._Au_matrix()
        Adelta = self._Adelta_matrix()
        Adelta_inv = np.diag(1.0 / np.diag(Adelta))
        Adelta_inv_Au = Adelta_inv @ Au

        K = np.zeros((24, 24), dtype=float)

        for xi, eta, w in GAUSS_POINTS:
            t1, t2, n_vec, J, Co, bc_mat, normals = self._geometry_at(xi, eta)

            Bm = self._Bm_at(xi, eta, t1, t2, n_vec, J, Co, bc_mat, normals)
            Bb = self._Bb_at(xi, eta, t1, t2, n_vec, J, Co, bc_mat, normals, Adelta_inv_Au)
            Bs = self._Bs_at(xi, eta, t1, t2, n_vec, J, Co, bc_mat, normals, Adelta_inv_Au)

            K += w * J * (
                Bm.T @ A @ Bm
                + Bm.T @ B @ Bb
                + Bb.T @ B @ Bm
                + Bb.T @ D @ Bb
                + Bs.T @ As @ Bs
            )

        K += self._k_fictitious()
        return 0.5 * (K + K.T)

    def recompute_stiffness(self) -> None:
        self._K24 = self._compute_stiffness()


# Backward/friendly aliases
comp_DKMQ24 = CompDKMQ24
comp_dkmq24 = CompDKMQ24
CompDkmq24 = CompDKMQ24

__all__ = ["CompDKMQ24", "comp_DKMQ24", "comp_dkmq24", "CompDkmq24"]

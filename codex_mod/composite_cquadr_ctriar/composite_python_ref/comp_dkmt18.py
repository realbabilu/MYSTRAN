"""
comp_dkmt18.py

Composite DKMT18 Maknun triangular shell element.

Basis:
  - dkmt18_maknun_snorm.DKMT18MaknunSNORM
  - laminate_utils.Laminate

This is the first composite-capable variant for CTRIAR-like usage.

Key difference from isotropic DKMT18:
  K = Bm.T A Bm
    + Bm.T B Bb + Bb.T B Bm
    + Bb.T D Bb
    + Bs.T As Bs
    + fictitious normal-rotation stiffness

So unsymmetric laminates with B != 0 are supported at static stiffness level.
"""

from __future__ import annotations

# Sandbox fallback: some uploaded workspaces store core.py as core(2).py.
from pathlib import Path as _Path
import sys as _sys
_HERE = _Path(__file__).resolve().parent
if "core" not in _sys.modules and not (_HERE / "core.py").exists() and (_HERE / "core(2).py").exists():
    import importlib.util as _importlib_util
    _spec = _importlib_util.spec_from_file_location("core", _HERE / "core(2).py")
    _core_mod = _importlib_util.module_from_spec(_spec)
    _sys.modules["core"] = _core_mod
    _spec.loader.exec_module(_core_mod)

import numpy as np

from laminate_utils import Laminate
from dkmt18_maknun_snorm import DKMT18MaknunSNORM, DKMT18MaknunOptions, _TRI3


class CompDKMT18(DKMT18MaknunSNORM):
    """
    Composite DKMT18 triangular shell element.

    Parameters
    ----------
    eid, nodes:
        Same as DKMT18MaknunSNORM.
    laminate:
        Laminate object with A/B/D/As matrices.
    nodal_normals:
        Same as DKMT18MaknunSNORM.  For curved shell, pass SNORM/directors.
    options:
        Same as DKMT18MaknunSNORM.
    """

    def __init__(
        self,
        eid: int,
        nodes: list,
        laminate: Laminate,
        kappa: float = 5.0 / 6.0,
        nodal_normals=None,
        options: DKMT18MaknunOptions | None = None,
    ):
        self.laminate = laminate
        super().__init__(
            eid=eid,
            nodes=nodes,
            E=laminate.E_ref,
            nu=laminate.nu_ref,
            h=laminate.h,
            kappa=kappa,
            nodal_normals=nodal_normals,
            options=options,
        )

    def _Hm(self) -> np.ndarray:
        return np.asarray(self.laminate.A, dtype=float)

    def _Hb(self) -> np.ndarray:
        return np.asarray(self.laminate.D, dtype=float)

    def _Hs(self) -> np.ndarray:
        if not self.options.use_dkmt_shear:
            return np.zeros((2, 2), dtype=float)
        return np.asarray(self.laminate.As, dtype=float)

    def _Hmb(self) -> np.ndarray:
        return np.asarray(self.laminate.B, dtype=float)

    def _compute_stiffness(self) -> np.ndarray:
        A = self._Hm()
        B = self._Hmb()
        D = self._Hb()
        As = self._Hs()

        Adelta_inv_Au = self._Adelta_inv_Au()
        K = np.zeros((18, 18), dtype=float)

        for xi, eta, w in _TRI3:
            t1, t2, n_vec, J, Co, bc_mat, normals = self._geometry_at(xi, eta)
            Bm = self._Bm_at(xi, eta, t1, t2, J, Co)
            Bb = self._Bb_at(xi, eta, t1, t2, n_vec, Co, bc_mat, normals, Adelta_inv_Au)
            Bs = self._Bs_at(t1, t2, Adelta_inv_Au)
            fac = w * J

            K += fac * (
                Bm.T @ A @ Bm
                + Bm.T @ B @ Bb
                + Bb.T @ B @ Bm
                + Bb.T @ D @ Bb
                + Bs.T @ As @ Bs
            )

        K += self._k_fictitious()
        return 0.5 * (K + K.T)


# Backward/friendly aliases
comp_DKMT18 = CompDKMT18
comp_dkmt18 = CompDKMT18
CompDkmt18 = CompDKMT18

__all__ = ["CompDKMT18", "comp_DKMT18", "comp_dkmt18", "CompDkmt18"]

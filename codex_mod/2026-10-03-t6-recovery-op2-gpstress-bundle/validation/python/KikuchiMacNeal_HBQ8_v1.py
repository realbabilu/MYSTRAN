"""
KikuchiMacNeal_HBQ8_v1.py
============================
HBQ8_ShellElement's architecture (Darilmaz & Kumbasar 2006, plain
direct-Cartesian Bm/Bb/Bs/Bdrill, no competing assumed-strain scheme --
the simplest of the four Q8 hybrids built so far) with FIELD
interpolation swapped for Kikuchi (1999) or MacNeal-Harder (1992)
modified shape functions.

Lowest-risk hybrid in this family: HBQ8's own _compute_Bm/_Bb/_Bs/
_Bdrill all use plain Jinv@dN Cartesian differentiation directly --
same flat architecture as MacNealQ8_1992_native itself. Unlike MITC8
(needed its ILS averaging bypassed) or ANS8 (needed its beta+delta
tying bypassed), there is no competing locking-treatment machinery
here to conflict with the modified shape functions -- the substitution
is direct.

Orientation handling: per KikuchiMacNeal_MITC8_v1's lesson (the same
"+Z tiebreak" fix that ANS8 needed actively BROKE MITC8, because each
element's own inherited frame-construction convention differs), this
is checked empirically below rather than assumed.
"""
import numpy as np
from HBQ8_ShellElement import HBQ8_ShellElement
from macneal_kikuchi_shapefun import kikuchi_shape_functions, macneal_harder_shape_functions


class KikuchiMacNeal_HBQ8_v1(HBQ8_ShellElement):

    def __init__(self, eid, nodes, E, nu, h, kshear=5./6., beta_drill=1e-3, variant='kikuchi'):
        if len(nodes) != 8:
            raise AssertionError(f"requires 8 nodes, got {len(nodes)}")
        self.variant = variant
        super().__init__(eid, nodes, E, nu, h, kshear=kshear, beta_drill=beta_drill)
        e1, e2, e3 = self.ex, self.ey, self.ez
        self.E1, self.E2, self.E3 = e1, e2, e3
        coords = np.array([[n.x, n.y, n.z] for n in self.nodes])
        centroid = coords[:4].mean(axis=0)
        rel = coords - centroid
        self._xy8 = np.column_stack([rel @ e1, rel @ e2])

    def _local_basis_at_point(self, r, s):
        # Same fix, same root cause as KikuchiMacNeal_ANS8_v1 and
        # KikuchiMacNeal_MITC8_v1's overrides (see those files'
        # comments for the full derivation): e1,e2 are the working,
        # unmodified self.E1/E2 (no +Z tiebreak applied here at all --
        # unlike ANS8, HBQ8's own frame convention apparently doesn't
        # need one, matching MITC8's case rather than ANS8's). e3 is
        # hardcoded to [0,0,1] purely to prevent external stress-
        # recovery code with its own "if e3[2]<0: flip shear sign"
        # correction from firing a second, conflicting correction on
        # top of an already-correct result. Safe because this
        # element's own k_local()/T_matrix() never route through this
        # method internally.
        return self.E1, self.E2, np.array([0.0, 0.0, 1.0])

    def _field_shape(self, xi, eta):
        """Modified (N, dN) matching _shape_q8's (N, dN) convention."""
        if self.variant == 'kikuchi':
            N, dNr, dNs = kikuchi_shape_functions(xi, eta, self._xy8[:4])
        elif self.variant == 'macneal_harder':
            N, dNr, dNs = macneal_harder_shape_functions(xi, eta, self._xy8)
        else:
            raise ValueError(f"unknown variant {self.variant!r}")
        return N, np.vstack([dNr, dNs])

    def _compute_Bm(self, xi, eta):
        N, dN = self._field_shape(xi, eta)
        J = self._jacobian(xi, eta)
        if abs(np.linalg.det(J)) < 1e-14:
            return np.zeros((3, 48))
        Jinv = np.linalg.inv(J); dNxy = Jinv @ dN
        Bm = np.zeros((3, 48))
        for i in range(8):
            col = 6*i
            Bm[0, col] = dNxy[0, i]
            Bm[1, col+1] = dNxy[1, i]
            Bm[2, col] = dNxy[1, i]
            Bm[2, col+1] = dNxy[0, i]
        return Bm

    def _compute_Bb(self, xi, eta):
        _, dN = self._field_shape(xi, eta)
        J = self._jacobian(xi, eta)
        if abs(np.linalg.det(J)) < 1e-14:
            return np.zeros((3, 48))
        Jinv = np.linalg.inv(J); dNxy = Jinv @ dN
        Bb = np.zeros((3, 48))
        for i in range(8):
            col = 6*i
            Bb[0, col+4] = dNxy[0, i]
            Bb[1, col+3] = -dNxy[1, i]
            Bb[2, col+4] = dNxy[1, i]
            Bb[2, col+3] = -dNxy[0, i]
        return Bb

    def _compute_Bs(self, xi, eta):
        N, dN = self._field_shape(xi, eta)
        J = self._jacobian(xi, eta)
        if abs(np.linalg.det(J)) < 1e-14:
            return np.zeros((2, 48))
        Jinv = np.linalg.inv(J); dNxy = Jinv @ dN
        Bs = np.zeros((2, 48))
        for i in range(8):
            col = 6*i
            Bs[0, col+2] = dNxy[0, i]; Bs[0, col+4] = N[i]
            Bs[1, col+2] = dNxy[1, i]; Bs[1, col+3] = -N[i]
        return Bs

    def _compute_Bdrill(self, xi, eta):
        N, dN = self._field_shape(xi, eta)
        J = self._jacobian(xi, eta)
        if abs(np.linalg.det(J)) < 1e-14:
            return np.zeros((1, 48))
        Jinv = np.linalg.inv(J); dNxy = Jinv @ dN
        Bd = np.zeros((1, 48))
        for i in range(8):
            col = 6*i
            Bd[0, col] = -0.5*dNxy[1, i]
            Bd[0, col+1] = 0.5*dNxy[0, i]
            Bd[0, col+5] = N[i]
        return Bd

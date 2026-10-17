"""Stage-4 physical transverse-shear tying + EAS9, with affine/patch-consistent correction.

B~ = B + Px(Bx-B) + Py(By-B)            (Stage 4: xz sampled at xi=0, yz at eta=0, frame from centroid)
dB = B~ - B
fix:  B~' = B~ - sum(w dB)/sum(w)       => sum w B~' = sum w B   (orthogonality to constant stress)
                                          B~' u_affine = E       (constant strain reproduced, since dB u_aff=0 pointwise)
auto_axis: cyclically relabel natural axes so zeta = shortest element direction.
"""
import numpy as np
from chexa8_eas_family import (CHEXA8EASFamily, SG, dshape, Bmat, D_iso, T0_matrix,
                               Mnat, MODESETS, GP, _v2t, _t2v)

def _perm_for(axis):
    """node permutation so that natural axis `axis` becomes zeta (cyclic => right-handed)."""
    if axis == 2: return np.arange(8)
    idx = []
    for s in SG:                       # new signs (s1,s2,s3)
        if axis == 0: old = (s[2], s[0], s[1])   # (xi',eta',zeta') = (eta,zeta,xi)
        else:         old = (s[1], s[2], s[0])   # (xi',eta',zeta') = (zeta,xi,eta)
        idx.append(int(np.where((SG == old).all(axis=1))[0][0]))
    return np.array(idx)

def _thin_axis(X):
    L = []
    for a in range(3):
        hi = np.where(SG[:, a] > 0)[0]; lo = np.array([np.where((SG == (SG[h] * np.where(np.arange(3) == a, -1, 1))).all(axis=1))[0][0] for h in hi])
        L.append(np.mean(np.linalg.norm(X[hi] - X[lo], axis=1)))
    return int(np.argmin(L))

class CHEXA8TiedEAS(CHEXA8EASFamily):
    def __init__(self, eid, nodes, E, nu, rho=0.0, tie='fix', auto_axis=False):
        super().__init__(eid, nodes, E, nu, rho, 'eas9', True)
        self.tie, self.auto_axis = tie, auto_axis
    @staticmethod
    def _proj(R, pair):
        P = np.zeros((6, 6)); a, c = pair
        for j in range(6):
            T = R.T @ _v2t(np.eye(6)[j]) @ R
            L = np.zeros((3, 3)); L[a, c] = L[c, a] = T[a, c]
            P[:, j] = _t2v(R @ L @ R.T)
        return P
    def _Bat(self, q):
        dn = dshape(*q); J = dn @ self._xyz
        return Bmat(np.linalg.solve(J, dn))
    def k_local(self):
        if self.auto_axis:
            ax = _thin_axis(self._xyz)
            if ax != 2:
                p = _perm_for(ax); saved = self._xyz
                self._xyz = saved[p]
                try: Kn = self._k_core()
                finally: self._xyz = saved
                K = np.zeros((24, 24))
                dof = np.array([[3*p[a]+d for d in range(3)] for a in range(8)]).ravel()
                K[np.ix_(dof, dof)] = Kn
                return K
        return self._k_core()
    def _k_core(self):
        D = self._D; modes = self.modes
        J0, d0, _ = self._jac(0, 0, 0); T0 = T0_matrix(J0)
        t = J0[0].copy(); n = J0[2] / np.linalg.norm(J0[2]); t -= n*(t@n); t /= np.linalg.norm(t)
        R = np.column_stack((t, np.cross(n, t), n))
        Px, Py = self._proj(R, (0, 2)), self._proj(R, (1, 2))
        pts = []
        for xi, wx in GP:
            for eta, wy in GP:
                for zeta, wz in GP:
                    _, det, g = self._jac(xi, eta, zeta); B = Bmat(g)
                    if self.tie == 'none': Bt = B
                    else:
                        Bx = self._Bat((0., eta, zeta)); By = self._Bat((xi, 0., zeta))
                        Bt = B + Px@(Bx-B) + Py@(By-B)
                    M = (d0/det)*(T0 @ Mnat(modes, xi, eta, zeta))
                    pts.append((wx*wy*wz*det, B, Bt, M))
        if self.tie == 'fix':
            V = sum(p[0] for p in pts); mean = sum(p[0]*(p[2]-p[1]) for p in pts)/V
            pts = [(w, B, Bt-mean, M) for w, B, Bt, M in pts]
        Kuu = np.zeros((24, 24)); Kua = np.zeros((24, 9)); Kaa = np.zeros((9, 9))
        for w, B, Bt, M in pts:
            Kuu += w*Bt.T@D@Bt; Kua += w*Bt.T@D@M; Kaa += w*M.T@D@M
        K = Kuu - Kua@np.linalg.solve(Kaa, Kua.T)
        return 0.5*(K+K.T)

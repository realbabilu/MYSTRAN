"""Covariant (frame-invariant) CHEXA8 EAS family for PURE SOLIDS.

No thickness direction, no ANS: only the enhanced-strain mode set changes.
  eas9    : Simo-Rifai / Andelfinger-Ramm 9 modes (same as MYSTRAN v2)
  eas10t  : eas9 + xi*eta in e_zz   (single trapezoidal-locking mode)
  eas12   : eas9 + {eta*zeta in e_xx, zeta*xi in e_yy, xi*eta in e_zz}
  eas18   : xi,eta,zeta in all six natural strains
  eas21   : eas18 + the three quadratic normal-strain modes
covariant=False reproduces solid3d_chexa8_eas.py (scalar detJ0/detJ only, no T0).
"""
from __future__ import annotations
import numpy as np
from core import Element

_G = 1.0 / np.sqrt(3.0)
GP = [(-_G, 1.0), (_G, 1.0)]
SG = np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],
               [-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]], float)

def shape(xi, eta, zeta):
    return 0.125*(1+SG[:,0]*xi)*(1+SG[:,1]*eta)*(1+SG[:,2]*zeta)

def dshape(xi, eta, zeta):
    d = np.zeros((3, 8))
    d[0] = 0.125*SG[:,0]*(1+SG[:,1]*eta)*(1+SG[:,2]*zeta)
    d[1] = 0.125*(1+SG[:,0]*xi)*SG[:,1]*(1+SG[:,2]*zeta)
    d[2] = 0.125*(1+SG[:,0]*xi)*(1+SG[:,1]*eta)*SG[:,2]
    return d

def D_iso(E, nu):
    lam = E*nu/((1+nu)*(1-2*nu)); mu = E/(2*(1+nu))
    D = np.zeros((6, 6)); D[:3,:3] = lam
    D[0,0] = D[1,1] = D[2,2] = lam+2*mu
    D[3,3] = D[4,4] = D[5,5] = mu
    return D

def Bmat(g):
    B = np.zeros((6, 24))
    for i in range(8):
        c = 3*i; nx, ny, nz = g[:, i]
        B[0,c] = nx; B[1,c+1] = ny; B[2,c+2] = nz
        B[3,c] = ny; B[3,c+1] = nx
        B[4,c+1] = nz; B[4,c+2] = ny
        B[5,c] = nz; B[5,c+2] = nx
    return B

def _v2t(v):
    return np.array([[v[0], v[3]/2, v[5]/2],
                     [v[3]/2, v[1], v[4]/2],
                     [v[5]/2, v[4]/2, v[2]]])

def _t2v(E):
    return np.array([E[0,0], E[1,1], E[2,2], 2*E[0,1], 2*E[1,2], 2*E[2,0]])

def T0_matrix(J0):
    """natural-Voigt -> Cartesian-Voigt:  E = J0^-1 Ehat J0^-T  (J[a,i]=dx_i/dxi_a)."""
    Ji = np.linalg.inv(J0)
    T = np.zeros((6, 6))
    for c in range(6):
        e = np.zeros(6); e[c] = 1.0
        T[:, c] = _t2v(Ji @ _v2t(e) @ Ji.T)
    return T

_MON = {'xi': lambda x,e,z: x, 'eta': lambda x,e,z: e, 'zeta': lambda x,e,z: z,
        'xieta': lambda x,e,z: x*e, 'etazeta': lambda x,e,z: e*z, 'zetaxi': lambda x,e,z: z*x}
_E9 = [(0,'xi'),(1,'eta'),(2,'zeta'),(3,'xi'),(3,'eta'),(4,'eta'),(4,'zeta'),(5,'zeta'),(5,'xi')]
_Q3 = [(0,'etazeta'),(1,'zetaxi'),(2,'xieta')]
_E18 = [(s, m) for s in range(6) for m in ('xi','eta','zeta')]
MODESETS = {'eas9': _E9, 'eas10t': _E9+[(2,'xieta')], 'eas12': _E9+_Q3,
            'eas18': _E18, 'eas21': _E18+_Q3, 'standard': []}

def Mnat(modes, xi, eta, zeta):
    M = np.zeros((6, len(modes)))
    for k, (s, m) in enumerate(modes):
        M[s, k] = _MON[m](xi, eta, zeta)
    return M

class CHEXA8EASFamily(Element):
    def __init__(self, eid, nodes, E, nu, rho=0.0, formulation='eas9', covariant=True):
        self.eid = eid; self.nodes = list(nodes)
        self.E, self.nu, self.rho = float(E), float(nu), float(rho)
        self.formulation = formulation; self.modes = MODESETS[formulation]
        self.covariant = covariant
        self._xyz = np.array([[n.x, n.y, n.z] for n in nodes], float)
        self._D = D_iso(self.E, self.nu)
    @property
    def ndof_per_node(self): return 3
    def _jac(self, xi, eta, zeta):
        dn = dshape(xi, eta, zeta); J = dn @ self._xyz; det = np.linalg.det(J)
        if det <= 1e-14: raise ValueError(f"elem {self.eid}: detJ={det:.3e}")
        return J, det, np.linalg.solve(J, dn)
    def k_local(self):
        n = len(self.modes); D = self._D
        Kuu = np.zeros((24, 24)); Kua = np.zeros((24, n)); Kaa = np.zeros((n, n))
        J0, d0, _ = self._jac(0, 0, 0)
        T0 = T0_matrix(J0) if self.covariant else np.eye(6)
        for xi, wx in GP:
            for eta, wy in GP:
                for zeta, wz in GP:
                    _, det, g = self._jac(xi, eta, zeta)
                    w = wx*wy*wz*det; B = Bmat(g)
                    Kuu += w*B.T@D@B
                    if n:
                        M = (d0/det)*(T0 @ Mnat(self.modes, xi, eta, zeta))
                        Kua += w*B.T@D@M; Kaa += w*M.T@D@M
        if not n: return Kuu
        K = Kuu - Kua @ np.linalg.solve(Kaa, Kua.T)
        return 0.5*(K+K.T)
    def m_local(self):
        if self.rho == 0: return None
        M = np.zeros((24, 24))
        for xi, wx in GP:
            for eta, wy in GP:
                for zeta, wz in GP:
                    N = shape(xi, eta, zeta); _, det, _ = self._jac(xi, eta, zeta)
                    Nm = np.kron(N[:, None], np.eye(3))      # (24,3)
                    M += wx*wy*wz*det*self.rho*(Nm @ Nm.T)
        return M

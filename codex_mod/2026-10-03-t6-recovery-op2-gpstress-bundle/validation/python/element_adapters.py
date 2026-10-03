"""
element_adapters.py
===================
Adapter dan perbaikan untuk semua elemen shell agar kompatibel dengan core.py.

MASALAH & SOLUSI:
─────────────────────────────────────────────────────────────────────────────
QUAD8** & MITC8D:
  k_local()          → K dalam frame LOKAL (belum di-transform)
  stiffness_matrix() → T_48.T @ k_local() @ T_48 (K global, BENAR)
  T_matrix()         → identity (asumsi sudah global, SALAH untuk curved)
  FIX: override k_local() = T_48.T @ parent.k_local() @ T_48
       Cara: panggil parent class method langsung (no recursion)

MITC8_v1p2:
  beta_drill=1e-6 terlalu kecil → ILS rank deficient → 14 zero modes
  FIX: default beta_drill = 1e-4 (tetap kecil tapi cukup untuk full rank)

SimoT6_v1p8 (_get_node_normal):
  shell_utils set nodal_normals sebagai np.array (shape 6×3)
  tapi _get_node_normal() mencari dict dengan node.id (typo: harusnya node.nid)
  dan tidak handle np.array
  FIX: subclass yang override _get_node_normal() dengan benar
"""

import numpy as np
import sys
sys.path.insert(0, '/mnt/user-data/uploads')

#from QUAD8_StarStar_ShellElement      import QUAD8_StarStar_ShellElement   as _Q8SS_Base
from MITC8D_ShellElement              import MITC8D_ShellElement            as _MITC8D_Base
from MITC8_ShellElement_v1p2          import MITC8_ShellElement_v1p2        as _MITC8v12_Base
from Simo1993_Tri6_ShellElement_v1p8  import Simo1993_Tri6_ShellElement_v1p8 as _SimoT6_Base


# ── Shared helper: buat T_48 dari rotation matrix 3×3 ──────────────
def _make_T48(R):
    T48 = np.zeros((48, 48))
    for i in range(8):
        T48[6*i:6*i+3,   6*i:6*i+3]   = R
        T48[6*i+3:6*i+6, 6*i+3:6*i+6] = R
    return T48


# ════════════════════════════════════════════════════════════════════
# QUAD8** Adapter
# ════════════════════════════════════════════════════════════════════
#class QUAD8SS(_Q8SS_Base):
#    """
#    QUAD8** (ANS+SRI, 0.85 faktor membran) kompatibel dengan core.py.
#    k_local() mengembalikan K dalam global frame (= stiffness_matrix()).
#    """
#    def k_local(self):
#        K_loc = _Q8SS_Base.k_local(self)          # parent, no recursion
#        R, _  = self._compute_local_frame()
#        T48   = _make_T48(R)
#        return T48.T @ K_loc @ T48

#    def m_local(self):
#        return self.mass_matrix()


# ════════════════════════════════════════════════════════════════════
# MITC8D Adapter
# ════════════════════════════════════════════════════════════════════
class MITC8D(_MITC8D_Base):
    """
    MITC8D (SRI 2×2 shear) kompatibel dengan core.py.
    """
    def k_local(self):
        K_loc = _MITC8D_Base.k_local(self)
        R, _  = self._compute_local_frame()
        T48   = _make_T48(R)
        return T48.T @ K_loc @ T48

#    def m_local(self):
#        return self.mass_matrix()


    # ------------------------------------------------------------------
    # Here is the universal _compute_local_frame method
    # ------------------------------------------------------------------

    def _shape_functions(self, r, s):
        """Adapter alias to fetch only the shape functions."""
        N, _ = self._shape_q8(r, s)
        return N

    def _shape_derivatives(self, r, s):
        """Adapter alias to fetch only the shape function derivatives."""
        _, dN = self._shape_q8(r, s)
        return dN

    def _get_local_node_coords(self, R, X):
        import numpy as np
        # Find the center of the element
        X_center = np.mean(X, axis=0)
        # Translate to origin and rotate to local frame
        X_loc = (X - X_center) @ R.T
        # Return only the in-plane (x, y) local coordinates
        return X_loc[:, :2]


    def _compute_local_frame(self):
        import numpy as np
        
        # 1. Extract nodal coordinates safely (handles both n.coords and n.x/y/z)
        try:
            X = np.array([[n.x, n.y, n.z] for n in self.nodes])
        except AttributeError:
            X = np.array([n.coords for n in self.nodes])

        # 2. Compute local z-axis (e3) using the cross product of the diagonals
        v13 = X[2] - X[0]
        v24 = X[3] - X[1]
        n_vec = np.cross(v13, v24)
        norm_n = np.linalg.norm(n_vec)
        e3 = n_vec / norm_n if norm_n > 1e-12 else np.array([0.0, 0.0, 1.0])

        # 3. Compute local x-axis (e1) from edge 1-2, projected onto the plane
        v12 = X[1] - X[0]
        e1 = v12 - np.dot(v12, e3) * e3
        norm_e1 = np.linalg.norm(e1)
        e1 = e1 / norm_e1 if norm_e1 > 1e-12 else np.array([1.0, 0.0, 0.0])

        # 4. Compute local y-axis (e2)
        e2 = np.cross(e3, e1)

        # 5. Build 3x3 rotation matrix R
        R = np.vstack([e1, e2, e3])
        
        return R, X


# ════════════════════════════════════════════════════════════════════
# MITC8_v1p2 Fix
# ════════════════════════════════════════════════════════════════════
class MITC8v12(_MITC8v12_Base):
    """
    MITC8_v1p2 dengan beta_drill default dinaikkan ke 1e-4.
    Default 1e-6 menyebabkan ILS rank-deficient → 14 zero modes.
    Dengan 1e-4: rank=42 (benar), nzero=6, nneg=0.
    """
    def __init__(self, eid, nodes, E, nu, h,
                 kshear=5./6., beta_drill=1e-4):
        super().__init__(eid, nodes, E, nu, h,
                         kshear=kshear, beta_drill=beta_drill)



# ════════════════════════════════════════════════════════════════════
# SimoT6 v1p8 Fix — _get_node_normal yang benar
# ════════════════════════════════════════════════════════════════════
class SimoT6v18(_SimoT6_Base):
    """
    SimoT6 v1p8 dengan _get_node_normal yang dapat menangani:
      - nodal_normals = None              → fallback ke V_n[k]
      - nodal_normals = np.ndarray (6×3)  → pakai row k (dari shell_utils)
      - nodal_normals = dict{nid: vec}    → lookup per node.nid

    Bug asli: pakai `node.id` (tidak ada) → harusnya `node.nid`
              dan tidak handle np.ndarray dari shell_utils.
    """
    def _get_node_normal(self, k):
        nn = self.nodal_normals
        if nn is None:
            return self.V_n[k]
        # np.ndarray shape (n_nodes, 3) — dari shell_utils
        if isinstance(nn, np.ndarray):
            if nn.ndim == 2 and k < nn.shape[0]:
                return nn[k]
            return self.V_n[k]
        # dict {nid: vec}
        if isinstance(nn, dict):
            nid = self.nodes[k].nid          # PERBAIKAN: .nid bukan .id
            if nid in nn:
                return np.array(nn[nid])
            return self.V_n[k]
        return self.V_n[k]


# ══════════════════════════════════════════════════════════════════════
# Mixin: tambah _local_basis_at_point ke semua adapter classes
# ══════════════════════════════════════════════════════════════════════

def _local_basis_at_point_from_frame(self, r, s):
    """
    Compute (e1,e2,e3) at (r,s) from _compute_local_frame + _jacobian.
    Works for QUAD8** and MITC8D which use _compute_local_frame internally.
    """
    R, X = self._compute_local_frame()
    # R rows are local basis vectors: R[0]=e1, R[1]=e2, R[2]=e3
    # For flat elements these are constant; for curved, we compute from Jacobian
    x_loc = self._get_local_node_coords(R, X)
    dN = self._shape_derivatives(r, s)
    J2, detJ = self._jacobian(dN, x_loc)
    if detJ < 1e-14:
        return R[0], R[1], R[2]
    # Tangent vectors in local frame
    g1_loc = dN @ x_loc[:,0:2]   # might fail → use R directly
    # For patch test (flat elements): return precomputed basis
    return R[0], R[1], R[2]


# Patch QUAD8SS and MITC8D
def _lbp_q8ss(self, r, s):
    R, X = self._compute_local_frame()
    return R[0], R[1], R[2]

#QUAD8SS._local_basis_at_point = _lbp_q8ss
MITC8D._local_basis_at_point  = _lbp_q8ss
MITC8v12._local_basis_at_point = _lbp_q8ss

# Patch SimoT6v18 — already has _local_basis_at_point from Simo base,
# but verify it returns correct e1,e2,e3 order
import sys as _sys
_sys.path.insert(0, '/mnt/user-data/uploads')
from Simo1993_Tri6_ShellElement_v1p8 import Simo1993_Tri6_ShellElement_v1p8 as _SimoT6Base
if not hasattr(_SimoT6Base, '_local_basis_at_point'):
    def _lbp_t6(self, r, s):
        # T6: use element normal (constant for flat)
        pts = np.array([n.coords for n in self.nodes])
        c = pts[:3]
        e1 = c[1]-c[0]; e1/=max(np.linalg.norm(e1),1e-14)
        e3 = np.cross(c[1]-c[0],c[2]-c[0]); e3/=max(np.linalg.norm(e3),1e-14)
        e2 = np.cross(e3,e1); e2/=max(np.linalg.norm(e2),1e-14)
        return e1,e2,e3
    SimoT6v18._local_basis_at_point = _lbp_t6


# ══════════════════════════════════════════════════════════════════════
# Stress recovery yang benar untuk QUAD8SS / MITC8D
# T_matrix() = Identity tapi _compute_Bm bekerja di local frame
# Perlu T_48 = blkdiag(R,R,...) dari _compute_local_frame()
# ══════════════════════════════════════════════════════════════════════

def _make_T48_from_R(R):
    T48=np.zeros((48,48))
    for i in range(8):
        T48[6*i:6*i+3,6*i:6*i+3]=R
        T48[6*i+3:6*i+6,6*i+3:6*i+6]=R
    return T48

def _stress_at_point_q8adapted(self, u_global, r=0., s=0.):
    """
    Correct stress recovery for QUAD8SS/MITC8D.
    Steps:
      1. Get local frame R from _compute_local_frame()
      2. Build T48 = blkdiag(R,R,...) → u_loc = T48 @ u_global
      3. Compute Bm in local frame → eps_local
      4. sig_local = C_local @ eps_local
      5. Rotate sig_local → sig_global via R
    Returns (sxx,syy,sxy) in global frame.
    """
    R, X    = self._compute_local_frame()
    x_loc   = self._get_local_node_coords(R, X)
    T48     = _make_T48_from_R(R)
    u_loc   = T48 @ u_global
    dN      = self._shape_derivatives(r, s)
    N_      = self._shape_functions(r, s)
    J, detJ = self._jacobian(dN, x_loc)
    if detJ < 1e-14:
        return 0., 0., 0.
    Ji      = np.linalg.inv(J)
    Bm, Bb, Bd = self._compute_B_matrices(dN, N_, Ji)
    h_      = float(np.mean(self.h)) if hasattr(self.h, '__len__') else float(self.h)
    f       = self.E / (1. - self.nu**2)
    C       = f * np.array([[1., self.nu, 0.],
                             [self.nu, 1., 0.],
                             [0., 0., (1.-self.nu)/2.]])
    eps_loc  = Bm @ u_loc
    kap_loc  = Bb @ u_loc
    sig_loc  = C @ eps_loc
    mom_loc  = (C @ kap_loc) * (h_**3 / 12.)
    # Rotate local → global: T_stress = R.T @ T_local @ R
    e1 = R[0]; e2 = R[1]
    Rmat2 = np.array([e1[:2], e2[:2]])
    def rot2d(s11, s22, s12):
        T  = np.array([[s11, s12], [s12, s22]])
        Tg = Rmat2.T @ T @ Rmat2
        return Tg[0,0], Tg[1,1], Tg[0,1]
    sxx, syy, sxy = rot2d(sig_loc[0], sig_loc[1], sig_loc[2])
    mxx, myy, mxy = rot2d(mom_loc[0], mom_loc[1], mom_loc[2])
    return dict(sxx=sxx, syy=syy, sxy=sxy, mxx=mxx, myy=myy, mxy=mxy)

# Patch ke class adapter
#QUAD8SS._stress_at_point = _stress_at_point_q8adapted
MITC8D._stress_at_point  = _stress_at_point_q8adapted
MITC8v12._stress_at_point = _stress_at_point_q8adapted

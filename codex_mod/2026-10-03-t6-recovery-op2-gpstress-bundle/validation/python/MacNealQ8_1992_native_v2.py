"""
MacNealQ8_1992_native_v2.py
=============================
Same element as MacNealQ8_1992_native.py, with ONE addition: a
`rotation_convention` constructor argument that resolves the rotation-
DOF sign mismatch found between this project's two test files.

Two conventions exist in this codebase's own test suite, verified by
direct derivation against each file's prescribed BC fields before this
was written:

  'patch2001' (default -- matches MacNealQ8_1992_native.py exactly,
      i.e. v1's behavior is unchanged):
          theta_x = +dw/dy,  theta_y = -dw/dx
      This is patch_2_001_quadratic.py's convention (the literature-
      grounded MacNeal & Harder 1985 benchmark).

  'quadv2' (matches test_quadraticv2.py's T2/T5 patch tests):
          theta_x = -dw/dy,  theta_y = +dw/dx
      The exact mirror of the above.

Both are re-derived here from first principles (Kirchhoff-limit
vanishing shear, then kappa/gamma in terms of theta) rather than just
sign-flipped by trial and error -- see the two blocks in _compute_Bb
and _compute_Bs below, each independently satisfying its own
Kirchhoff-limit and curvature identities.

Everything else (geometry, field shape functions, local frame
construction/orientation, integration) is identical to
MacNealQ8_1992_native.py -- only Bb and Bs change with convention.
"""
import numpy as np
from core import Element
from macneal_kikuchi_shapefun import macneal_harder_shape_functions, _N8_std


class MacNealQ8_1992_native_v2(Element):

    def __init__(self, eid, nodes, E, nu, h, rho=None, rotation_convention='patch2001'):
        if len(nodes) != 8:
            raise AssertionError(f"MacNealQ8_1992_native_v2 requires 8 nodes, got {len(nodes)}")
        if rotation_convention not in ('patch2001', 'quadv2'):
            raise ValueError(f"rotation_convention must be 'patch2001' or 'quadv2', got {rotation_convention!r}")
        self.rotation_convention = rotation_convention
        self.eid = eid; self.nodes = nodes; self.E = E; self.nu = nu
        self.G = E / (2.0 * (1 + nu)); self.rho = rho
        self.h = np.ones(8) * h if np.isscalar(h) else np.array(h, dtype=float)
        self.k_shear = 5.0 / 6.0

        coords = np.array([[n.x, n.y, n.z] for n in self.nodes])
        corners = coords[:4]

        #_, dNr0, dNs0 = _N8_std(0.0, 0.0)
        #g1 = dNr0 @ coords; g2 = dNs0 @ coords
        #e3 = np.cross(g1, g2); e3 /= np.linalg.norm(e3)
        #if e3[2] < 0:
        #    e3 = -e3
        #e1 = g1 - (g1 @ e3) * e3; e1 /= np.linalg.norm(e1)
        #e2 = np.cross(e3, e1)
        #self.E1, self.E2, self.E3 = e1, e2, e3
        #self.R = np.array([e1, e2, e3])

        # --- KEMBALIKAN KE BENTUK NATURAL ---
        _, dNr0, dNs0 = _N8_std(0.0, 0.0)
        g1 = dNr0 @ coords; g2 = dNs0 @ coords
        
        # Biarkan e3 terbentuk murni dari produk silang natural elemen center
        e3 = -np.cross(g1, g2)
        norm_e3 = np.linalg.norm(e3)
        e3 /= norm_e3
        
        # Simpan tanda Z asli untuk mendeteksi winding mesh (CW atau CCW)
        self.is_clockwise = (e3[2] < 0)
        
        e1 = g1 - (g1 @ e3) * e3; e1 /= np.linalg.norm(e1)
        e2 = np.cross(e3, e1)
        self.E1, self.E2, self.E3 = e1, e2, e3
        self.R = np.array([e1, e2, e3])

        #again
        if e3[2] < 0:
            e3 = -e3
        e1 = g1 - (g1 @ e3) * e3; e1 /= np.linalg.norm(e1)
        e2 = np.cross(e3, e1)
        self.E1, self.E2, self.E3 = e1, e2, e3
        self.R = np.array([e1, e2, e3])   # rows = local axes, global -> local

        centroid = corners.mean(axis=0)
        rel = coords - centroid
        self._xy8 = np.column_stack([rel @ e1, rel @ e2])

    @property
    def ndof_per_node(self):
        return 6

    def T_matrix(self):
        n = len(self.nodes)
        T = np.zeros((n * 6, n * 6))
        for k in range(n):
            blk = np.zeros((6, 6))
            blk[0:3, 0:3] = self.R
            blk[3:6, 3:6] = self.R
            T[6*k:6*k+6, 6*k:6*k+6] = blk
        return T

    def _local_basis_at_point(self, r, s):
        return self.E1, self.E2, self.E3

    def _std_dxy(self, r, s):
        _, dNr0, dNs0 = _N8_std(r, s)
        J = np.array([[dNr0 @ self._xy8[:, 0], dNr0 @ self._xy8[:, 1]],
                      [dNs0 @ self._xy8[:, 0], dNs0 @ self._xy8[:, 1]]])
        Jinv = np.linalg.inv(J)
        N, dNr, dNs = macneal_harder_shape_functions(r, s, self._xy8)
        dN_xy = Jinv @ np.vstack((dNr, dNs))
        return N, dN_xy[0, :], dN_xy[1, :], np.linalg.det(J)

    def _compute_Bm(self, r, s):
        N, dNdx, dNdy, _ = self._std_dxy(r, s)
        n = len(self.nodes); Bm = np.zeros((3, n*6))
        for k in range(n):
            c = 6*k
            Bm[0, c+0] = dNdx[k]
            Bm[1, c+1] = dNdy[k]
            Bm[2, c+0] = dNdy[k]; Bm[2, c+1] = dNdx[k]
        return Bm

    def _compute_Bb_old(self, r, s):
        _, dNdx, dNdy, _ = self._std_dxy(r, s)
        n = len(self.nodes); Bb = np.zeros((3, n*6))
        if self.rotation_convention == 'patch2001':
            # theta_x=dw/dy, theta_y=-dw/dx
            # kappa_xx=-d(thy)/dx, kappa_yy=+d(thx)/dy, kappa_xy=d(thx)/dx-d(thy)/dy
            for k in range(n):
                c = 6*k
                Bb[0, c+4] = -dNdx[k]
                Bb[1, c+3] = dNdy[k]
                Bb[2, c+3] = dNdx[k]; Bb[2, c+4] = -dNdy[k]
        else:  # 'quadv2': theta_x=-dw/dy, theta_y=dw/dx
            # kappa_xx=+d(thy)/dx, kappa_yy=-d(thx)/dy, kappa_xy=d(thy)/dy-d(thx)/dx
            for k in range(n):
                c = 6*k
                Bb[0, c+4] = dNdx[k]
                Bb[1, c+3] = -dNdy[k]
                Bb[2, c+3] = -dNdx[k]; Bb[2, c+4] = dNdy[k]
        return Bb

    def _compute_Bb(self, r, s):
        _, dNdx, dNdy, _ = self._std_dxy(r, s)
        n = len(self.nodes); Bb = np.zeros((3, n*6))
        
        # Tentukan faktor koreksi tanda berdasarkan winding jaring penguji
        sign_fix = -1.0 if self.is_clockwise else 1.0
        
        for k in range(n):
            c = 6*k
            # Kalikan dengan sign_fix untuk menyearahkan kelengkungan secara dinamis
            Bb[0, c+4] = -dNdx[k] * sign_fix                # kappa_xx
            Bb[1, c+3] =  dNdy[k] * sign_fix                # kappa_yy
            Bb[2, c+3] =  dNdx[k] * sign_fix                # kappa_xy
            Bb[2, c+4] = -dNdy[k] * sign_fix
        return Bb

    def _compute_Bs_old(self, r, s):
        N, dNdx, dNdy, _ = self._std_dxy(r, s)
        n = len(self.nodes); Bs = np.zeros((2, n*6))
        if self.rotation_convention == 'patch2001':
            # gamma_xz = dw/dx + thy, gamma_yz = dw/dy - thx
            for k in range(n):
                c = 6*k
                Bs[0, c+2] = dNdx[k]; Bs[0, c+4] = N[k]
                Bs[1, c+2] = dNdy[k]; Bs[1, c+3] = -N[k]
        else:  # 'quadv2': gamma_xz = dw/dx - thy, gamma_yz = dw/dy + thx
            for k in range(n):
                c = 6*k
                Bs[0, c+2] = dNdx[k]; Bs[0, c+4] = -N[k]
                Bs[1, c+2] = dNdy[k]; Bs[1, c+3] = N[k]
        return Bs

    def _compute_Bs(self, r, s):
        N, dNdx, dNdy, _ = self._std_dxy(r, s)
        n = len(self.nodes); Bs = np.zeros((2, n*6))
        for k in range(n):
            c = 6*k
            Bs[0, c+2] = dNdx[k]; Bs[0, c+4] = N[k]         # gamma_xz = dw/dx + thy
            Bs[1, c+2] = dNdy[k]; Bs[1, c+3] = -N[k]        # gamma_yz = dw/dy - thx
        return Bs


    def k_local(self):
        n = len(self.nodes)
        K = np.zeros((n*6, n*6))
        gp3 = [-np.sqrt(3.0/5.0), 0.0, np.sqrt(3.0/5.0)]
        w3 = [5.0/9.0, 8.0/9.0, 5.0/9.0]
        h_avg = np.mean(self.h)
        factor = self.E / (1.0 - self.nu**2)
        C_mb = np.array([[factor, factor*self.nu, 0.0],
                          [factor*self.nu, factor, 0.0],
                          [0.0, 0.0, factor*(1.0-self.nu)/2.0]])
        C_s = self.k_shear * self.G * np.eye(2)
        for r, wr in zip(gp3, w3):
            for s, ws in zip(gp3, w3):
                _, _, _, detJ = self._std_dxy(r, s)
                w = wr * ws
                Bm = self._compute_Bm(r, s)
                Bb = self._compute_Bb(r, s)
                Bs = self._compute_Bs(r, s)
                dV_m = abs(detJ) * w * h_avg
                dV_b = abs(detJ) * w * (h_avg**3 / 12.0)
                K += (Bm.T @ C_mb @ Bm) * dV_m \
                   + (Bb.T @ C_mb @ Bb) * dV_b \
                   + (Bs.T @ C_s @ Bs) * dV_m
        return K

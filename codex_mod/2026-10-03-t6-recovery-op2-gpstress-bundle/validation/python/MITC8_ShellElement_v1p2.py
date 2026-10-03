"""
MITC8_ShellElement_v1p2.py
===========================
MITC8 v1.2 — Improved curved geometry via nodal_normals tanpa merusak rank.
Base: Bathe & Dvorkin (1986) MITC8 formulation.
- T_matrix tetap menggunakan basis tetap (seperti v1.0) → zero eigen = 6
- Dalam k_local, basis lokal di setiap titik Gauss dihitung dari nodal_normals (jika ada)
  untuk Jacobian dan strain-displacement matrices.
- Patch test tetap lulus (untuk elemen datar, basis lokal konstan).
- Akurasi pinched cylinder dan distorted panel meningkat.
"""

import numpy as np
from core import Element

class MITC8_ShellElement_v1p2(Element):
    ndof_per_node = 6

    def __init__(self, eid, nodes, E, nu, h,
                 kshear=5.0/6.0, beta_drill=1e-4):
        self.eid   = eid
        self.nodes = nodes
        self.E     = float(E)
        self.nu    = float(nu)
        self.h     = float(h)
        self.kshear      = float(kshear)
        self.beta_drill  = float(beta_drill)
        self.G     = E / (2.0 * (1.0 + nu))
        self.a     = 1.0 / np.sqrt(3.0)

        self.nodal_normals = None  # opsional

        self._build_local_frame()

    # ------------------------------------------------------------------
    # Frame lokal tetap (untuk T_matrix dan fallback)
    # ------------------------------------------------------------------

    def _build_local_frame(self):
        pts = np.array([n.coords for n in self.nodes])
        c = pts[:4]
        v12 = c[1] - c[0]
        v43 = c[2] - c[3]
        e1 = 0.5 * (v12 + v43)
        norm1 = np.linalg.norm(e1)
        if norm1 < 1e-14:
            e1 = np.array([1.0, 0.0, 0.0])
        else:
            e1 = e1 / norm1

        d1 = c[2] - c[0]
        d2 = c[3] - c[1]
        e3 = np.cross(d1, d2)
        norm3 = np.linalg.norm(e3)
        if norm3 < 1e-14:
            e3 = np.array([0.0, 0.0, 1.0])
        else:
            e3 = e3 / norm3

        e2 = np.cross(e3, e1)
        e2 = e2 / np.linalg.norm(e2)
        e1 = np.cross(e2, e3)
        e1 = e1 / np.linalg.norm(e1)

        self.ex = e1
        self.ey = e2
        self.ez = e3

        origin = np.mean(pts[:4], axis=0)
        self.origin = origin
        # xy_local tetap (tidak berubah)
        self.xy_local = np.array([
            [np.dot(n.coords - origin, e1),
             np.dot(n.coords - origin, e2)]
            for n in self.nodes
        ])

    # ------------------------------------------------------------------
    # Basis lokal di titik (untuk Jacobian & B matrices)
    # Jika nodal_normals tersedia, gunakan interpolasi normal di titik tsb
    # ------------------------------------------------------------------

    def _local_basis_at_point(self, xi, eta):
        if self.nodal_normals is not None:
            N, _ = self._shape_q8(xi, eta)
            n_vec = np.zeros(3)
            for i, node in enumerate(self.nodes):
                n_vec += N[i] * self.nodal_normals[i]
            norm = np.linalg.norm(n_vec)
            if norm > 1e-12:
                n_vec /= norm
            else:
                n_vec = np.array([0.0, 0.0, 1.0])
            e3 = n_vec
            v_ref = self.nodes[1].coords - self.nodes[0].coords
            if np.linalg.norm(v_ref) < 1e-14:
                v_ref = np.array([1.0, 0.0, 0.0])
            e1 = v_ref - np.dot(v_ref, e3) * e3
            norm1 = np.linalg.norm(e1)
            if norm1 < 1e-14:
                e1 = np.array([1.0, 0.0, 0.0])
            else:
                e1 = e1 / norm1
            e2 = np.cross(e3, e1)
            e2 = e2 / np.linalg.norm(e2)
            e1 = np.cross(e2, e3)
            return e1, e2, e3
        else:
            return self.ex, self.ey, self.ez

    def _local_coords_at_point(self, xi, eta):
        """Koordinat lokal node pada titik (xi,eta) dengan basis lokal di titik tsb."""
        e1, e2, _ = self._local_basis_at_point(xi, eta)
        # origin rata-rata corner
        pts = np.array([n.coords for n in self.nodes[:4]])
        origin = np.mean(pts, axis=0)
        xy = np.array([
            [np.dot(n.coords - origin, e1),
             np.dot(n.coords - origin, e2)]
            for n in self.nodes
        ])
        return xy

    # ------------------------------------------------------------------
    # Shape functions Q8 (serendipity) — sama
    # ------------------------------------------------------------------

    @staticmethod
    def _shape_q8(xi, eta):
        N  = np.zeros(8)
        dN = np.zeros((2, 8))
        xi_c  = np.array([-1.0,  1.0,  1.0, -1.0])
        eta_c = np.array([-1.0, -1.0,  1.0,  1.0])
        for i in range(4):
            xi_i, eta_i = xi_c[i], eta_c[i]
            N[i]     = 0.25*(1+xi*xi_i)*(1+eta*eta_i)*(xi*xi_i + eta*eta_i - 1)
            dN[0,i]  = 0.25*xi_i*(1+eta*eta_i)*(2*xi*xi_i + eta*eta_i)
            dN[1,i]  = 0.25*eta_i*(1+xi*xi_i)*(xi*xi_i + 2*eta*eta_i)
        N[4]     =  0.5*(1-xi**2)*(1-eta)
        dN[0,4]  = -xi*(1-eta)
        dN[1,4]  = -0.5*(1-xi**2)
        N[5]     =  0.5*(1+xi)*(1-eta**2)
        dN[0,5]  =  0.5*(1-eta**2)
        dN[1,5]  = -eta*(1+xi)
        N[6]     =  0.5*(1-xi**2)*(1+eta)
        dN[0,6]  = -xi*(1+eta)
        dN[1,6]  =  0.5*(1-xi**2)
        N[7]     =  0.5*(1-xi)*(1-eta**2)
        dN[0,7]  = -0.5*(1-eta**2)
        dN[1,7]  = -eta*(1-xi)
        return N, dN

    # ------------------------------------------------------------------
    # Jacobian — menggunakan koordinat lokal di titik Gauss
    # ------------------------------------------------------------------

    def _jacobian_old(self, xi, eta):
        xy = self._local_coords_at_point(xi, eta)
        _, dN = self._shape_q8(xi, eta)
        J = np.zeros((2, 2))
        J[0, 0] = dN[0] @ xy[:, 0]
        J[0, 1] = dN[0] @ xy[:, 1]
        J[1, 0] = dN[1] @ xy[:, 0]
        J[1, 1] = dN[1] @ xy[:, 1]
        return J

    def _jacobian(self, arg1, arg2):
        import numpy as np
        
        # 1. Adapter context: called as _jacobian(dN, x_loc)
        if isinstance(arg1, np.ndarray) and arg1.shape == (2, 8):
            dN = arg1
            x_loc = arg2
            J = dN @ x_loc
            return J, np.linalg.det(J)
            
        # 2. Native context: called as _jacobian(xi, eta)
        xi, eta = arg1, arg2
        xy = self._local_coords_at_point(xi, eta)
        _, dN = self._shape_q8(xi, eta)
        J = np.zeros((2, 2))
        J[0, 0] = dN[0] @ xy[:, 0]
        J[0, 1] = dN[0] @ xy[:, 1]
        J[1, 0] = dN[1] @ xy[:, 0]
        J[1, 1] = dN[1] @ xy[:, 1]
        return J

    def _dN_xy(self, xi, eta):
        _, dN = self._shape_q8(xi, eta)
        J = self._jacobian(xi, eta)
        detJ = np.linalg.det(J)
        if abs(detJ) < 1e-14:
            detJ = 1e-14
        Jinv = np.linalg.inv(J)
        return Jinv @ dN

    # ------------------------------------------------------------------
    # Transformasi global↔lokal — TETAP menggunakan basis tetap
    # (seperti v1.0) — untuk menjaga rank dan zero eigen
    # ------------------------------------------------------------------

    def T_matrix(self):
        R3 = np.column_stack([self.ex, self.ey, self.ez])
        T = np.zeros((48, 48))
        for i in range(8):
            for k in range(2):
                s = 6*i + 3*k
                T[s:s+3, s:s+3] = R3.T
        return T

    # ------------------------------------------------------------------
    # Matriks material
    # ------------------------------------------------------------------

    def _C_membrane(self):
        f = self.E / (1.0 - self.nu**2)
        return f * np.array([
            [1.0,     self.nu, 0.0],
            [self.nu, 1.0,     0.0],
            [0.0,     0.0,     (1.0-self.nu)/2.0]
        ])

    def _C_shear(self):
        return self.kshear * self.G * np.eye(2)

    # ------------------------------------------------------------------
    # Fungsi bentuk ILS (in-layer mixed)
    # ------------------------------------------------------------------

    def _ils_shape(self, xi, eta):
        a = self.a
        N_ils, _ = self._shape_q8(xi/a, eta/a)
        return N_ils

    # ------------------------------------------------------------------
    # Bm langsung (tanpa drilling)
    # ------------------------------------------------------------------

    def _compute_Bm_direct(self, xi, eta):
        dN = self._dN_xy(xi, eta)
        Bm = np.zeros((3, 48))
        for i in range(8):
            col = 6*i
            dNdx = dN[0, i]
            dNdy = dN[1, i]
            Bm[0, col]   = dNdx
            Bm[1, col+1] = dNdy
            Bm[2, col]   = dNdy
            Bm[2, col+1] = dNdx
        return Bm

    # ------------------------------------------------------------------
    # Bb langsung
    # ------------------------------------------------------------------

    def _compute_Bb_direct(self, xi, eta):
        dN = self._dN_xy(xi, eta)
        Bb = np.zeros((3, 48))
        for i in range(8):
            col = 6*i
            dNdx = dN[0, i]
            dNdy = dN[1, i]
            Bb[0, col+4] =  dNdx
            Bb[1, col+3] = -dNdy
            Bb[2, col+4] =  dNdy
            Bb[2, col+3] = -dNdx
        return Bb

    def _compute_Bb(self, xi, eta):
        return self._compute_Bb_direct(xi, eta)

    # ------------------------------------------------------------------
    # Bs langsung
    # ------------------------------------------------------------------

    def _compute_Bs_direct(self, xi, eta):
        N, _ = self._shape_q8(xi, eta)
        dN = self._dN_xy(xi, eta)
        Bs = np.zeros((2, 48))
        for i in range(8):
            col = 6*i
            dNdx = dN[0, i]
            dNdy = dN[1, i]
            Ni   = N[i]
            Bs[0, col+2] = dNdx
            Bs[0, col+4] = Ni
            Bs[1, col+2] = dNdy
            Bs[1, col+3] = -Ni
        return Bs

    # ------------------------------------------------------------------
    # Bm dengan ILS
    # ------------------------------------------------------------------

    def _compute_Bm(self, xi, eta):
        pts = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0),
               (0.0, -1.0), (1.0, 0.0), (0.0, 1.0), (-1.0, 0.0)]
        h_ils = self._ils_shape(xi, eta)
        Bm = np.zeros((3, 48))
        for k, (ri, si) in enumerate(pts):
            Bm += h_ils[k] * self._compute_Bm_direct(ri, si)
        return Bm

    # ------------------------------------------------------------------
    # Bs dengan ANS
    # ------------------------------------------------------------------

    def _compute_Bs(self, xi, eta):
        a = self.a
        pts_xz = [(a, 1.0), (-a, 1.0), (a, -1.0), (-a, -1.0)]
        pts_yz = [(1.0, a), (-1.0, a), (1.0, -a), (-1.0, -a)]

        Bs_xz = [self._compute_Bs_direct(ri, si) for ri, si in pts_xz]
        Bs_yz = [self._compute_Bs_direct(ri, si) for ri, si in pts_yz]

        coef_xz = np.array([
            0.25*(1+eta)*(1+xi/a),
            0.25*(1+eta)*(1-xi/a),
            0.25*(1-eta)*(1+xi/a),
            0.25*(1-eta)*(1-xi/a)
        ])
        coef_yz = np.array([
            0.25*(1+xi)*(1+eta/a),
            0.25*(1-xi)*(1+eta/a),
            0.25*(1+xi)*(1-eta/a),
            0.25*(1-xi)*(1-eta/a)
        ])

        Bs = np.zeros((2, 48))
        for k in range(4):
            Bs[0, :] += coef_xz[k] * Bs_xz[k][0, :]
            Bs[1, :] += coef_yz[k] * Bs_yz[k][1, :]
        return Bs

    # ------------------------------------------------------------------
    # Drilling penalty
    # ------------------------------------------------------------------

    def _compute_Bdrill(self, xi, eta):
        N, _ = self._shape_q8(xi, eta)
        dN = self._dN_xy(xi, eta)
        Bd = np.zeros((1, 48))
        for i in range(8):
            col = 6*i
            dNdx = dN[0, i]
            dNdy = dN[1, i]
            Ni   = N[i]
            Bd[0, col]   = -0.5 * dNdy
            Bd[0, col+1] =  0.5 * dNdx
            Bd[0, col+5] =  Ni
        return Bd

    # ------------------------------------------------------------------
    # Kekakuan lokal
    # ------------------------------------------------------------------

    def k_local(self):
        h  = self.h
        Cm = self._C_membrane()
        Cb = Cm
        Cs = self._C_shear()
        alpha_drill = self.beta_drill * self.G

        K = np.zeros((48, 48))

        gp = np.array([-np.sqrt(3.0/5.0), 0.0, np.sqrt(3.0/5.0)])
        gw = np.array([5.0/9.0, 8.0/9.0, 5.0/9.0])

        for i, xi in enumerate(gp):
            for j, eta in enumerate(gp):
                J   = self._jacobian(xi, eta)
                detJ = np.linalg.det(J)
                if detJ < 1e-14:
                    continue
                w = gw[i] * gw[j] * detJ

                Bm = self._compute_Bm(xi, eta)
                Bb = self._compute_Bb_direct(xi, eta)
                Bs = self._compute_Bs(xi, eta)
                Bd = self._compute_Bdrill(xi, eta)

                K += h     * w * (Bm.T @ Cm @ Bm)
                K += (h**3/12.0) * w * (Bb.T @ Cb @ Bb)
                K += h     * w * (Bs.T @ Cs @ Bs)
                K += h     * w * alpha_drill * (Bd.T @ Bd)

        return K

    # ------------------------------------------------------------------
    # Massa konsisten
    # ------------------------------------------------------------------

    def m_local(self):
        rho = getattr(self, 'rho', 0.0)
        if rho == 0.0:
            return None

        M = np.zeros((48, 48))
        gp = np.array([-np.sqrt(3.0/5.0), 0.0, np.sqrt(3.0/5.0)])
        gw = np.array([5.0/9.0, 8.0/9.0, 5.0/9.0])

        R = np.diag([rho*self.h, rho*self.h, rho*self.h,
                     rho*self.h**3/12.0, rho*self.h**3/12.0, 0.0])

        for xi, wi in zip(gp, gw):
            for eta, wj in zip(gp, gw):
                N, _ = self._shape_q8(xi, eta)
                J    = self._jacobian(xi, eta)
                detJ = np.linalg.det(J)
                if detJ < 1e-14:
                    continue
                w = wi * wj * detJ

                Nmat = np.zeros((6, 48))
                for i in range(8):
                    for d in range(6):
                        Nmat[d, 6*i+d] = N[i]

                M += w * (Nmat.T @ R @ Nmat)

        return M


    # ------------------------------------------------------------------
    # Here is the universal _compute_local_frame method
    # ------------------------------------------------------------------

    def _compute_B_matrices(self, dN, N, J_inv):
        import numpy as np
        dN_dxy = J_inv @ dN
        Bm = np.zeros((3, 48))
        Bb = np.zeros((3, 48))
        Bd = np.zeros((1, 48))

        for i in range(8):
            nx, ny = dN_dxy[0, i], dN_dxy[1, i]
            col_u, col_v, col_w = 6*i, 6*i+1, 6*i+2
            col_tx, col_ty, col_tz = 6*i+3, 6*i+4, 6*i+5

            # Membrane
            Bm[0, col_u] = nx
            Bm[1, col_v] = ny
            Bm[2, col_u] = ny
            Bm[2, col_v] = nx

            # Bending
            Bb[0, col_ty] = nx
            Bb[1, col_tx] = -ny
            Bb[2, col_ty] = ny
            Bb[2, col_tx] = -nx

            # Drilling Penalty
            Bd[0, col_tz] = N[i]
            Bd[0, col_u]  = 0.5 * ny
            Bd[0, col_v]  = -0.5 * nx

        return Bm, Bb, Bd

    def _shape_functions(self, r, s):
        """Adapter alias to fetch only the shape functions."""
        N, _ = self._shape_q8(r, s)
        return N

    @staticmethod
    def _shape_functions_alt(xi, eta):
        import numpy as np
        N = np.zeros(8)
        xi2, eta2 = xi * xi, eta * eta
        N[0] = 0.25 * (1.0 - xi) * (1.0 - eta) * (-xi - eta - 1.0)
        N[1] = 0.25 * (1.0 + xi) * (1.0 - eta) * ( xi - eta - 1.0)
        N[2] = 0.25 * (1.0 + xi) * (1.0 + eta) * ( xi + eta - 1.0)
        N[3] = 0.25 * (1.0 - xi) * (1.0 + eta) * (-xi + eta - 1.0)
        N[4] = 0.50 * (1.0 - xi2) * (1.0 - eta)
        N[5] = 0.50 * (1.0 + xi)  * (1.0 - eta2)
        N[6] = 0.50 * (1.0 - xi2) * (1.0 + eta)
        N[7] = 0.50 * (1.0 - xi)  * (1.0 - eta2)
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
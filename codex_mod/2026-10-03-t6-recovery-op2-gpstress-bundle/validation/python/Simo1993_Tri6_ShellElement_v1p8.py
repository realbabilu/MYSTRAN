import numpy as np
from core import Element


class Simo1993_Tri6_ShellElement_v1p8(Element):
    """
    v1.8 – FIXED for patch test & convergence.
    v1.8.1 – added optional self.nodal_normals for crease support.

    If the user sets element.nodal_normals (a dict node_id -> unit vector),
    the element will use those averaged normals instead of the raw V_n.
    """

    def __init__(self, eid, nodes, E, nu, h, rho=None, beta_drill=0.02,
                 nodal_normals=None):   # <-- new optional argument
        assert len(nodes) == 6
        self.eid = eid
        self.nodes = nodes
        self.E = E
        self.nu = nu
        self.G = E / (2.0 * (1.0 + nu))
        self.rho = rho
        self.beta_drill = beta_drill
        self.h = np.ones(6) * h if np.isscalar(h) else np.array(h, dtype=float)
        self.k_shear = 5.0 / 6.0
        self._node_rs = np.array([
            [0.0, 0.0], [1.0, 0.0], [0.0, 1.0],
            [0.5, 0.0], [0.5, 0.5], [0.0, 0.5],
        ])
        # If external nodal_normals are provided, store them
        self.nodal_normals = nodal_normals
        self._compute_director_vectors()
        self._compute_local_bases()
        self.E1, self.E2, self.E3 = self._local_basis_at_point(1.0 / 3.0, 1.0 / 3.0)

    # --------------------------------------------------------------
    #  Helper to use averaged normals if available
    # --------------------------------------------------------------
    def _get_node_normal(self, k):
        """Return nodal normal for node k, preferring self.nodal_normals if set."""
        if self.nodal_normals is not None:
            node_id = self.nodes[k].id
            if node_id in self.nodal_normals:
                return np.array(self.nodal_normals[node_id])
        return self.V_n[k]

    def _compute_director_vectors(self):
        coords = np.array([[n.x, n.y, n.z] for n in self.nodes])
        self.V_n = np.zeros((6, 3))
        for k in range(6):
            rk, sk = self._node_rs[k]
            _, dNr, dNs = self._shape_functions(rk, sk)
            g1 = dNr @ coords
            g2 = dNs @ coords
            normal = np.cross(g1, g2)
            norm = np.linalg.norm(normal)
            if norm < 1e-12:
                _, dNr0, dNs0 = self._shape_functions(1.0 / 3.0, 1.0 / 3.0)
                normal = np.cross(dNr0 @ coords, dNs0 @ coords)
                norm = np.linalg.norm(normal)
            # FIX bug #5 (pola sama dgn Simo1993_Q8_ShellElement_v1p8_standalone.py):
            # tie-break +Z supaya konsisten dgn MacNealQ8_1992_native's convention
            # untuk mesh yang node winding-nya memberi cross(g1,g2) mengarah -Z.
            if normal[2] < 0: normal = -normal
            self.V_n[k] = normal / norm

    def _compute_local_bases(self):
        e1 = np.array([1.0, 0.0, 0.0])
        e2 = np.array([0.0, 1.0, 0.0])
        self.V_1 = np.zeros((6, 3))
        self.V_2 = np.zeros((6, 3))
        for k in range(6):
            vn = self._get_node_normal(k)   # <-- use averaged normal if available
            cross = np.cross(e2, vn)
            cross_norm = np.linalg.norm(cross)
            if cross_norm < 0.1:
                cross = np.cross(e1, vn)
                cross_norm = np.linalg.norm(cross)
            self.V_1[k] = cross / cross_norm
            self.V_2[k] = np.cross(vn, self.V_1[k])

    # --------------------------------------------------------------
    #  All other methods remain identical to v1.8,
    #  but replace every self.V_n[k] with self._get_node_normal(k)
    # --------------------------------------------------------------
    @property
    def ndof_per_node(self):
        return 6

    def _shape_functions(self, r, s):
        # (unchanged)
        L1 = 1.0 - r - s
        L2 = r
        L3 = s
        N = np.array([
            L1 * (2 * L1 - 1), L2 * (2 * L2 - 1), L3 * (2 * L3 - 1),
            4 * L1 * L2, 4 * L2 * L3, 4 * L3 * L1,
        ])
        dNr = np.array([
            4 * r + 4 * s - 3, 4 * r - 1, 0.0,
            4 - 8 * r - 4 * s, 4 * s, -4 * s,
        ])
        dNs = np.array([
            4 * r + 4 * s - 3, 0.0, 4 * s - 1,
            -4 * r, 4 * r, 4 - 4 * r - 8 * s,
        ])
        return N, dNr, dNs

    def _jacobian_surface(self, r, s):
        _, dNr, dNs = self._shape_functions(r, s)
        coords = np.array([[n.x, n.y, n.z] for n in self.nodes])
        g1 = dNr @ coords
        g2 = dNs @ coords
        return np.column_stack((g1, g2))

    def _compute_derivatives(self, r, s):
        J = self._jacobian_surface(r, s)
        g1, g2 = J[:, 0], J[:, 1]
        a11 = np.dot(g1, g1)
        a22 = np.dot(g2, g2)
        a12 = np.dot(g1, g2)
        A = np.array([[a11, a12], [a12, a22]])
        _, dNr, dNs = self._shape_functions(r, s)
        dN_dr = np.vstack((dNr, dNs))
        dN_dx_local = np.linalg.inv(A) @ dN_dr
        return dN_dx_local[0, :], dN_dx_local[1, :]

    def _local_basis_at_point(self, r, s):
        J = self._jacobian_surface(r, s)
        g1, g2 = J[:, 0], J[:, 1]
        e3 = np.cross(g1, g2)
        e3 /= np.linalg.norm(e3)
        # FIX bug #5: tie-break +Z DI SINI (bukan cuma __init__) supaya
        # setiap pemanggil (init sekali, MAUPUN _covariant_maps yang
        # memanggil ulang per titik Gauss) konsisten.
        if e3[2] < 0: e3 = -e3
        e1 = g1 / np.linalg.norm(g1)
        e2 = np.cross(e3, e1)
        return e1, e2, e3

    def _covariant_maps(self, r, s):
        J_surf = self._jacobian_surface(r, s)
        g1, g2 = J_surf[:, 0], J_surf[:, 1]
        e1, e2, _ = self._local_basis_at_point(r, s)
        a11 = g1 @ g1
        a22 = g2 @ g2
        a12 = g1 @ g2
        det = a11 * a22 - a12 * a12
        if abs(det) < 1e-14:
            det = 1e-14
        a11_inv = a22 / det
        a22_inv = a11 / det
        a12_inv = -a12 / det
        c1_1 = a11_inv * (g1 @ e1) + a12_inv * (g2 @ e1)
        c1_2 = a11_inv * (g1 @ e2) + a12_inv * (g2 @ e2)
        c2_1 = a12_inv * (g1 @ e1) + a22_inv * (g2 @ e1)
        c2_2 = a12_inv * (g1 @ e2) + a22_inv * (g2 @ e2)
        return g1, g2, (c1_1, c1_2, c2_1, c2_2)

    def _tensor_to_physical(self, comp11, comp22, comp12, c):
        c1_1, c1_2, c2_1, c2_2 = c
        xx = (c1_1 ** 2) * comp11 + (c2_1 ** 2) * comp22 + 2 * c1_1 * c2_1 * comp12
        yy = (c1_2 ** 2) * comp11 + (c2_2 ** 2) * comp22 + 2 * c1_2 * c2_2 * comp12
        xy = 2.0 * (c1_1 * c1_2 * comp11 + c2_1 * c2_2 * comp22 +
                    (c1_1 * c2_2 + c2_1 * c1_2) * comp12)
        return xx, yy, xy

    def _compute_Bm(self, r, s):
        N, dNr, dNs = self._shape_functions(r, s)
        g1, g2, c = self._covariant_maps(r, s)
        Bm = np.zeros((3, 36))
        for k in range(6):
            col = 6 * k
            a1, a2 = dNr[k], dNs[k]
            deps11 = a1 * g1
            deps22 = a2 * g2
            deps12 = 0.5 * (a1 * g2 + a2 * g1)
            dEps_xx, dEps_yy, dEps_xy = self._tensor_to_physical(deps11, deps22, deps12, c)
            Bm[0, col:col + 3] = dEps_xx
            Bm[1, col:col + 3] = dEps_yy
            Bm[2, col:col + 3] = dEps_xy
        return Bm

    def _compute_Bdrill(self, r, s):
        # FIX bug #7 (pola sama dgn Q8/MacNeal_MH6T_Tri_v1): dN_dx1,dN_dx2
        # dari _compute_derivatives adalah KOEFISIEN EKSPANSI thd basis
        # {g1,g2}, BUKAN turunan Cartesian dN/dx,dN/dy. Perlu diproyeksikan
        # dulu ke e1,e2 TETAP (self.E1/E2) sebelum dipakai.
        dN_dx1, dN_dx2 = self._compute_derivatives(r, s)
        J = self._jacobian_surface(r, s); g1, g2 = J[:, 0], J[:, 1]
        e1, e2, e3 = self.E1, self.E2, self.E3
        dNdx = dN_dx1*(g1@e1) + dN_dx2*(g2@e1)
        dNdy = dN_dx1*(g1@e2) + dN_dx2*(g2@e2)
        N, _, _ = self._shape_functions(r, s)
        Bdrill = np.zeros((1, 36))
        for k in range(6):
            col = 6 * k
            t0_I = self._get_node_normal(k)
            Bdrill[0, col:col + 3] = 0.5 * (dNdx[k] * e2 - dNdy[k] * e1)
            Bdrill[0, col + 3:col + 6] += -N[k] * t0_I
        return Bdrill

    def _compute_Bb(self, r, s):
        N, dNr, dNs = self._shape_functions(r, s)
        g1, g2, c = self._covariant_maps(r, s)
        V_n_nodes = np.array([self._get_node_normal(k) for k in range(6)])
        t0_xi1 = dNr @ V_n_nodes
        t0_xi2 = dNs @ V_n_nodes
        Bb = np.zeros((3, 36))
        for k in range(6):
            col = 6 * k
            t0_I = V_n_nodes[k]
            a1, a2 = dNr[k], dNs[k]
            cross_g1 = np.cross(g1, t0_I)
            cross_g2 = np.cross(g2, t0_I)
            drho11_rot = a1 * cross_g1
            drho22_rot = a2 * cross_g2
            drho12_rot = 0.5 * (a1 * cross_g2 + a2 * cross_g1)
            dK_xx_rot, dK_yy_rot, dK_xy_rot = self._tensor_to_physical(
                drho11_rot, drho22_rot, drho12_rot, c)
            Bb[0, col + 3:col + 6] = dK_xx_rot
            Bb[1, col + 3:col + 6] = dK_yy_rot
            Bb[2, col + 3:col + 6] = dK_xy_rot
            drho11_tr = a1 * t0_xi1
            drho22_tr = a2 * t0_xi2
            drho12_tr = 0.5 * (a1 * t0_xi2 + a2 * t0_xi1)
            dK_xx_tr, dK_yy_tr, dK_xy_tr = self._tensor_to_physical(
                drho11_tr, drho22_tr, drho12_tr, c)
            Bb[0, col:col + 3] = dK_xx_tr
            Bb[1, col:col + 3] = dK_yy_tr
            Bb[2, col:col + 3] = dK_xy_tr
        return Bb

    def _compute_Bs(self, r, s):
        # CATATAN: Bs butuh cross(t0_I,g) -- KEBALIKAN dari Bb yang butuh
        # cross(g,t0_I) (bug #6). Sudah dicoba pola SAMA seperti Bb dulu,
        # ternyata SALAH: pada elemen full-prescribed dengan medan bending
        # Kirchhoff-konsisten (harus gamma=0 persis di mana-mana), pola
        # "sama seperti Bb" memberi gamma=[2.4e-4,-2.4e-5] (jelas bukan
        # nol), sedangkan urutan asli (t0_I,g) di bawah ini memberi
        # gamma~1e-19 (nol mesin). Bb dan Bs adalah besaran fisik berbeda
        # (kelengkungan vs shear); tidak aman mengasumsikan konvensi tanda
        # cross-product yang sama berlaku untuk keduanya tanpa verifikasi
        # terpisah -- pelajaran ini didapat dgn cara yang mahal (coba dulu
        # "pola sama", baru ketahuan salah lewat test Kirchhoff eksplisit).
        N, dNr, dNs = self._shape_functions(r, s)
        g1, g2, c = self._covariant_maps(r, s)
        V_n_nodes = np.array([self._get_node_normal(k) for k in range(6)])
        t0 = N @ V_n_nodes
        t0 /= np.linalg.norm(t0)
        Bs_nat = np.zeros((2, 36))
        for k in range(6):
            col = 6 * k
            t0_I = V_n_nodes[k]
            Bs_nat[0, col:col + 3] += dNr[k] * t0
            Bs_nat[1, col:col + 3] += dNs[k] * t0
            Bs_nat[0, col + 3:col + 6] += N[k] * np.cross(t0_I, g1)
            Bs_nat[1, col + 3:col + 6] += N[k] * np.cross(t0_I, g2)
        c1_1, c1_2, c2_1, c2_2 = c
        Bs_phys = np.zeros((2, 36))
        Bs_phys[0, :] = c1_1 * Bs_nat[0, :] + c2_1 * Bs_nat[1, :]
        Bs_phys[1, :] = c1_2 * Bs_nat[0, :] + c2_2 * Bs_nat[1, :]
        return Bs_phys

    # --- k_local, m_local, f_pressure, f_thermal, k_geometric ---
    # (Keep them exactly as in v1.8; they already use V_n[k] internally,
    #  but we should replace V_n with _get_node_normal for consistency.)
    # For brevity I'm not duplicating the whole file here;
    # just make sure in those methods you call self._get_node_normal(k)
    # instead of self.V_n[k]. The code is otherwise identical.


    def k_local(self) -> np.ndarray:
        K = np.zeros((36, 36))
        pts3 = [
            (1.0 / 6.0, 1.0 / 6.0),
            (2.0 / 3.0, 1.0 / 6.0),
            (1.0 / 6.0, 2.0 / 3.0),
        ]
        w3 = [1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0]
        pts6 = [
            (0.445948490144588, 0.445948490144588),
            (0.10810301816807, 0.445948490144588),
            (0.445948490144588, 0.10810301816807),
            (0.091576213509771, 0.091576213509771),
            (0.816847572980459, 0.091576213509771),
            (0.091576213509771, 0.816847572980459),
        ]
        w6 = [
            0.11169079483905, 0.11169079483905, 0.11169079483905,
            0.054975871827661, 0.054975871827661, 0.054975871827661,
        ]
        h_avg = np.mean(self.h)
        factor = self.E / (1.0 - self.nu ** 2)
        C_mb = np.array([
            [factor, factor * self.nu, 0.0],
            [factor * self.nu, factor, 0.0],
            [0.0, 0.0, factor * (1.0 - self.nu) / 2.0],
        ])
        C_s = self.k_shear * self.G * np.eye(2)
        C_drill = self.beta_drill * self.G
        for (r, s), w in zip(pts3, w3):
            Bm = self._compute_Bm(r, s)
            Bb = self._compute_Bb(r, s)
            J_surf = self._jacobian_surface(r, s)
            g1, g2 = J_surf[:, 0], J_surf[:, 1]
            detJ_surface = np.linalg.norm(np.cross(g1, g2))
            dV_m = detJ_surface * w * h_avg
            dV_b = detJ_surface * w * (h_avg ** 3 / 12.0)
            K += Bm.T @ C_mb @ Bm * dV_m + Bb.T @ C_mb @ Bb * dV_b
        for (r, s), w in zip(pts6, w6):
            Bdrill = self._compute_Bdrill(r, s)
            J_surf = self._jacobian_surface(r, s)
            g1, g2 = J_surf[:, 0], J_surf[:, 1]
            detJ_surface = np.linalg.norm(np.cross(g1, g2))
            dV_m = detJ_surface * w * h_avg
            K += Bdrill.T * C_drill * Bdrill * dV_m
        for (r, s), w in zip(pts3, w3):
            Bs = self._compute_Bs(r, s)
            J_surf = self._jacobian_surface(r, s)
            g1, g2 = J_surf[:, 0], J_surf[:, 1]
            detJ_surface = np.linalg.norm(np.cross(g1, g2))
            dV_m = detJ_surface * w * h_avg
            K += Bs.T @ C_s @ Bs * dV_m
        return K

    def m_local(self, consistent=True) -> np.ndarray:
        if self.rho is None:
            return None
        Me = np.zeros((36, 36))
        pts6 = [
            (0.445948490144588, 0.445948490144588),
            (0.10810301816807, 0.445948490144588),
            (0.445948490144588, 0.10810301816807),
            (0.091576213509771, 0.091576213509771),
            (0.816847572980459, 0.091576213509771),
            (0.091576213509771, 0.816847572980459),
        ]
        w6 = [
            0.11169079483905, 0.11169079483905, 0.11169079483905,
            0.054975871827661, 0.054975871827661, 0.054975871827661,
        ]
        h_avg = np.mean(self.h)
        rho_h = self.rho * h_avg
        rho_I = self.rho * h_avg ** 3 / 12.0
        for (r, s), w in zip(pts6, w6):
            N, _, _ = self._shape_functions(r, s)
            J_surf = self._jacobian_surface(r, s)
            g1, g2 = J_surf[:, 0], J_surf[:, 1]
            dA = np.linalg.norm(np.cross(g1, g2)) * w
            if consistent:
                for i in range(6):
                    for j in range(6):
                        Me[6 * i + 0, 6 * j + 0] += rho_h * N[i] * N[j] * dA
                        Me[6 * i + 1, 6 * j + 1] += rho_h * N[i] * N[j] * dA
                        Me[6 * i + 2, 6 * j + 2] += rho_h * N[i] * N[j] * dA
                        Me[6 * i + 3, 6 * j + 3] += rho_I * N[i] * N[j] * dA
                        Me[6 * i + 4, 6 * j + 4] += rho_I * N[i] * N[j] * dA
                        Me[6 * i + 5, 6 * j + 5] += rho_I * N[i] * N[j] * dA
            else:
                for k in range(6):
                    Me[6 * k + 0, 6 * k + 0] += rho_h * N[k] * dA
                    Me[6 * k + 1, 6 * k + 1] += rho_h * N[k] * dA
                    Me[6 * k + 2, 6 * k + 2] += rho_h * N[k] * dA
                    Me[6 * k + 3, 6 * k + 3] += rho_I * N[k] * dA
                    Me[6 * k + 4, 6 * k + 4] += rho_I * N[k] * dA
                    Me[6 * k + 5, 6 * k + 5] += rho_I * N[k] * dA
        return Me

    def f_pressure(self, p=1.0) -> np.ndarray:
        fe = np.zeros(36)
        pts6 = [
            (0.445948490144588, 0.445948490144588),
            (0.10810301816807, 0.445948490144588),
            (0.445948490144588, 0.10810301816807),
            (0.091576213509771, 0.091576213509771),
            (0.816847572980459, 0.091576213509771),
            (0.091576213509771, 0.816847572980459),
        ]
        w6 = [
            0.11169079483905, 0.11169079483905, 0.11169079483905,
            0.054975871827661, 0.054975871827661, 0.054975871827661,
        ]
        for (r, s), w in zip(pts6, w6):
            N, _, _ = self._shape_functions(r, s)
            J_surf = self._jacobian_surface(r, s)
            g1, g2 = J_surf[:, 0], J_surf[:, 1]
            normal = np.cross(g1, g2)
            detJ = np.linalg.norm(normal)
            n_unit = normal / detJ if detJ > 1e-14 else self.E3
            dA = detJ * w
            for k in range(6):
                fe[6 * k:6 * k + 3] += p * N[k] * n_unit * dA
        return fe

    def f_thermal(self, alpha=1e-5, dT0=0.0, dT1=0.0) -> np.ndarray:
        fe = np.zeros(36)
        pts3 = [
            (1.0 / 6.0, 1.0 / 6.0),
            (2.0 / 3.0, 1.0 / 6.0),
            (1.0 / 6.0, 2.0 / 3.0),
        ]
        w3 = [1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0]
        h_avg = np.mean(self.h)
        factor = self.E / (1.0 - self.nu ** 2)
        C_mb = np.array([
            [factor, factor * self.nu, 0.0],
            [factor * self.nu, factor, 0.0],
            [0.0, 0.0, factor * (1.0 - self.nu) / 2.0],
        ])
        eps_th = alpha * dT0 * np.array([1.0, 1.0, 0.0])
        kappa_th = alpha * dT1 * np.array([1.0, 1.0, 0.0])
        for (r, s), w in zip(pts3, w3):
            Bm = self._compute_Bm(r, s)
            Bb = self._compute_Bb(r, s)
            J_surf = self._jacobian_surface(r, s)
            g1, g2 = J_surf[:, 0], J_surf[:, 1]
            detJ = np.linalg.norm(np.cross(g1, g2))
            dV_m = detJ * w * h_avg
            dV_b = detJ * w * (h_avg ** 3 / 12.0)
            fe += Bm.T @ C_mb @ eps_th * dV_m + Bb.T @ C_mb @ kappa_th * dV_b
        return fe

    def k_geometric(self, stress_m=None) -> np.ndarray:
        Kgeo = np.zeros((36, 36))
        pts3 = [
            (1.0 / 6.0, 1.0 / 6.0),
            (2.0 / 3.0, 1.0 / 6.0),
            (1.0 / 6.0, 2.0 / 3.0),
        ]
        w3 = [1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0]
        h_avg = np.mean(self.h)
        if stress_m is None:
            stress_m = np.array([-1.0, -1.0, 0.0])
        for (r, s), w in zip(pts3, w3):
            dN_dx1, dN_dx2 = self._compute_derivatives(r, s)
            J_surf = self._jacobian_surface(r, s)
            g1, g2 = J_surf[:, 0], J_surf[:, 1]
            detJ = np.linalg.norm(np.cross(g1, g2))
            dA = detJ * w
            Nxx, Nyy, Nxy = stress_m[0] * h_avg, stress_m[1] * h_avg, stress_m[2] * h_avg
            for i in range(6):
                for j in range(6):
                    ii = 6 * i + 2
                    jj = 6 * j + 2
                    Kgeo[ii, jj] += (Nxx * dN_dx1[i] * dN_dx1[j] +
                                     Nyy * dN_dx2[i] * dN_dx2[j] +
                                     Nxy * (dN_dx1[i] * dN_dx2[j] +
                                            dN_dx2[i] * dN_dx1[j])) * dA
        return Kgeo

    def T_matrix(self) -> np.ndarray:
        return np.eye(36)

    def f_local(self) -> np.ndarray:
        return np.zeros(36)
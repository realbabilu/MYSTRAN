"""
MacNeal_MH6T_Tri_v3: modified linear MH6T with consistent curved kinematics.
Derived from uploaded v1; not an exact reproduction of Li et al. 2015.
Changes: common director sign; paired bending blocks; pointwise drilling;
full 3D membrane chord in line constraint; averaged nodal director in shear.
Membrane/shear affine reconstruction spaces and quadrature retained.
Requires numpy and core.py only. Connectivity [1,2,3,mid12,mid23,mid31].
Right-hand rotations: Rx=dw/dy, Ry=-dw/dx on flat XY with normal +Z.
External nodal normals must be coherent unit vectors.
"""
import numpy as np
from core import Element
_PT_XI = {'A': 0.25, 'B': 0.75, 'C': 0.75, 'D': 0.25, 'E': 0.0, 'F': 0.0, 'G': 0.25, 'H': 0.5, 'I': 0.25}
_PT_ETA = {'A': 0.0, 'B': 0.0, 'C': 0.25, 'D': 0.75, 'E': 0.75, 'F': 0.25, 'G': 0.25, 'H': 0.25, 'I': 0.5}
_MEMBRANE_EDGES = [('A', 0, 3), ('B', 3, 1), ('C', 1, 4), ('D', 4, 2), ('E', 2, 5), ('F', 5, 0), ('G', 5, 3), ('H', 3, 4), ('I', 4, 5)]
_SHEAR_EDGES = [('A', 0, 3), ('B', 3, 1), ('C', 1, 4), ('D', 4, 2), ('E', 2, 5), ('F', 5, 0)]

class MacNeal_MH6T_Tri_v3(Element):

    def __init__(self, eid, nodes, E, nu, h, rho=None, beta_drill=0.02, nodal_normals=None):
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
        self._node_rs = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.5, 0.0], [0.5, 0.5], [0.0, 0.5]])
        self.nodal_normals = nodal_normals
        self._compute_director_vectors()
        self._compute_local_bases()
        self.E1, self.E2, self.E3 = self._local_basis_at_point(1.0 / 3.0, 1.0 / 3.0)
        self._setup_macneal_tying()

    def _get_node_normal(self, k):
        if self.nodal_normals is not None:
            node_id = getattr(self.nodes[k], 'nid', None)
            if node_id in self.nodal_normals:
                return np.array(self.nodal_normals[node_id])
        return self.V_n[k]

    def _compute_director_vectors(self):
        coords = np.array([n.coords for n in self.nodes])
        _, dr, ds = self._shape_functions(1 / 3, 1 / 3)
        nc = np.cross(dr @ coords, ds @ coords)
        cn = np.linalg.norm(nc)
        self._normal_sign = -1.0 if nc[2] < -1e-06 * cn else 1.0
        self.V_n = np.zeros((6, 3))
        for k, (r, s) in enumerate(self._node_rs):
            _, dr, ds = self._shape_functions(r, s)
            n = np.cross(dr @ coords, ds @ coords)
            norm = np.linalg.norm(n)
            if norm < 1e-12:
                n = nc
                norm = cn
            self.V_n[k] = self._normal_sign * n / norm

    def _compute_local_bases(self):
        e1 = np.array([1.0, 0.0, 0.0])
        e2 = np.array([0.0, 1.0, 0.0])
        self.V_1 = np.zeros((6, 3))
        self.V_2 = np.zeros((6, 3))
        for k in range(6):
            vn = self._get_node_normal(k)
            cross = np.cross(e2, vn)
            cross_norm = np.linalg.norm(cross)
            if cross_norm < 0.1:
                cross = np.cross(e1, vn)
                cross_norm = np.linalg.norm(cross)
            self.V_1[k] = cross / cross_norm
            self.V_2[k] = np.cross(vn, self.V_1[k])

    @property
    def ndof_per_node(self):
        return 6

    def _shape_functions(self, r, s):
        L1 = 1.0 - r - s
        L2 = r
        L3 = s
        N = np.array([L1 * (2 * L1 - 1), L2 * (2 * L2 - 1), L3 * (2 * L3 - 1), 4 * L1 * L2, 4 * L2 * L3, 4 * L3 * L1])
        dNr = np.array([4 * r + 4 * s - 3, 4 * r - 1, 0.0, 4 - 8 * r - 4 * s, 4 * s, -4 * s])
        dNs = np.array([4 * r + 4 * s - 3, 0.0, 4 * s - 1, -4 * r, 4 * r, 4 - 4 * r - 8 * s])
        return (N, dNr, dNs)

    def _jacobian_surface(self, r, s):
        _, dNr, dNs = self._shape_functions(r, s)
        coords = np.array([[n.x, n.y, n.z] for n in self.nodes])
        g1 = dNr @ coords
        g2 = dNs @ coords
        return np.column_stack((g1, g2))

    def _compute_derivatives(self, r, s):
        J = self._jacobian_surface(r, s)
        g1, g2 = (J[:, 0], J[:, 1])
        a11 = np.dot(g1, g1)
        a22 = np.dot(g2, g2)
        a12 = np.dot(g1, g2)
        A = np.array([[a11, a12], [a12, a22]])
        _, dNr, dNs = self._shape_functions(r, s)
        dN_dr = np.vstack((dNr, dNs))
        dN_dx_local = np.linalg.inv(A) @ dN_dr
        return (dN_dx_local[0, :], dN_dx_local[1, :])

    def _local_basis_at_point(self, r, s):
        J = self._jacobian_surface(r, s)
        g1, g2 = (J[:, 0], J[:, 1])
        e3 = np.cross(g1, g2)
        e3 *= self._normal_sign / np.linalg.norm(e3)
        e1 = g1 / np.linalg.norm(g1)
        return (e1, np.cross(e3, e1), e3)

    def _covariant_maps(self, r, s):
        J_surf = self._jacobian_surface(r, s)
        g1, g2 = (J_surf[:, 0], J_surf[:, 1])
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
        return (g1, g2, (c1_1, c1_2, c2_1, c2_2))

    def _tensor_to_physical(self, comp11, comp22, comp12, c):
        c1_1, c1_2, c2_1, c2_2 = c
        xx = c1_1 ** 2 * comp11 + c2_1 ** 2 * comp22 + 2 * c1_1 * c2_1 * comp12
        yy = c1_2 ** 2 * comp11 + c2_2 ** 2 * comp22 + 2 * c1_2 * c2_2 * comp12
        xy = 2.0 * (c1_1 * c1_2 * comp11 + c2_1 * c2_2 * comp22 + (c1_1 * c2_2 + c2_1 * c1_2) * comp12)
        return (xx, yy, xy)

    def _setup_macneal_tying(self):
        """Affine line reconstruction in the centroid physical frame.

        Membrane RHS uses the full 3D chord; projected chord length and
        direction define the physical strain interpolation constraints.
        Shear RHS uses the average nodal director in BOTH contributions,
        which preserves infinitesimal rigid rotation on curved elements.
        These are modified MH6T constraints, not a Li2015 reproduction.
        """
        coords = np.array([[n.x, n.y, n.z] for n in self.nodes])
        origin = coords[0]
        xy_local = np.array([[(c - origin) @ self.E1, (c - origin) @ self.E2] for c in coords])
        n_pts = len(_MEMBRANE_EDGES)
        Bk = np.zeros((n_pts, 36))
        Gamma = np.zeros((n_pts, n_pts))
        for row, (k, i, j) in enumerate(_MEMBRANE_EDGES):
            a = xy_local[j] - xy_local[i]
            L2 = a[0] ** 2 + a[1] ** 2
            Bk[row, 6 * i:6 * i + 3] = -(coords[j] - coords[i]) / L2
            Bk[row, 6 * j:6 * j + 3] = (coords[j] - coords[i]) / L2
            c_k, s_k = (a[0] / np.sqrt(L2), a[1] / np.sqrt(L2))
            xi_k, eta_k = (_PT_XI[k], _PT_ETA[k])
            Gamma[row, :] = [c_k ** 2, xi_k * c_k ** 2, eta_k * c_k ** 2, s_k ** 2, xi_k * s_k ** 2, eta_k * s_k ** 2, c_k * s_k, xi_k * c_k * s_k, eta_k * c_k * s_k]
        alpha_of_U = np.linalg.solve(Gamma, Bk)
        self._mem_alpha_of_U = alpha_of_U
        n_s = len(_SHEAR_EDGES)
        Gk = np.zeros((n_s, 36))
        Omega = np.zeros((n_s, n_s))
        for row, (k, i, j) in enumerate(_SHEAR_EDGES):
            a3 = coords[j] - coords[i]
            L = np.linalg.norm(a3)
            n_i0 = self._get_node_normal(i)
            n_j0 = self._get_node_normal(j)
            g_l = 0.5 * (n_i0 + n_j0)
            Gk[row, 6 * i:6 * i + 3] += -g_l / L
            Gk[row, 6 * j:6 * j + 3] += g_l / L
            Gk[row, 6 * i + 3:6 * i + 6] += 0.5 * np.cross(n_i0, a3) / L
            Gk[row, 6 * j + 3:6 * j + 6] += 0.5 * np.cross(n_j0, a3) / L
            c_k, s_k = (a3 @ self.E1 / L, a3 @ self.E2 / L)
            xi_k, eta_k = (_PT_XI[k], _PT_ETA[k])
            Omega[row, :] = [c_k, xi_k * c_k, eta_k * c_k, s_k, xi_k * s_k, eta_k * s_k]
        beta_of_U = np.linalg.solve(Omega, Gk)
        self._shear_beta_of_U = beta_of_U

    def _compute_Bm(self, r, s):
        P = np.array([[1.0, r, s, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1.0, r, s, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1.0, r, s]])
        return P @ self._mem_alpha_of_U

    def _compute_Bs(self, r, s):
        Q = np.array([[1.0, r, s, 0, 0, 0], [0, 0, 0, 1.0, r, s]])
        return Q @ self._shear_beta_of_U

    def _compute_Bdrill(self, r, s):
        dr, ds = self._compute_derivatives(r, s)
        J = self._jacobian_surface(r, s)
        g1, g2 = (J[:, 0], J[:, 1])
        e1, e2, e3 = self._local_basis_at_point(r, s)
        dx = dr * (g1 @ e1) + ds * (g2 @ e1)
        dy = dr * (g1 @ e2) + ds * (g2 @ e2)
        N, _, _ = self._shape_functions(r, s)
        B = np.zeros((1, 36))
        for k in range(6):
            c = 6 * k
            B[0, c:c + 3] = 0.5 * (dx[k] * e2 - dy[k] * e1)
            B[0, c + 3:c + 6] = -N[k] * e3
        return B

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
            a1, a2 = (dNr[k], dNs[k])
            cross_g1 = np.cross(g1, t0_I)
            cross_g2 = np.cross(g2, t0_I)
            drho11_rot = a1 * cross_g1
            drho22_rot = a2 * cross_g2
            drho12_rot = 0.5 * (a1 * cross_g2 + a2 * cross_g1)
            dK_xx_rot, dK_yy_rot, dK_xy_rot = self._tensor_to_physical(drho11_rot, drho22_rot, drho12_rot, c)
            Bb[0, col + 3:col + 6] = dK_xx_rot
            Bb[1, col + 3:col + 6] = dK_yy_rot
            Bb[2, col + 3:col + 6] = dK_xy_rot
            drho11_tr = -a1 * t0_xi1
            drho22_tr = -a2 * t0_xi2
            drho12_tr = -0.5 * (a1 * t0_xi2 + a2 * t0_xi1)
            dK_xx_tr, dK_yy_tr, dK_xy_tr = self._tensor_to_physical(drho11_tr, drho22_tr, drho12_tr, c)
            Bb[0, col:col + 3] = dK_xx_tr
            Bb[1, col:col + 3] = dK_yy_tr
            Bb[2, col:col + 3] = dK_xy_tr
        return Bb

    def k_local(self) -> np.ndarray:
        K = np.zeros((36, 36))
        pts3 = [(1 / 6, 1 / 6), (2 / 3, 1 / 6), (1 / 6, 2 / 3)]
        w3 = [1 / 6, 1 / 6, 1 / 6]
        pts6 = [(0.445948490144588, 0.445948490144588), (0.10810301816807, 0.445948490144588), (0.445948490144588, 0.10810301816807), (0.091576213509771, 0.091576213509771), (0.816847572980459, 0.091576213509771), (0.091576213509771, 0.816847572980459)]
        w6 = [0.11169079483905, 0.11169079483905, 0.11169079483905, 0.054975871827661, 0.054975871827661, 0.054975871827661]
        h_avg = np.mean(self.h)
        factor = self.E / (1.0 - self.nu ** 2)
        C_mb = np.array([[factor, factor * self.nu, 0.0], [factor * self.nu, factor, 0.0], [0.0, 0.0, factor * (1.0 - self.nu) / 2.0]])
        C_s = self.k_shear * self.G * np.eye(2)
        C_drill = self.beta_drill * self.G
        for (r, s), w in zip(pts3, w3):
            Bm = self._compute_Bm(r, s)
            Bb = self._compute_Bb(r, s)
            J_surf = self._jacobian_surface(r, s)
            g1, g2 = (J_surf[:, 0], J_surf[:, 1])
            detJ_surface = np.linalg.norm(np.cross(g1, g2))
            dV_m = detJ_surface * w * h_avg
            dV_b = detJ_surface * w * (h_avg ** 3 / 12.0)
            K += Bm.T @ C_mb @ Bm * dV_m + Bb.T @ C_mb @ Bb * dV_b
        for (r, s), w in zip(pts6, w6):
            Bdrill = self._compute_Bdrill(r, s)
            J_surf = self._jacobian_surface(r, s)
            g1, g2 = (J_surf[:, 0], J_surf[:, 1])
            detJ_surface = np.linalg.norm(np.cross(g1, g2))
            dV_m = detJ_surface * w * h_avg
            K += Bdrill.T * C_drill * Bdrill * dV_m
        for (r, s), w in zip(pts3, w3):
            Bs = self._compute_Bs(r, s)
            J_surf = self._jacobian_surface(r, s)
            g1, g2 = (J_surf[:, 0], J_surf[:, 1])
            detJ_surface = np.linalg.norm(np.cross(g1, g2))
            dV_m = detJ_surface * w * h_avg
            K += Bs.T @ C_s @ Bs * dV_m
        return K

    def m_local(self, consistent=True) -> np.ndarray:
        if self.rho is None:
            return None
        Me = np.zeros((36, 36))
        pts6 = [(0.445948490144588, 0.445948490144588), (0.10810301816807, 0.445948490144588), (0.445948490144588, 0.10810301816807), (0.091576213509771, 0.091576213509771), (0.816847572980459, 0.091576213509771), (0.091576213509771, 0.816847572980459)]
        w6 = [0.11169079483905, 0.11169079483905, 0.11169079483905, 0.054975871827661, 0.054975871827661, 0.054975871827661]
        h_avg = np.mean(self.h)
        rho_h = self.rho * h_avg
        rho_I = self.rho * h_avg ** 3 / 12.0
        for (r, s), w in zip(pts6, w6):
            N, _, _ = self._shape_functions(r, s)
            J_surf = self._jacobian_surface(r, s)
            g1, g2 = (J_surf[:, 0], J_surf[:, 1])
            dA = np.linalg.norm(np.cross(g1, g2)) * w
            for i in range(6):
                for j in range(6):
                    Me[6 * i + 0, 6 * j + 0] += rho_h * N[i] * N[j] * dA
                    Me[6 * i + 1, 6 * j + 1] += rho_h * N[i] * N[j] * dA
                    Me[6 * i + 2, 6 * j + 2] += rho_h * N[i] * N[j] * dA
                    Me[6 * i + 3, 6 * j + 3] += rho_I * N[i] * N[j] * dA
                    Me[6 * i + 4, 6 * j + 4] += rho_I * N[i] * N[j] * dA
                    Me[6 * i + 5, 6 * j + 5] += rho_I * N[i] * N[j] * dA
        return Me

    def f_pressure(self, p=1.0) -> np.ndarray:
        fe = np.zeros(36)
        pts6 = [(0.445948490144588, 0.445948490144588), (0.10810301816807, 0.445948490144588), (0.445948490144588, 0.10810301816807), (0.091576213509771, 0.091576213509771), (0.816847572980459, 0.091576213509771), (0.091576213509771, 0.816847572980459)]
        w6 = [0.11169079483905, 0.11169079483905, 0.11169079483905, 0.054975871827661, 0.054975871827661, 0.054975871827661]
        for (r, s), w in zip(pts6, w6):
            N, _, _ = self._shape_functions(r, s)
            J_surf = self._jacobian_surface(r, s)
            g1, g2 = (J_surf[:, 0], J_surf[:, 1])
            normal = np.cross(g1, g2)
            detJ = np.linalg.norm(normal)
            n_unit = normal / detJ if detJ > 1e-14 else self.E3
            dA = detJ * w
            for k in range(6):
                fe[6 * k:6 * k + 3] += p * N[k] * n_unit * dA
        return fe

    def T_matrix(self) -> np.ndarray:
        return np.eye(36)

    def f_local(self) -> np.ndarray:
        return np.zeros(36)

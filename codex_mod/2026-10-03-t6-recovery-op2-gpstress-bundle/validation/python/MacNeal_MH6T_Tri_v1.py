"""
MacNeal_MH6T_Tri_v1.py
=======================
A 6-node triangular shell element using MacNeal's LINE-INTEGRATION
assumed-strain method for membrane and transverse shear -- a genuinely
different technique from MITC6_Tri_v1's point-tying, giving a third,
methodologically distinct competitor in the battle (alongside the plain
displacement-based Simo1993_Tri6 and the point-tied MITC6_Tri_v1).

Formulation source
-------------------
Li, Xiang, Izzuddin, Vu-Quoc, Zhuo & Zhang (2015), "A 6-node
co-rotational triangular elasto-plastic shell element", Comput. Mech.,
Section 3.2 and Appendix 1 -- specifically the LINEAR-ELASTIC assumed
membrane strain (Eqs. 26-38, Table 1, Fig. 2a) and assumed transverse
shear strain (Eqs. 39-43, Fig. 2b), which trace back to MacNeal (1982),
"Derivation of element stiffness matrices by assumed strain
distribution".

Only the linear-elastic membrane + shear treatment is taken from that
paper -- NOT its co-rotational large-displacement / elasto-plastic
machinery, which is out of scope for this linear-static battle harness.
Bending and drilling use the SAME Simo/MITC-style covariant director
kinematics (cross-product-based curvature/drilling formulas) as
Simo1993_Tri6_ShellElement_v1p8 / MITC6_Tri_v1 -- they are NOT derived
from Li et al. This file's OWN copies of these two routines have since
been corrected for two bugs found via direct numerical verification
against the boundary-only patch test on the SAP2000 2-001 mesh
(disp/mom/stress all 0.0000% after fixing, vs 20-206% before):
  - _compute_Bb: cross(t0_I,g1/g2) gave curvature with the opposite
    sign from the standard Mindlin convention (thx=dw/dy, thy=-dw/dx);
    fixed by swapping the cross-product argument order.
  - _compute_Bdrill: dN_dx1/dN_dx2 from _compute_derivatives are
    expansion COEFFICIENTS in the (non-orthonormal) {g1,g2} basis, not
    physical Cartesian dN/dx,dN/dy -- they were used as if they already
    were; fixed by projecting through g1,g2 onto (self.E1,self.E2)
    before use.
Whether the CURRENT Simo1993_Tri6_ShellElement_v1p8 / MITC6_Tri_v1
files carry the same two bugs has NOT been checked here (out of scope
for this file) -- so "identical to" no longer holds without also
verifying/patching those siblings the same way. Also not yet checked
for this file specifically: whether _compute_director_vectors /
_local_basis_at_point need the same "+Z" orientation tie-break applied
to Simo1993_Q8_ShellElement_v1p8_standalone.py's e3/V_n construction
(bug found there for meshes whose node winding gives cross(g1,g2)
pointing -Z). The SAP2000 2-001 T6 sub-mesh happens to give E3=+Z
naturally for every sub-triangle here, so this hasn't been exercised --
it remains an open, unverified item for meshes with different winding.

Flat-facet simplification
--------------------------
Li et al.'s method is explicitly derived assuming the element's four
sub-triangular facets are flat (stated directly under Eq. 34). This
implementation uses a single ELEMENT-CONSTANT local Cartesian frame
(computed once, at the centroid) for the membrane line-integration and
the sub-triangle normal g_l in the shear line-integration, rather than
tracking each sub-triangle's own facet plane -- a reasonable
simplification consistent with the paper's own flat-facet premise, and
exact for genuinely flat elements (all patch/plate tests below).
"""

import numpy as np
from core import Element

# ---- Table 1 of Li et al. (2015): natural coords of points A-I ----
_PT_XI = {'A': 0.25, 'B': 0.75, 'C': 0.75, 'D': 0.25, 'E': 0.0,
          'F': 0.0, 'G': 0.25, 'H': 0.5, 'I': 0.25}
_PT_ETA = {'A': 0.0, 'B': 0.0, 'C': 0.25, 'D': 0.75, 'E': 0.75,
           'F': 0.25, 'G': 0.25, 'H': 0.25, 'I': 0.5}
# Eq. 34/37: (k, i, j) ordered triplets, 1-indexed node numbers -> 0-indexed
_MEMBRANE_EDGES = [
    ('A', 0, 3), ('B', 3, 1), ('C', 1, 4), ('D', 4, 2), ('E', 2, 5),
    ('F', 5, 0), ('G', 5, 3), ('H', 3, 4), ('I', 4, 5),
]
# Eq. 39/42: (k, i, j) for the six shear line-integration edges (A-F only)
_SHEAR_EDGES = [
    ('A', 0, 3), ('B', 3, 1), ('C', 1, 4), ('D', 4, 2), ('E', 2, 5), ('F', 5, 0),
]


class MacNeal_MH6T_Tri_v1(Element):

    def __init__(self, eid, nodes, E, nu, h, rho=None, beta_drill=0.02,
                 nodal_normals=None):
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
        self.nodal_normals = nodal_normals
        self._compute_director_vectors()
        self._compute_local_bases()
        self.E1, self.E2, self.E3 = self._local_basis_at_point(1.0 / 3.0, 1.0 / 3.0)
        self._setup_macneal_tying()

    # ---- machinery shared with Simo1993_Tri6_v1.8 / MITC6_Tri_v1 ----
    def _get_node_normal(self, k):
        if self.nodal_normals is not None:
            node_id = getattr(self.nodes[k], 'nid', None)
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
            # FIX bug #5 (KELALAIAN -- baru ketahuan lewat mesh v4/CW,
            # tidak pernah ketahuan di mesh asli kita yang selalu CCW
            # natural +Z): tie-break +Z, sama seperti file T6/Q8 lain.
            if normal[2] < 0: normal = -normal
            self.V_n[k] = normal / norm

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
        a11 = np.dot(g1, g1); a22 = np.dot(g2, g2); a12 = np.dot(g1, g2)
        A = np.array([[a11, a12], [a12, a22]])
        _, dNr, dNs = self._shape_functions(r, s)
        dN_dr = np.vstack((dNr, dNs))
        dN_dx_local = np.linalg.inv(A) @ dN_dr
        return dN_dx_local[0, :], dN_dx_local[1, :]

    def _local_basis_at_point(self, r, s):
        J = self._jacobian_surface(r, s)
        g1, g2 = J[:, 0], J[:, 1]
        e3 = np.cross(g1, g2); e3 /= np.linalg.norm(e3)
        # FIX bug #5 (kelalaian sama seperti di atas): tie-break +Z DI
        # SINI juga supaya konsisten dgn V_n (self.E1/E2/E3 dihitung
        # sekali di __init__ lewat method ini).
        if e3[2] < 0: e3 = -e3
        e1 = g1 / np.linalg.norm(g1)
        e2 = np.cross(e3, e1)
        return e1, e2, e3

    def _covariant_maps(self, r, s):
        J_surf = self._jacobian_surface(r, s)
        g1, g2 = J_surf[:, 0], J_surf[:, 1]
        e1, e2, _ = self._local_basis_at_point(r, s)
        a11 = g1 @ g1; a22 = g2 @ g2; a12 = g1 @ g2
        det = a11 * a22 - a12 * a12
        if abs(det) < 1e-14:
            det = 1e-14
        a11_inv = a22 / det; a22_inv = a11 / det; a12_inv = -a12 / det
        c1_1 = a11_inv * (g1 @ e1) + a12_inv * (g2 @ e1)
        c1_2 = a11_inv * (g1 @ e2) + a12_inv * (g2 @ e2)
        c2_1 = a12_inv * (g1 @ e1) + a22_inv * (g2 @ e1)
        c2_2 = a12_inv * (g1 @ e2) + a22_inv * (g2 @ e2)
        return g1, g2, (c1_1, c1_2, c2_1, c2_2)

    def _tensor_to_physical(self, comp11, comp22, comp12, c):
        c1_1, c1_2, c2_1, c2_2 = c
        xx = (c1_1**2)*comp11 + (c2_1**2)*comp22 + 2*c1_1*c2_1*comp12
        yy = (c1_2**2)*comp11 + (c2_2**2)*comp22 + 2*c1_2*c2_2*comp12
        xy = 2.0*(c1_1*c1_2*comp11 + c2_1*c2_2*comp22 + (c1_1*c2_2+c2_1*c1_2)*comp12)
        return xx, yy, xy

    # -----------------------------------------------------------
    #  MacNeal line-integration setup (Li et al. 2015, Eqs. 26-43)
    #  Precomputed once at construction -- purely geometric, element-
    #  constant flat-facet frame (see module docstring).
    # -----------------------------------------------------------
    def _setup_macneal_tying(self):
        coords = np.array([[n.x, n.y, n.z] for n in self.nodes])
        origin = coords[0]
        xy_local = np.array([
            [(c - origin) @ self.E1, (c - origin) @ self.E2] for c in coords
        ])  # (6,2) flat local coordinates

        # ---- Membrane: Eq. 34-38 ----
        n_pts = len(_MEMBRANE_EDGES)
        Bk = np.zeros((n_pts, 36))
        Gamma = np.zeros((n_pts, n_pts))
        for row, (k, i, j) in enumerate(_MEMBRANE_EDGES):
            a = xy_local[j] - xy_local[i]
            L2 = a[0]**2 + a[1]**2
            Bk[row, 6*i:6*i+3] = -(a[0]*self.E1 + a[1]*self.E2) / L2
            Bk[row, 6*j:6*j+3] = (a[0]*self.E1 + a[1]*self.E2) / L2
            c_k, s_k = a[0]/np.sqrt(L2), a[1]/np.sqrt(L2)
            xi_k, eta_k = _PT_XI[k], _PT_ETA[k]
            Gamma[row, :] = [
                c_k**2, xi_k*c_k**2, eta_k*c_k**2,
                s_k**2, xi_k*s_k**2, eta_k*s_k**2,
                c_k*s_k, xi_k*c_k*s_k, eta_k*c_k*s_k,
            ]
        alpha_of_U = np.linalg.solve(Gamma, Bk)  # (9,36): dalpha/dU
        self._mem_alpha_of_U = alpha_of_U

        # ---- Shear: Eq. 39-43, with PER-SUB-TRIANGLE facet normals ----
        # Li et al. Eq. 40a-d: g1,g2,g3 are the normals of the three
        # "corner" sub-triangles (1,4,6), (4,2,5), (5,3,6) -- NOT one
        # element-constant frame. Using a single flat E3 for every edge
        # (the earlier simplification) is exact only when the element
        # itself is flat, and was confirmed to break down on genuinely
        # warped elements (Problem 2-004 twisted beam: Fy error went
        # from +3.0%/+0.5% for SimoT6/MITC6_v1 to -41%/-24% here once
        # the underlying mesh geometry was corrected to be properly
        # twisted). Eq. 39's quadruplet (k,i,j,l) maps each shear edge
        # to the correct sub-triangle normal g_l.
        def _tri_normal(p, q, r_pt):
            v1 = coords[q] - coords[p]
            v2 = coords[r_pt] - coords[p]
            n = np.cross(v1, v2)
            nrm = np.linalg.norm(n)
            if nrm <= 1e-14:
                return self.E3
            n = n / nrm
            # FIX (kelalaian sama seperti V_n/_local_basis_at_point):
            # samakan sisi dengan self.E3 (sudah di-tie-break +Z), bukan
            # cross-product mentah tanpa acuan. Tanpa ini, g1/g2/g3 di
            # sini bisa berlawanan arah dgn self.E3 untuk mesh dgn
            # winding CW (natural -Z) -- ketahuan lewat uji Kirchhoff-
            # consistency (gamma harus ~0 utk medan lentur murni) di
            # mesh v4/DUEL3.dat, TIDAK ketahuan di mesh asli yg selalu
            # CCW natural +Z (di situ self.E3 dan cross mentah kebetulan
            # sudah searah, jadi bug ini tersembunyi).
            if n @ self.E3 < 0: n = -n
            return n

        g1 = _tri_normal(0, 3, 5)   # sub-triangle (1,4,6)
        g2 = _tri_normal(3, 1, 4)   # sub-triangle (4,2,5)
        g3 = _tri_normal(4, 2, 5)   # sub-triangle (5,3,6)
        _EDGE_G = {'A': g1, 'B': g2, 'C': g2, 'D': g3, 'E': g3, 'F': g1}

        n_s = len(_SHEAR_EDGES)
        Gk = np.zeros((n_s, 36))
        Omega = np.zeros((n_s, n_s))
        for row, (k, i, j) in enumerate(_SHEAR_EDGES):
            a3 = coords[j] - coords[i]              # a_ij0, full 3D
            L = np.linalg.norm(a3)
            n_i0 = self._get_node_normal(i)
            n_j0 = self._get_node_normal(j)
            g_l = _EDGE_G[k]
            # translational part: (t_j - t_i) . g_l / L
            Gk[row, 6*i:6*i+3] += -g_l / L
            Gk[row, 6*j:6*j+3] += g_l / L
            # rotational part: (r_j-r_j0 + r_i-r_i0).a_ij0 / (2L)
            #   r_i - r_i0 ~= theta_i x n_i0  =>  (theta_i x n_i0).a3
            #   = theta_i . (n_i0 x a3)   [scalar triple product]
            Gk[row, 6*i+3:6*i+6] += 0.5 * np.cross(n_i0, a3) / L
            Gk[row, 6*j+3:6*j+6] += 0.5 * np.cross(n_j0, a3) / L

            c_k, s_k = a3 @ self.E1 / L, a3 @ self.E2 / L
            xi_k, eta_k = _PT_XI[k], _PT_ETA[k]
            Omega[row, :] = [c_k, xi_k*c_k, eta_k*c_k, s_k, xi_k*s_k, eta_k*s_k]
        beta_of_U = np.linalg.solve(Omega, Gk)  # (6,36): dbeta/dU
        self._shear_beta_of_U = beta_of_U

    def _compute_Bm(self, r, s):
        P = np.array([
            [1.0, r, s, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 1.0, r, s, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 1.0, r, s],
        ])
        return P @ self._mem_alpha_of_U  # (3,36), physical xx/yy/xy in E1,E2 frame

    def _compute_Bs(self, r, s):
        Q = np.array([
            [1.0, r, s, 0, 0, 0],
            [0, 0, 0, 1.0, r, s],
        ])
        return Q @ self._shear_beta_of_U  # (2,36), physical xz/yz in E1,E2 frame

    # -----------------------------------------------------------
    #  Bending + drilling: unchanged (same as Simo1993_Tri6_v1.8 /
    #  MITC6_Tri_v1) -- only membrane/shear are MacNeal-treated here.
    # -----------------------------------------------------------
    def _compute_Bdrill(self, r, s):
        # FIX (bug #7, pola SAMA dgn Q8 KikuchiMacNeal_Q8_ShellElement_v1):
        # dN_dx1,dN_dx2 dari _compute_derivatives adalah KOEFISIEN EKSPANSI
        # terhadap basis {g1,g2} (bukan basis ortonormal), BUKAN turunan
        # Cartesian dN/dx,dN/dy langsung. Dipakai sebelumnya seolah-olah
        # sudah Cartesian -- salah. Perlu diproyeksikan dulu ke e1,e2:
        # dN/dx_fisik = dNdx1*(g1.e1)+dNdx2*(g2.e1), dN/dy_fisik sama pola.
        # Diverifikasi numerik: dgn Bdrill dimatikan, max stress err patch
        # test membran turun dari 20-35% ke 7e-13% (noise) -- konfirmasi
        # Bdrill sumber TUNGGAL bug ini. e1,e2 dipakai FROZEN (self.E1/E2)
        # juga, konsisten dgn fix Q8.
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
            # FIX Bug #3: cross(t0_I, g1/g2) memberi kelengkungan dengan
            # tanda TERBALIK relatif konvensi Mindlin (thx=dw/dy, thy=-dw/dx;
            # kxx=-d(thy)/dx, kyy=+d(thx)/dy) yang dipakai field_bending()
            # test dan MacNealQ8_1992_native (yang PASS 0.0000% di test
            # yang sama). Dibuktikan lewat debug_t6_per_element.py: SEMUA
            # 12 elemen T6 (winding CCW benar, E3=+Z benar, kontribusi
            # translasional t0_xi1/t0_xi2 = 0 identik untuk mesh datar)
            # menunjukkan mxx/myy/mxy = -1.0000x referensi secara SERAGAM
            # -- sidik jari konvensi cross-product yang berlawanan arah,
            # bukan bug per-elemen. cross(a,b)=-cross(b,a), jadi tukar
            # urutan argumen adalah fix minimal yang tepat sasaran.
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

    # -----------------------------------------------------------
    #  k_local, m_local, f_pressure, T_matrix, f_local -- identical
    #  integration scheme to Simo1993_Tri6_v1.8 / MITC6_Tri_v1
    # -----------------------------------------------------------
    def k_local(self) -> np.ndarray:
        K = np.zeros((36, 36))
        pts3 = [(1/6, 1/6), (2/3, 1/6), (1/6, 2/3)]
        w3 = [1/6, 1/6, 1/6]
        pts6 = [
            (0.445948490144588, 0.445948490144588),
            (0.10810301816807, 0.445948490144588),
            (0.445948490144588, 0.10810301816807),
            (0.091576213509771, 0.091576213509771),
            (0.816847572980459, 0.091576213509771),
            (0.091576213509771, 0.816847572980459),
        ]
        w6 = [0.11169079483905, 0.11169079483905, 0.11169079483905,
              0.054975871827661, 0.054975871827661, 0.054975871827661]
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
        w6 = [0.11169079483905, 0.11169079483905, 0.11169079483905,
              0.054975871827661, 0.054975871827661, 0.054975871827661]
        h_avg = np.mean(self.h)
        rho_h = self.rho * h_avg
        rho_I = self.rho * h_avg ** 3 / 12.0
        for (r, s), w in zip(pts6, w6):
            N, _, _ = self._shape_functions(r, s)
            J_surf = self._jacobian_surface(r, s)
            g1, g2 = J_surf[:, 0], J_surf[:, 1]
            dA = np.linalg.norm(np.cross(g1, g2)) * w
            for i in range(6):
                for j in range(6):
                    Me[6*i+0, 6*j+0] += rho_h * N[i]*N[j]*dA
                    Me[6*i+1, 6*j+1] += rho_h * N[i]*N[j]*dA
                    Me[6*i+2, 6*j+2] += rho_h * N[i]*N[j]*dA
                    Me[6*i+3, 6*j+3] += rho_I * N[i]*N[j]*dA
                    Me[6*i+4, 6*j+4] += rho_I * N[i]*N[j]*dA
                    Me[6*i+5, 6*j+5] += rho_I * N[i]*N[j]*dA
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
        w6 = [0.11169079483905, 0.11169079483905, 0.11169079483905,
              0.054975871827661, 0.054975871827661, 0.054975871827661]
        for (r, s), w in zip(pts6, w6):
            N, _, _ = self._shape_functions(r, s)
            J_surf = self._jacobian_surface(r, s)
            g1, g2 = J_surf[:, 0], J_surf[:, 1]
            normal = np.cross(g1, g2)
            detJ = np.linalg.norm(normal)
            n_unit = normal / detJ if detJ > 1e-14 else self.E3
            dA = detJ * w
            for k in range(6):
                fe[6*k:6*k+3] += p * N[k] * n_unit * dA
        return fe

    def T_matrix(self) -> np.ndarray:
        return np.eye(36)

    def f_local(self) -> np.ndarray:
        return np.zeros(36)

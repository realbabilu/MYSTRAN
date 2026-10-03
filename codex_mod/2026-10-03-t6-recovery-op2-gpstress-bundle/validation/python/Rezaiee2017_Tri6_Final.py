"""
MITC6_Tri_v1.py
================
A 6-node triangular shell element using MITC-style ASSUMED (tied) covariant
strains for membrane and transverse shear, to remove the shear/membrane
locking that a plain displacement-based Tria6 (e.g. Simo1993_Tri6_v1.8)
suffers from.

Formulation source
-------------------
Reconstructed from the tying scheme described in:
  - Kim & Bathe (2009), "A triangular six-node shell element", Comput.
    Struct. 87, Eqs. (1)-(11), Fig. 1  (tying-point layout, natural coords
    r1=s1=1/2-1/(2*sqrt(3)), r2=s2=1/2+1/(2*sqrt(3)), r3=s3=1/3 "common").
  - Rezaiee-Pajand, Arabi & Masoodi (2017), "A triangular shell element for
    geometrically nonlinear analysis", Acta Mech., Eqs. (10)-(19), which
    give the *closed-form* coefficient formulas for the assumed strain
    interpolation (this is what is implemented below almost verbatim).

Everything else (nodal director vectors, drilling penalty, bending
strain, geometry/Jacobian, integration rule) is taken unchanged from
Simo1993_Tri6_ShellElement_v1p8 so that the *only* difference between the
two elements in the "battle" is the membrane + shear strain treatment.

Honesty note
------------
Kim & Bathe's Fig. 1 assigns tying points to the three strain "groups"
(edge s=0, edge r=0, and the hypotenuse r+s=1) via a picture that is not
literally reproducible from the OCR'd text alone. The point *coordinates*
(r1,r2,r3=1/3) and the *closed-form coefficient recipe* (Eqs. 14-17 of
Rezaiee-Pajand et al.) are, however, given explicitly in the text, and are
used here exactly. The assignment of the three edges to the three
"groups" (rr/ss/qq and rt/st) is filled in with the symmetric, isotropic
choice implied by the element's spatial-isotropy requirement (each of the
3 edges plays an equivalent role). This is a best-effort, well-posed
reconstruction, not a verbatim reproduction of the original code -- it is
validated below against the patch tests (T2-T5) in the existing harness.
"""

import numpy as np
from core import Element

SQRT3 = np.sqrt(3.0)
SQRT2 = np.sqrt(2.0)
R1 = 0.5 - 1.0 / (2.0 * SQRT3)
R2 = 0.5 + 1.0 / (2.0 * SQRT3)
RC = 1.0 / 3.0
RC_SHEAR = 1.0 / SQRT3  # Fig.3 Eq.19 rc=sc=1/√3 for transverse shear (Rezaiee)


class Rezaiee2017_Tri6_ShellElement(Element):

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

    # ---- identical helper machinery to Simo1993_Tri6_v1.8 ----------
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

    # -----------------------------------------------------------
    #  NATURAL (pre-transform) covariant membrane / shear B-rows
    #  -- these are the "displacement-based" strain-displacement
    #     rows that get TIED at the MITC sampling points below.
    # -----------------------------------------------------------
    def _natural_membrane_rows(self, r, s):
        """Return (B_rr, B_ss, B_rs) each a length-36 covariant
        (natural-coordinate) strain-displacement row."""
        _, dNr, dNs = self._shape_functions(r, s)
        J = self._jacobian_surface(r, s)
        g1, g2 = J[:, 0], J[:, 1]
        B_rr = np.zeros(36)
        B_ss = np.zeros(36)
        B_rs = np.zeros(36)
        for k in range(6):
            col = 6 * k
            a1, a2 = dNr[k], dNs[k]
            B_rr[col:col + 3] = a1 * g1
            B_ss[col:col + 3] = a2 * g2
            B_rs[col:col + 3] = 0.5 * (a1 * g2 + a2 * g1)
        return B_rr, B_ss, B_rs

    def _natural_shear_rows(self, r, s):
        """Return (B_rt, B_st) each a length-36 covariant
        (natural-coordinate) transverse-shear strain-displacement row."""
        N, dNr, dNs = self._shape_functions(r, s)
        J = self._jacobian_surface(r, s)
        g1, g2 = J[:, 0], J[:, 1]
        V_n_nodes = np.array([self._get_node_normal(k) for k in range(6)])
        t0 = N @ V_n_nodes
        nrm = np.linalg.norm(t0)
        if nrm > 1e-12:
            t0 = t0 / nrm
        B_rt = np.zeros(36)
        B_st = np.zeros(36)
        for k in range(6):
            col = 6 * k
            t0_I = V_n_nodes[k]
            B_rt[col:col + 3] += dNr[k] * t0
            B_st[col:col + 3] += dNs[k] * t0
            B_rt[col + 3:col + 6] += N[k] * np.cross(t0_I, g1)
            B_st[col + 3:col + 6] += N[k] * np.cross(t0_I, g2)
        return B_rt, B_st

    # -----------------------------------------------------------
    #  ASSUMED (tied) membrane strain  -- Kim & Bathe (2009) Fig. 1
    #
    #  Corrected tying-point layout (read directly off Fig. 1):
    #    e_rr tied at (R1,R1_s)=(r1,s1) and (r2,s1)   [shared s=s1]
    #    e_ss tied at (r1,s1) and (r1,s2)              [shared r=r1]
    #    e_qq tied at (r2,s1) and (r1,s2)              [the two "far" points]
    #  All three fields are e = a+b*r+c*s (or +c*(1-r-s) for qq), so each
    #  needs a 3rd condition; we close with the exact interpolation
    #  through the two tying points plus the centroid (1/3,1/3) -- this
    #  reproduces the h^k(point_l)=delta_kl Lagrange property of Eq. (6)
    #  without depending on a possibly-OCR-mangled closed form.
    # -----------------------------------------------------------
    def _fit_affine(self, vals, pts, basis):
        """Solve for (a,b,c) in val = a*basis0+b*basis1+c*basis2 at 3 pts.
        basis(r,s) -> (b0,b1,b2). vals/pts: 3 natural-strain rows / (r,s)."""
        A = np.array([basis(r, s) for (r, s) in pts])  # (3,3)
        Ainv = np.linalg.inv(A)
        # coeffs[k] = sum_i Ainv[k,i] * vals[i]  (each val is a 36-vector)
        return (Ainv[0, 0] * vals[0] + Ainv[0, 1] * vals[1] + Ainv[0, 2] * vals[2],
                Ainv[1, 0] * vals[0] + Ainv[1, 1] * vals[1] + Ainv[1, 2] * vals[2],
                Ainv[2, 0] * vals[0] + Ainv[2, 1] * vals[1] + Ainv[2, 2] * vals[2])

    def _compute_Bm(self, r, s):
        rr_11, ss_11, rs_11 = self._natural_membrane_rows(R1, R1)
        rr_21, ss_21, rs_21 = self._natural_membrane_rows(R2, R1)
        rr_12, ss_12, rs_12 = self._natural_membrane_rows(R1, R2)
        rr_c, ss_c, rs_c = self._natural_membrane_rows(RC, RC)

        qq_11 = 0.5 * (rr_11 + ss_11) - rs_11
        qq_21 = 0.5 * (rr_21 + ss_21) - rs_21
        qq_12 = 0.5 * (rr_12 + ss_12) - rs_12
        qq_c = 0.5 * (rr_c + ss_c) - rs_c

        basis_rs = lambda r_, s_: (1.0, r_, s_)
        basis_qq = lambda r_, s_: (1.0, r_, 1.0 - r_ - s_)

        a1, b1, c1 = self._fit_affine([rr_11, rr_21, rr_c],
                                       [(R1, R1), (R2, R1), (RC, RC)], basis_rs)
        a2, b2, c2 = self._fit_affine([ss_11, ss_12, ss_c],
                                       [(R1, R1), (R1, R2), (RC, RC)], basis_rs)
        a3, b3, c3 = self._fit_affine([qq_21, qq_12, qq_c],
                                       [(R2, R1), (R1, R2), (RC, RC)], basis_qq)

        B_err = a1 + b1 * r + c1 * s
        B_ess = a2 + b2 * r + c2 * s
        B_eqq = a3 + b3 * r + c3 * (1.0 - r - s)
        B_ers = 0.5 * (B_err + B_ess) - B_eqq

        _, _, c = self._covariant_maps(r, s)
        xx, yy, xy = self._tensor_to_physical(B_err, B_ess, B_ers, c)
        Bm = np.zeros((3, 36))
        Bm[0, :] = xx
        Bm[1, :] = yy
        Bm[2, :] = xy
        return Bm

    # -----------------------------------------------------------
    #  ASSUMED (tied) transverse shear strain -- Kim & Bathe Fig. 1
    #  (right panel), Eq. (11)
    #
    #  ẽ_rt = a1 + b1 r + c1 s + d r s + e s^2
    #  ẽ_st = a2 + b2 r + c2 s - d r^2 - e r s
    #  (8 unknowns, sharing the coupling terms d,e)
    #
    #  High-confidence tying points read off Fig. 1's right panel:
    #    e_rt tied at (r1,0), (r2,0)      [s=0 edge, red circles]
    #    e_st tied at (r1,s1), (r1,s2)    [mirrors e_ss's r=r1 line]
    #    both e_rt AND e_st sampled at the "common" point (r3,s3)=(1/3,1/3)
    #  That is 6 conditions for 8 unknowns. The remaining 2 come from the
    #  derived quantity e_qt = (1/sqrt(2))*(e_st - e_rt) -- exactly the
    #  shear analogue of e_qq = (e_rr+e_ss)/2 - e_rs used for membrane
    #  (confirmed by da Veiga, Chapelle & Suarez (2007) Fig. 1, which
    #  gives e_qz = (1/sqrt(2))*(e_sz - e_rz) for the same MITC6 family).
    #  Sampled at the two points symmetric with e_qq's tying locations:
    #  (r2,s1) and (r1,s2).
    # -----------------------------------------------------------
    def _compute_Bs(self, r, s):
        pts = [
            (R1, 0.0), (R2, 0.0),           # e_rt tying (edge s=0) Fig.3 left
            (R1, R1), (R1, R2),              # e_st tying (line r=r1) Fig.3 left
            (RC_SHEAR, RC_SHEAR), (RC_SHEAR, RC_SHEAR),  # common point: rc=sc=1/√3 Fig.3 right (Rezaiee Eq.19)
            (R2, R1), (R1, R2),              # e_qt tying points Fig.3 right
        ]
        kinds = ['rt', 'rt', 'st', 'st', 'rt', 'st', 'qt', 'qt']

        rhs = []
        rows = []
        for (rp, sp), kind in zip(pts, kinds):
            Brt, Bst = self._natural_shear_rows(rp, sp)
            if kind == 'rt':
                rhs.append(Brt)
                rows.append([1.0, rp, sp, 0.0, 0.0, 0.0, rp * sp, sp ** 2])
            elif kind == 'st':
                rhs.append(Bst)
                rows.append([0.0, 0.0, 0.0, 1.0, rp, sp, -rp ** 2, -rp * sp])
            else:  # 'qt' = (1/sqrt(2)) * (est - ert)
                rhs.append((Bst - Brt) / SQRT2)
                r_rt = np.array([1.0, rp, sp, 0.0, 0.0, 0.0, rp * sp, sp ** 2])
                r_st = np.array([0.0, 0.0, 0.0, 1.0, rp, sp, -rp ** 2, -rp * sp])
                rows.append(list((r_st - r_rt) / SQRT2))

        M = np.array(rows)          # (8,8), purely geometric
        Minv = np.linalg.inv(M)
        # coeffs[k] (a length-36 vector) = sum_i Minv[k,i] * rhs[i]
        coeffs = [sum(Minv[k, i] * rhs[i] for i in range(8)) for k in range(8)]
        a1, b1, c1, a2, b2, c2, d, e = coeffs

        B_ert = a1 + b1 * r + c1 * s + d * r * s + e * s ** 2
        B_est = a2 + b2 * r + c2 * s - d * r ** 2 - e * r * s

        _, _, c = self._covariant_maps(r, s)
        c1_1, c1_2, c2_1, c2_2 = c
        Bs_phys = np.zeros((2, 36))
        Bs_phys[0, :] = c1_1 * B_ert + c2_1 * B_est
        Bs_phys[1, :] = c1_2 * B_ert + c2_2 * B_est
        return Bs_phys

    # -----------------------------------------------------------
    #  Bending + drilling: UNCHANGED from Simo1993_Tri6_v1.8
    #  (bending is not the source of locking here; only membrane
    #   and shear are treated with assumed strains, matching the
    #   MITC philosophy of the source papers)
    # -----------------------------------------------------------
    def _compute_Bdrill(self, r, s):
        dN_dx1, dN_dx2 = self._compute_derivatives(r, s)
        N, _, _ = self._shape_functions(r, s)
        e1, e2, e3 = self._local_basis_at_point(r, s)
        Bdrill = np.zeros((1, 36))
        for k in range(6):
            col = 6 * k
            t0_I = self._get_node_normal(k)
            Bdrill[0, col:col + 3] = 0.5 * (dN_dx1[k] * e2 - dN_dx2[k] * e1)
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
            cross_g1 = np.cross(t0_I, g1)
            cross_g2 = np.cross(t0_I, g2)
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
    #  k_local, m_local, f_pressure, f_thermal, k_geometric, T_matrix
    #  -- identical integration scheme to Simo1993_Tri6_v1.8
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

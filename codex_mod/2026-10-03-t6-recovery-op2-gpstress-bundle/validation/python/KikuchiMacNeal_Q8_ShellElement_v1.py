"""
KikuchiMacNeal_Q8_ShellElement_v1.py
=====================================
CQUAD8 contender: Simo1993_Q8_ShellElement_v1p8's covariant/ANS shell
formulation (kept verbatim -- transverse shear tying, drilling, mass,
pressure, thermal, geometric stiffness), with the FIELD interpolation
(the shape functions that weight nodal DOFs when building Bm, Bb, Bs,
Bdrill) swapped for the modified 8-node shape functions of either:

    variant='kikuchi'         Kikuchi, Okabe & Fujio (1999), closed form
    variant='macneal_harder'  MacNeal & Harder (1992), general numeric

Architecture (see macneal_kikuchi_shapefun.py docstring): geometry /
Jacobian mapping (g1, g2, the covariant metric, the director field's
own flat-normal construction) is LEFT UNTOUCHED -- it still uses the
standard Q8 serendipity map, exactly like Simo1993_Q8. Only the
"which nodal DOF gets how much weight at this integration point" part
of Bm/Bb/Bs/Bdrill is modified. This mirrors the two papers exactly:
they fix serendipity's inability to reproduce Cartesian quadratics in
the FIELD, not the geometry.
"""
import numpy as np
from Simo1993_Q8_ShellElement_v1p8_standalone import Simo1993_Q8_ShellElement_v1p8
from macneal_kikuchi_shapefun import (kikuchi_shape_functions,
                                       macneal_harder_shape_functions)


class KikuchiMacNeal_Q8_ShellElement_v1(Simo1993_Q8_ShellElement_v1p8):

    def __init__(self, eid, nodes, E, nu, h, rho=None, beta_drill=0.02,
                 variant='kikuchi', freeze_bb=True):  # variant='kikuchi' or 'macneal_harder'
        if len(nodes) != 8:
            raise AssertionError("KikuchiMacNeal_Q8 only supports 8-node Q8 "
                                  f"(got {len(nodes)}); use the T6 sibling for triangles.")
        self.variant = variant
        # freeze_bb: TERBUKTI WAJIB (default True), bukan lagi eksperimen.
        # Sama seperti Bm, Bb butuh frame TETAP (self.E1/E2) -- BUKAN
        # _covariant_maps per-titik-Gauss yang berputar -- karena _xy8
        # (dasar koreksi Kikuchi/MH) dibangun dari frame tetap itu juga.
        # Diverifikasi numerik: dengan freeze_bb=True, patch test bending
        # 2-001 (mesh SAP2000 terdistorsi ekstrem) turun dari residual
        # 6.1412% ke PERSIS 0.0000% (disp DAN moment), match sempurna
        # dengan MacNealQ8_1992_native. freeze_bb=False disisakan hanya
        # untuk kebutuhan riset/regresi, bukan direkomendasikan.
        self.freeze_bb = freeze_bb
        super().__init__(eid, nodes, E, nu, h, rho=rho, beta_drill=beta_drill)
        # Local planar (x,y) node coordinates -- projected onto the
        # element's own tangent basis (E1,E2) at its center, exactly the
        # "x,y plane parallel to the diagonals" convention both papers use.
        coords = np.array([[n.x, n.y, n.z] for n in self.nodes])
        centroid = coords.mean(axis=0)
        rel = coords - centroid
        self._xy8 = np.column_stack([rel @ self.E1, rel @ self.E2])

    def _field_shape_functions(self, r, s):
        """Modified (N, dNr, dNs) for interpolating/differentiating FIELD
        quantities (translations, director components) -- NOT geometry."""
        if self.variant == 'kikuchi':
            return kikuchi_shape_functions(r, s, self._xy8[:4])
        elif self.variant == 'macneal_harder':
            return macneal_harder_shape_functions(r, s, self._xy8)
        else:
            raise ValueError(f"unknown variant {self.variant!r}")

    def _compute_derivatives_field(self, r, s):
        """Cartesian (local E1,E2) derivatives of the FIELD shape functions,
        using the standard-geometry metric tensor A (unchanged) but the
        MODIFIED shape function's parametric derivatives as the numerator
        -- mirrors the parent's _compute_derivatives but decouples the two."""
        J = self._jacobian_surface(r, s)          # standard geometry, untouched
        g1, g2 = J[:, 0], J[:, 1]
        a11 = g1 @ g1; a22 = g2 @ g2; a12 = g1 @ g2
        A = np.array([[a11, a12], [a12, a22]])
        _, dNr, dNs = self._field_shape_functions(r, s)
        dN_dr = np.vstack((dNr, dNs))
        dN_dx_local = np.linalg.inv(A) @ dN_dr
        return dN_dx_local[0, :], dN_dx_local[1, :]

    # ---- Bm, Bb, Bs, Bdrill: same algebra as parent, field N/dNr/dNs swapped ----

    def _compute_Bm(self, r, s):
        _, dNr, dNs = self._field_shape_functions(r, s)
        # FIX bug #2: WAJIB pakai frame TETAP (self.E1/E2), bukan
        # _covariant_maps(r,s) yang berputar per titik Gauss. self._xy8
        # (dipakai kikuchi_shape_functions/macneal_harder_shape_functions
        # untuk menurunkan koreksi Mi) dibangun dari self.E1/E2 tetap --
        # kalau proyeksi fisiknya pakai frame berbeda (berputar), garansi
        # kelengkapan Cartesian-quadratic dari kedua paper tidak berlaku.
        # Lihat _covariant_maps_frozen di parent (v1p8) untuk detail.
        g1, g2, c = self._covariant_maps_frozen(r, s)
        n = len(self.nodes)
        Bm = np.zeros((3, n * 6))
        for k in range(n):
            col = 6 * k
            a1, a2 = dNr[k], dNs[k]
            deps11 = a1 * g1; deps22 = a2 * g2; deps12 = 0.5 * (a1 * g2 + a2 * g1)
            dEps_xx, dEps_yy, dEps_xy = self._tensor_to_physical(deps11, deps22, deps12, c)
            Bm[0, col:col+3] = dEps_xx
            Bm[1, col:col+3] = dEps_yy
            Bm[2, col:col+3] = dEps_xy
        return Bm

    def _compute_Bb(self, r, s):
        _, dNr, dNs = self._field_shape_functions(r, s)
        # freeze_bb: lihat catatan di __init__ -- default masih covariant
        # per-titik (perilaku lama tidak berubah kecuali diminta eksplisit).
        g1, g2, c = (self._covariant_maps_frozen(r, s) if self.freeze_bb
                     else self._covariant_maps(r, s))
        t0_xi1 = dNr @ self.V_n
        t0_xi2 = dNs @ self.V_n
        n = len(self.nodes)
        Bb = np.zeros((3, n * 6))
        for k in range(n):
            col = 6 * k
            t0_I = self.V_n[k]
            a1, a2 = dNr[k], dNs[k]
            cross_g1 = np.cross(g1, t0_I); cross_g2 = np.cross(g2, t0_I)
            drho11_rot = a1 * cross_g1; drho22_rot = a2 * cross_g2
            drho12_rot = 0.5 * (a1 * cross_g2 + a2 * cross_g1)
            dK_xx_rot, dK_yy_rot, dK_xy_rot = self._tensor_to_physical(
                drho11_rot, drho22_rot, drho12_rot, c)
            Bb[0, col+3:col+6] = dK_xx_rot
            Bb[1, col+3:col+6] = dK_yy_rot
            Bb[2, col+3:col+6] = dK_xy_rot
            drho11_tr = a1 * t0_xi1; drho22_tr = a2 * t0_xi2
            drho12_tr = 0.5 * (a1 * t0_xi2 + a2 * t0_xi1)
            dK_xx_tr, dK_yy_tr, dK_xy_tr = self._tensor_to_physical(
                drho11_tr, drho22_tr, drho12_tr, c)
            Bb[0, col:col+3] = dK_xx_tr
            Bb[1, col:col+3] = dK_yy_tr
            Bb[2, col:col+3] = dK_xy_tr
        return Bb

    def _compute_Bs(self, r, s):
        # FIX bug #8: index c1_1/c1_2/c2_1/c2_2 tertukar di proyeksi
        # kovarian->fisik (pola SAMA dgn bug #1 di _tensor_to_physical,
        # tapi Bs punya implementasi proyeksi sendiri, tidak lewat situ).
        # Bs_nat[0,:] adalah kontribusi arah-alpha=1 (dNr), Bs_nat[1,:]
        # arah-alpha=2 (dNs). Proyeksi fisik yang benar: physical_i =
        # sum_alpha c_{alpha,i} * (kontribusi alpha) -- jadi physical_x
        # (i=1) harus pakai (c1_1,c2_1), BUKAN (c1_1,c1_2). Diverifikasi
        # numerik: sebelum fix, gamma_local utk w=a*x_local murni memberi
        # [0.00009929,0.00000728] (harus [0.0001,0]); sesudah fix, EKSAK
        # [0.00010000,0.00000000], cocok MacNealQ8_1992_native.
        # freeze_bb (frame proyeksi TETAP) juga dipakai, konsisten dgn Bm/Bb.
        N, dNr, dNs = self._field_shape_functions(r, s)
        g1, g2, c = self._covariant_maps_frozen(r, s)
        t0 = N @ self.V_n
        t0 /= np.linalg.norm(t0)
        n = len(self.nodes)
        Bs_nat = np.zeros((2, n * 6))
        for k in range(n):
            col = 6 * k
            t0_I = self.V_n[k]
            Bs_nat[0, col:col+3] += dNr[k] * t0
            Bs_nat[1, col:col+3] += dNs[k] * t0
            Bs_nat[0, col+3:col+6] += N[k] * np.cross(t0_I, g1)
            Bs_nat[1, col+3:col+6] += N[k] * np.cross(t0_I, g2)
        c1_1, c1_2, c2_1, c2_2 = c
        Bs_phys = np.zeros((2, n * 6))
        Bs_phys[0, :] = c1_1 * Bs_nat[0, :] + c2_1 * Bs_nat[1, :]
        Bs_phys[1, :] = c1_2 * Bs_nat[0, :] + c2_2 * Bs_nat[1, :]
        return Bs_phys

    def _compute_Bdrill(self, r, s):
        # FIX bug #7: dN_dx1,dN_dx2 dari _compute_derivatives_field adalah
        # KOEFISIEN EKSPANSI terhadap basis {g1,g2} (bukan basis ortonormal),
        # BUKAN turunan Cartesian dN/dx,dN/dy langsung. Sebelumnya dipakai
        # seolah-olah sudah Cartesian -- salah magnitude sampai puluhan kali
        # lipat (diverifikasi numerik: diff bisa >100 dibanding dNdx yang
        # benar). Perlu diproyeksikan dulu: dN/dx_fisik = dNdx1*(g1.e1) +
        # dNdx2*(g2.e1), dN/dy_fisik = dNdx1*(g1.e2) + dNdx2*(g2.e2).
        # Setelah proyeksi, diverifikasi cocok EKSAK dengan
        # MacNealQ8_1992_native (diff=0.000000 di semua node/titik uji).
        # e1,e2 dipakai FROZEN (self.E1/E2) juga, konsisten dgn fix Bm.
        dN_dx1, dN_dx2 = self._compute_derivatives_field(r, s)
        J = self._jacobian_surface(r, s); g1, g2 = J[:, 0], J[:, 1]
        e1, e2, e3 = self.E1, self.E2, self.E3
        dNdx = dN_dx1*(g1@e1) + dN_dx2*(g2@e1)
        dNdy = dN_dx1*(g1@e2) + dN_dx2*(g2@e2)
        N, _, _ = self._field_shape_functions(r, s)
        n = len(self.nodes)
        Bdrill = np.zeros((1, n * 6))
        for k in range(n):
            col = 6 * k
            t0_I = self.V_n[k]
            Bdrill[0, col:col+3] = 0.5 * (dNdx[k] * e2 - dNdy[k] * e1)
            Bdrill[0, col+3:col+6] += -N[k] * t0_I
        return Bdrill

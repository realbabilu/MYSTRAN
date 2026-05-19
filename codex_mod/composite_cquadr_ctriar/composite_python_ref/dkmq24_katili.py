"""Katili-style DKMQ24 quadrilateral shell element.

This file is based on dkmq24_element_v3.py and keeps the v3 mesh-level
averaged-normal/SNORM support.  It is the recommended baseline over v2.
"""

"""
dkmq24_element.py — Implementasi Elemen Shell DKMQ24
=====================================================
Berdasarkan:
  - Katili et al. (2015): "The development of DKMQ plate bending element for
    thick to thin shell analysis based on the Naghdi/Reissner/Mindlin shell theory"
    Finite Elements in Analysis and Design, 100, 12-27.
  - Maknun et al. (2016): "Application of DKMQ24 shell element for twist of
    thin-walled beams: comparison with Vlassov theory"

Karakteristik elemen:
  - 4-node quadrilateral shell element
  - 6 DOF per node: (U, V, W, θX, θY, θZ) — global
  - 24 DOF per element
  - Interpolasi translasi: bilinear (Ni)
  - Interpolasi rotasi: incomplete quadratic (Ni + Pk)
  - Bebas shear-locking via faktor ϕk (Discrete Kirchhoff-Mindlin method)
  - Mendukung non-coplanar nodes (warping + bending-membrane coupling)
  - Fictitious stiffness untuk θZ (MacNeal stabilization)
"""

import numpy as np
try:
    from core import Element, Node
except Exception:
    try:
        from core_sparse import Element, Node
    except Exception:
        # Last-resort light stubs for import-time introspection only.
        class Element:
            pass
        class Node:
            pass


# ─────────────────────────────────────────────────────────
#  FUNGSI INTERPOLASI
# ─────────────────────────────────────────────────────────

def shape_N(xi: float, eta: float) -> np.ndarray:
    """Bilinear shape functions N1..N4 (Table 1, Katili 2015)."""
    return 0.25 * np.array([
        (1 - xi) * (1 - eta),
        (1 + xi) * (1 - eta),
        (1 + xi) * (1 + eta),
        (1 - xi) * (1 + eta),
    ])

def shape_dN(xi: float, eta: float) -> np.ndarray:
    """dNi/dxi dan dNi/deta, shape (2,4)."""
    return 0.25 * np.array([
        [-(1 - eta),  (1 - eta), (1 + eta), -(1 + eta)],   # dN/dxi
        [-(1 - xi), -(1 + xi),  (1 + xi),   (1 - xi)],    # dN/deta
    ])

def shape_P(xi: float, eta: float) -> np.ndarray:
    """Incomplete quadratic functions P5..P8 (Table 1, Katili 2015)."""
    return 0.5 * np.array([
        (1 - xi**2) * (1 - eta),
        (1 + xi)    * (1 - eta**2),
        (1 - xi**2) * (1 + eta),
        (1 - xi)    * (1 - eta**2),
    ])

def shape_dP(xi: float, eta: float) -> np.ndarray:
    """dPk/dxi dan dPk/deta, shape (2,4)."""
    return 0.5 * np.array([
        [-2*xi*(1-eta),   (1-eta**2),    -2*xi*(1+eta),   -(1-eta**2)],
        [-(1-xi**2),     -2*eta*(1+xi),   (1-xi**2),      -2*eta*(1-xi)],
    ])

# ─────────────────────────────────────────────────────────
#  GAUSS QUADRATURE  (2×2 standard)
# ─────────────────────────────────────────────────────────

_GP = np.array([-1.0/np.sqrt(3), 1.0/np.sqrt(3)])
_GW = np.array([1.0, 1.0])
GAUSS_POINTS  = [(xi, eta, wx*we)
                 for xi, wx in zip(_GP, _GW)
                 for eta, we in zip(_GP, _GW)]


# ─────────────────────────────────────────────────────────
#  DKMQ24 ELEMENT
# ─────────────────────────────────────────────────────────

class DKMQ24(Element):
    """
    DKMQ24 Shell Element (Katili 2015).

    Parameters
    ----------
    eid      : element id
    nodes    : list[Node], length 4, CCW order
    E        : Young's modulus
    nu       : Poisson's ratio
    h        : shell thickness
    kappa    : shear correction factor (default 5/6)
    use_mindlin : include transverse shear (ϕk factor); False → thin-shell DKQ24
    """

    def __init__(self, eid: int, nodes: list, E: float, nu: float, h: float,
                 kappa: float = 5.0/6.0, use_mindlin: bool = True):
        self.eid         = eid
        self.nodes       = nodes          # 4 Node objects
        self.E           = float(E)
        self.nu          = float(nu)
        self.h           = float(h)
        self.kappa       = float(kappa)
        self.use_mindlin = use_mindlin

        # Koordinat nodal (4×3)
        self._xyz = np.array([[n.x, n.y, n.z] for n in nodes])

        # Pasang indeks sisi: sisi k menghubungkan node (i,j)
        # Sisi 5:(0-1), 6:(1-2), 7:(2-3), 8:(3-0)
        self._side_pairs = [(0,1),(1,2),(2,3),(3,0)]

        # Assembly-averaged nodal normals (shape 4×3).
        # None → _nodal_normals() falls back to local single-element estimate.
        # Set by compute_averaged_nodal_normals() before stiffness is needed.
        self._averaged_normals: np.ndarray | None = None

        # Precompute stiffness
        self._K24 = self._compute_stiffness()

    # ── Interface wajib ────────────────────────────────────

    @property
    def ndof_per_node(self) -> int:
        return 6

    def k_local(self) -> np.ndarray:
        return self._K24

    def T_matrix(self) -> np.ndarray:
        return np.eye(24)

    def k_global(self) -> np.ndarray:
        return self._K24

    def global_dof_indices(self) -> list:
        dofs = []
        for node in self.nodes:
            dofs.extend(node.dofs[:6])
        return dofs

    # ── Normal vektor nodal ───────────────────────────────

    def _nodal_normals(self) -> np.ndarray:
        """
        Return unit normals at the 4 corner nodes (shape 4×3).

        Priority:
          1. Assembly-averaged normals injected by
             compute_averaged_nodal_normals() — geometrically consistent
             across the mesh, critical for curved/twisted shells.
          2. Local single-element fallback: cross-product of the two
             element edges meeting at each node.  Sufficient for flat
             patches but can deviate from the true surface normal at
             shared nodes when adjacent elements are non-coplanar.
        """
        if self._averaged_normals is not None:
            return self._averaged_normals   # already unit-normalised

        # ── Fallback: local single-element estimate ────────
        xyz = self._xyz
        normals = np.zeros((4, 3))
        pairs = [
            (xyz[1] - xyz[0], xyz[3] - xyz[0]),  # node 0
            (xyz[2] - xyz[1], xyz[0] - xyz[1]),  # node 1
            (xyz[3] - xyz[2], xyz[1] - xyz[2]),  # node 2
            (xyz[0] - xyz[3], xyz[2] - xyz[3]),  # node 3
        ]
        for i, (a, b) in enumerate(pairs):
            n = np.cross(a, b)
            nm = np.linalg.norm(n)
            normals[i] = n / nm if nm > 1e-15 else np.array([0., 0., 1.])
        return normals

    # ── Geometri di titik Gauss ───────────────────────────

    def _geometry_at(self, xi: float, eta: float):
        """
        Return (a1, a2, n, J, Co, bc_mat, Fn_mat) di titik (xi,eta).

        a1, a2  : covariant basis vectors (3,)
        n       : unit normal (3,)
        J       : sqrt(det[a]) — Jacobian faktor
        Co      : koordinat transformasi (2,2), eq.(14) Katili 2015
        bc_mat  : curvature coupling matrix (2,2), eq.(33)
        """
        N  = shape_N(xi, eta)
        dN = shape_dN(xi, eta)   # (2,4)

        xyz      = self._xyz       # (4,3)
        normals  = self._nodal_normals()  # (4,3)

        # Covariant basis: a1 = dxp/dxi, a2 = dxp/deta
        a1 = dN[0] @ xyz   # (3,)
        a2 = dN[1] @ xyz   # (3,)

        # Unit normal
        axb  = np.cross(a1, a2)
        J    = np.linalg.norm(axb)
        n_vec = axb / J if J > 1e-15 else np.array([0., 0., 1.])

        # Metric tensor [a] (2×2)
        a11 = a1 @ a1
        a12 = a1 @ a2
        a21 = a12
        a22 = a2 @ a2
        a_mat = np.array([[a11, a12],[a21, a22]])
        det_a = a11*a22 - a12*a21

        # Contravariant basis  a^1, a^2
        inv_a = np.linalg.inv(a_mat)
        a1c = inv_a[0,0]*a1 + inv_a[0,1]*a2   # a^1
        a2c = inv_a[1,0]*a1 + inv_a[1,1]*a2   # a^2

        # Local frame t1, t2 — Katili 2015 eq.(13)
        # t1 = n × k / |n × k|,  t2 = n × t1
        # Frame ini memberikan Co off-diagonal untuk curved element → warping coupling benar
        k_hat = np.array([0., 0., 1.])
        t1 = np.cross(n_vec, k_hat)
        nm_t1 = np.linalg.norm(t1)
        if nm_t1 < 1e-10:
            # n nearly parallel to k: use j-hat as reference
            t1 = np.cross(n_vec, np.array([0., 1., 0.]))
            nm_t1 = np.linalg.norm(t1)
        t1 /= nm_t1
        t2 = np.cross(n_vec, t1)
        nm_t2 = np.linalg.norm(t2)
        if nm_t2 < 1e-15:
            t2 = np.array([0., 1., 0.])
        else:
            t2 /= nm_t2

        # Co matrix (2×2): Co_ij = a^i · t_j  (eq.14)
        Co = np.array([
            [a1c @ t1, a1c @ t2],
            [a2c @ t1, a2c @ t2],
        ])

        # Second fundamental form coupling: bc matrix (eq.33)
        # n,xi = dN/dxi @ normals,  n,eta = dN/deta @ normals
        n_xi  = dN[0] @ normals   # (3,)
        n_eta = dN[1] @ normals   # (3,)

        # bn = F0^{-1} Fn  (eq.16), simplified to 2×2 block
        # bc = bn_hat @ Co  (eq.33)
        bn11 = a1c @ n_xi
        bn12 = a1c @ n_eta
        bn21 = a2c @ n_xi
        bn22 = a2c @ n_eta

        bn_hat = np.array([[ bn22, -bn12],
                           [-bn21,  bn11]])
        bc_mat = bn_hat @ Co  # (2×2)

        return t1, t2, n_vec, J, Co, bc_mat, normals

    # ── Shear-correction factor ϕk per sisi ───────────────

    def _phi_side(self, i: int, j: int) -> float:
        """
        ϕk = (Db/Ds) * (12/Lk^2) = 2/(K(1-ν)) * (h²/Lk²)
        eq.(45) Katili 2015.
        """
        Lk = np.linalg.norm(self._xyz[j] - self._xyz[i])
        if not self.use_mindlin or Lk < 1e-15:
            return 0.0
        phi = 2.0 / (self.kappa * (1.0 - self.nu)) * (self.h**2 / Lk**2)
        return phi

    # ── RN matrix (spin-rotation coupling)  ───────────────

    @staticmethod
    def _RN(nx: float, ny: float, nz: float) -> np.ndarray:
        """
        [RN]_i: 3×3 rotation-coupling matrix  (eq.21 Katili 2015).
        Maps (θX,θY,θZ) → βi component in normal direction.
        """
        return np.array([
            [ 0,   nz, -ny],
            [-nz,  0,   nx],
            [ ny, -nx,  0 ],
        ])

    # ── Assembly [Au] matrix (4×24) ───────────────────────

    def _Au_matrix(self) -> np.ndarray:
        """
        [Au] maps nodal DOFs {un}(24) → Δβs_n (4×1) numerator
        before division by (2/3)(1+ϕk).
        Eq.(58-59) Katili 2015.
        """
        xyz = self._xyz
        normals = self._nodal_normals()
        Au = np.zeros((4, 24))

        for k, (i, j) in enumerate(self._side_pairs):
            # Unit tangent vector t_sk on side k
            xji = xyz[j] - xyz[i]
            Lk  = np.linalg.norm(xji)
            if Lk < 1e-15:
                continue
            tsk = xji / Lk

            # Average normal on side k
            nk = 0.5 * (normals[i] + normals[j])
            nm = np.linalg.norm(nk)
            if nm > 1e-15:
                nk /= nm

            # RN matrices at nodes i and j
            RNi = self._RN(*normals[i])
            RNj = self._RN(*normals[j])

            # Contribution from node i (cols 6i..6i+5)
            # From eq.(58): (1/Lk) n_k · (uj - ui)  + (1/2) t_sk·RN_i θ_i  + ...
            ci = i * 6
            cj = j * 6

            # Translation part (U,V,W): (1/Lk) * nk
            Au[k, ci:ci+3]     = -nk / Lk   # node i: -nk/Lk * (uj-ui) from ui side
            Au[k, cj:cj+3]     =  nk / Lk   # node j:  nk/Lk

            # Rotation part at node i (θX,θY,θZ):  (1/2) tsk · RN_i
            Au[k, ci+3:ci+6]   = 0.5 * (RNi.T @ tsk)

            # Rotation part at node j (θX,θY,θZ):  (1/2) tsk · RN_j
            Au[k, cj+3:cj+6]   = 0.5 * (RNj.T @ tsk)

        return Au

    # ── AΔ matrix (4×4 diagonal) ──────────────────────────

    def _Adelta_matrix(self) -> np.ndarray:
        """Diagonal [AΔ] = (2/3)(1+ϕk) per side. Eq.(60) Katili 2015."""
        phi = [self._phi_side(i, j) for i, j in self._side_pairs]
        return np.diag([2.0/3.0*(1.0 + phi[k]) for k in range(4)])

    # ── Bm: membrane strain matrix (3×24) ─────────────────

    def _Bm_at(self, xi: float, eta: float,
               t1, t2, n_vec, J, Co, bc_mat, normals) -> np.ndarray:
        """
        Membrane strain-displacement matrix [Bm] (3×24).
        e = {ex, ey, exy} = [Bm] {un}
        Eq.(25-26) Katili 2015.
        """
        dN = shape_dN(xi, eta)   # (2,4)

        # Ni,x = Ni,xi * Co11 + Ni,eta * Co21
        # Ni,y = Ni,xi * Co12 + Ni,eta * Co22
        Nix = dN[0]*Co[0,0] + dN[1]*Co[1,0]   # (4,)
        Niy = dN[0]*Co[0,1] + dN[1]*Co[1,1]   # (4,)

        Bm = np.zeros((3, 24))
        for i in range(4):
            col = i * 6
            # Eq.(26): [t1]·Ni,x for ex,  [t2]·Ni,y for ey,  mixed for exy
            # Translation DOFs (0,1,2) only for membrane
            for d in range(3):
                Bm[0, col+d] = t1[d] * Nix[i]
                Bm[1, col+d] = t2[d] * Niy[i]
                Bm[2, col+d] = t1[d] * Niy[i] + t2[d] * Nix[i]
        return Bm

    # ── Bb: bending strain matrix (3×24) ──────────────────

    def _Bb_at(self, xi: float, eta: float,
               t1, t2, n_vec, J, Co, bc_mat, normals,
               Adelta_inv_Au: np.ndarray) -> np.ndarray:
        """
        Bending curvature-displacement matrix [Bb] (3×24).
        χ = [Bb] {un}
        Eq.(62) Katili 2015: [Bb] = [Bbβ] + [BbΔβ][AΔ]^{-1}[Au]
        """
        dN = shape_dN(xi, eta)   # (2,4)
        dP = shape_dP(xi, eta)   # (2,4)

        # Local derivatives
        Nix = dN[0]*Co[0,0] + dN[1]*Co[1,0]
        Niy = dN[0]*Co[0,1] + dN[1]*Co[1,1]
        Pkx = dP[0]*Co[0,0] + dP[1]*Co[1,0]
        Pky = dP[0]*Co[0,1] + dP[1]*Co[1,1]

        # bc coupling terms (for non-planar shells)
        bc = bc_mat   # (2×2)
        # Nbc1i = Ni,xi*bc11 + Ni,eta*bc21
        # Nbc2i = Ni,xi*bc12 + Ni,eta*bc22
        Nbc1 = dN[0]*bc[0,0] + dN[1]*bc[1,0]   # (4,)
        Nbc2 = dN[0]*bc[0,1] + dN[1]*bc[1,1]   # (4,)

        # [Bbβ]: contribution from nodal DOFs directly (4 nodes × 6 DOF)
        Bbβ = np.zeros((3, 24))
        for i in range(4):
            col  = i * 6
            RNi  = self._RN(*normals[i])
            V1i  = np.array([t1 @ RNi[:,0], t1 @ RNi[:,1], t1 @ RNi[:,2]])
            V2i  = np.array([t2 @ RNi[:,0], t2 @ RNi[:,1], t2 @ RNi[:,2]])

            # Translation DOFs — coupling only if non-coplanar
            for d in range(3):
                Bbβ[0, col+d] = t1[d]*Nbc1[i]
                Bbβ[1, col+d] = t2[d]*Nbc2[i]
                Bbβ[2, col+d] = t1[d]*Nbc2[i] + t2[d]*Nbc1[i]

            # Rotation DOFs (θX,θY,θZ)
            # χx = ... + V1i · Ni,x
            Bbβ[0, col+3:col+6] = V1i * Nix[i]
            Bbβ[1, col+3:col+6] = V2i * Niy[i]
            Bbβ[2, col+3:col+6] = V1i * Niy[i] + V2i * Nix[i]

        # [BbΔβ]: contribution from side DOFs (4 sides × 1 DOF)
        # Pk,x and Pk,y for k=5..8
        BbDβ = np.zeros((3, 4))
        for k in range(4):
            i, j = self._side_pairs[k]
            xji  = self._xyz[j] - self._xyz[i]
            Lk   = np.linalg.norm(xji)
            if Lk < 1e-15:
                continue
            tsk = xji / Lk
            t1_tsk = t1 @ tsk
            t2_tsk = t2 @ tsk
            BbDβ[0, k] = t1_tsk * Pkx[k]
            BbDβ[1, k] = t2_tsk * Pky[k]
            BbDβ[2, k] = t1_tsk * Pky[k] + t2_tsk * Pkx[k]

        # Full [Bb] = [Bbβ] + [BbΔβ] [AΔ]^{-1}[Au]
        Bb = Bbβ + BbDβ @ Adelta_inv_Au
        return Bb

    # ── Bs: shear strain matrix (2×24) ────────────────────

    def _Bs_at(self, xi: float, eta: float,
               t1, t2, n_vec, J, Co, bc_mat, normals,
               Adelta_inv_Au: np.ndarray) -> np.ndarray:
        """
        Transverse shear strain-displacement matrix [Bs] (2×24).
        γ = [Bs] {un}
        Eq.(63) Katili 2015: [Bs] = [Bsγ] [Aϕ] [AΔ]^{-1}[Au]
        """
        # Nγ: shear interpolation (2×4), eq.(37)
        Ngamma = np.array([
            [0.5*(1-eta), 0, 0.5*(1+eta), 0],
            [0, 0.5*(1+xi), 0, 0.5*(1-xi)],
        ])

        # Aγ: side length factors (4×4 diagonal), eq.(38)
        Ls = np.zeros(4)
        for k, (i, j) in enumerate(self._side_pairs):
            Ls[k] = np.linalg.norm(self._xyz[j] - self._xyz[i])
        signs = [1, 1, -1, -1]   # sign convention from eq.(38)
        Ag = np.diag([signs[k]*Ls[k]/2 for k in range(4)])

        # Aϕ: shear factor matrix (4×4 diagonal), eq.(46)
        phi = [self._phi_side(i, j) for i, j in self._side_pairs]
        Aphi = np.diag([2.0/3.0*phi[k] for k in range(4)])

        # [Bsγ] = Co^T [Nγ] [Ag]  (2×4)
        Bsg = Co.T @ Ngamma @ Ag   # (2×4)

        # [Bs] = [Bsγ][Aϕ][AΔ]^{-1}[Au]  (2×24)
        Bs = Bsg @ Aphi @ Adelta_inv_Au
        return Bs

    # ── Constitutive matrices ─────────────────────────────

    def _Hm(self) -> np.ndarray:
        """Membrane constitutive matrix [Hm] (3×3). Eq.(66)."""
        Dm = self.E * self.h / (1 - self.nu**2)
        return Dm * np.array([
            [1,      self.nu, 0             ],
            [self.nu, 1,      0             ],
            [0,      0,      (1-self.nu)/2  ],
        ])

    def _Hb(self) -> np.ndarray:
        """Bending constitutive matrix [Hb] (3×3). Eq.(68)."""
        Db = self.E * self.h**3 / (12*(1 - self.nu**2))
        return Db * np.array([
            [1,      self.nu, 0             ],
            [self.nu, 1,      0             ],
            [0,      0,      (1-self.nu)/2  ],
        ])

    def _Hs(self) -> np.ndarray:
        """Shear constitutive matrix [Hs] (2×2). Eq.(70)."""
        Ds = self.kappa * self.E * self.h / (2*(1 + self.nu))
        return Ds * np.eye(2)

    # ── Fictitious stiffness kθz (MacNeal stabilization) ──

    def _k_fictitious(self) -> np.ndarray:
        """
        Fictitious stiffness for drilling DOF θZ to avoid spurious modes.
        Eqs.(73-76) Katili 2015.
        Factor = 10^{-3} (empirical, Katili 2015 §2.8).
        """
        alpha = 1e-3 * self.E * self.h**3 / 12.0
        normals = self._nodal_normals()

        kthz_grad = np.zeros((24, 24))
        kthz_mac  = np.zeros((24, 24))

        for xi, eta, w in GAUSS_POINTS:
            N  = shape_N(xi, eta)
            dN = shape_dN(xi, eta)
            t1, t2, n_vec, J, Co, bc_mat, _ = self._geometry_at(xi, eta)

            # θZ at integration point
            # θz = Σ Ni * [ni] · {θXi, θYi, θZi}
            # Gradient of θz in local x,y
            # dθz/dx = Σ Ni,x * ni·{...},  dθz/dy = Σ Ni,y * ni·{...}

            Nix = dN[0]*Co[0,0] + dN[1]*Co[1,0]
            Niy = dN[0]*Co[0,1] + dN[1]*Co[1,1]

            # Build (2×24) matrix [Gθz] s.t. {dθz/dx, dθz/dy} = [Gθz]{un}
            Gθz = np.zeros((2, 24))
            for i in range(4):
                col = i * 6
                ni  = normals[i]
                Gθz[0, col+3:col+6] = Nix[i] * ni
                Gθz[1, col+3:col+6] = Niy[i] * ni

            # Build (24,) vector [hθz] s.t. θz = [hθz]·{un}
            Hθz = np.zeros(24)
            for i in range(4):
                col = i * 6
                Hθz[col+3:col+6] = N[i] * normals[i]

            kthz_grad += w * J * (Gθz.T @ Gθz)
            kthz_mac  += w * J * np.outer(Hθz, Hθz)

        kthz_grad *= alpha
        kthz_mac  *= 1e-3 * self.E * self.h / (2*(1+self.nu))

        return kthz_grad + kthz_mac

    # ── Main stiffness assembly ───────────────────────────

    def _compute_stiffness(self) -> np.ndarray:
        """
        Assemble element stiffness matrix K (24×24).
        K = km + kb + ks + kθz
        Eq.(77) Katili 2015.
        """
        Hm = self._Hm()
        Hb = self._Hb()
        Hs = self._Hs()

        # Precompute [AΔ]^{-1} [Au]  (4×24) — constant for this element
        Au     = self._Au_matrix()
        Adelta = self._Adelta_matrix()
        Adelta_inv = np.diag(1.0 / np.diag(Adelta))
        Adelta_inv_Au = Adelta_inv @ Au    # (4×24)

        K = np.zeros((24, 24))

        for xi, eta, w in GAUSS_POINTS:
            t1, t2, n_vec, J, Co, bc_mat, normals = self._geometry_at(xi, eta)

            Bm = self._Bm_at(xi, eta, t1, t2, n_vec, J, Co, bc_mat, normals)
            Bb = self._Bb_at(xi, eta, t1, t2, n_vec, J, Co, bc_mat, normals, Adelta_inv_Au)
            Bs = self._Bs_at(xi, eta, t1, t2, n_vec, J, Co, bc_mat, normals, Adelta_inv_Au)

            K += w * J * (Bm.T @ Hm @ Bm +
                          Bb.T @ Hb @ Bb +
                          Bs.T @ Hs @ Bs)

        # Add fictitious stiffness
        K += self._k_fictitious()

        return K

    # ── Public: recompute after averaged normals are injected ──

    def recompute_stiffness(self) -> None:
        """
        Re-run stiffness assembly after _averaged_normals has been set.
        Must be called after compute_averaged_nodal_normals() has patched
        this element.
        """
        self._K24 = self._compute_stiffness()


# ─────────────────────────────────────────────────────────
#  MESH-LEVEL UTILITY: Assembly-averaged nodal normals
# ─────────────────────────────────────────────────────────

def compute_averaged_nodal_normals(elements: list) -> None:
    """
    Two-pass assembly-averaging of nodal normals across all DKMQ24 elements.

    Pass 1 — accumulate weighted face normals at every node:
        For each element, compute the element centroid normal (cross product
        of diagonals, area-weighted) and add it to each of its 4 nodes.
        Area-weighting means larger elements contribute proportionally more,
        which is the standard angle/area-weighted averaging used in mesh
        processing and FEM pre-processors.

    Pass 2 — normalise and inject into each element:
        Each element's _averaged_normals (4×3) is filled from the global
        node table, then unit-normalised.  Elements call recompute_stiffness()
        to rebuild K with the improved normals.

    Parameters
    ----------
    elements : list[DKMQ24]
        All shell elements in the mesh.  Elements must share Node objects
        (same Python object, identified by id()) for the connectivity to
        be detected correctly.

    Notes
    -----
    - Pure-Python / NumPy; no external dependencies.
    - O(N_elem) time, O(N_node) memory beyond the element objects.
    - Safe to call multiple times (idempotent): re-averages from scratch.
    - For boundary nodes that appear in only one element the averaged
      normal equals the local single-element estimate, so the fallback
      path is numerically unchanged.
    - Sign consistency: all element normals are oriented by the element's
      own a1×a2 (CCW node ordering), so they should all point to the
      same side of the shell.  If your mesh has mixed orientations this
      step will produce incorrect results — fix element orientations first.
    """
    from collections import defaultdict

    # ── Pass 1: accumulate area-weighted face normals per node id ──
    node_normal_acc = defaultdict(lambda: np.zeros(3))   # node_id → sum of weighted normals

    for elem in elements:
        xyz = elem._xyz   # (4×3)

        # Element centroid normal via cross product of diagonals.
        # d1 = node2 - node0,  d2 = node3 - node1
        # This is more robust than using edge pairs at each corner.
        d1 = xyz[2] - xyz[0]
        d2 = xyz[3] - xyz[1]
        face_n = np.cross(d1, d2)
        area_weight = np.linalg.norm(face_n)   # ≈ 2 × element area
        if area_weight < 1e-15:
            continue
        face_n_unit = face_n / area_weight

        # Additionally compute per-corner normals from element edges
        # so that strongly non-planar quads don't introduce a single
        # face normal that poorly represents the corner neighbourhood.
        corner_n = np.zeros((4, 3))
        local_pairs = [
            (xyz[1] - xyz[0], xyz[3] - xyz[0]),
            (xyz[2] - xyz[1], xyz[0] - xyz[1]),
            (xyz[3] - xyz[2], xyz[1] - xyz[2]),
            (xyz[0] - xyz[3], xyz[2] - xyz[3]),
        ]
        for c, (a, b) in enumerate(local_pairs):
            cn = np.cross(a, b)
            cnm = np.linalg.norm(cn)
            corner_n[c] = cn / cnm if cnm > 1e-15 else face_n_unit

        # Weighted accumulation: use corner-local normal weighted by
        # the subtended solid angle (approximated by the corner triangle area).
        # Simpler area-weight: weight each corner equally by element area.
        for local_idx, node in enumerate(elem.nodes):
            nid = id(node)
            node_normal_acc[nid] += area_weight * corner_n[local_idx]

    # ── Pass 2: normalise and inject ──────────────────────
    # Build a quick lookup: node python-id → averaged unit normal
    node_normal_final = {}
    for nid, acc in node_normal_acc.items():
        nm = np.linalg.norm(acc)
        node_normal_final[nid] = acc / nm if nm > 1e-15 else np.array([0., 0., 1.])

    # Patch each element and recompute stiffness
    for elem in elements:
        avg = np.zeros((4, 3))
        for local_idx, node in enumerate(elem.nodes):
            avg[local_idx] = node_normal_final[id(node)]
        elem._averaged_normals = avg
        elem.recompute_stiffness()


# ─────────────────────────────────────────────────────────
#  Public aliases and direct normal injection utility
# ─────────────────────────────────────────────────────────

DKMQ24Katili = DKMQ24


def inject_nodal_normals(elements: list, normal_by_node: dict) -> None:
    """
    Inject externally supplied nodal normals into DKMQ24 elements and recompute.

    Parameters
    ----------
    elements:
        List of DKMQ24 elements.
    normal_by_node:
        Dictionary keyed either by node object, id(node), or node.nid.
        Each value is a 3-vector unit normal or non-unit normal.

    Notes
    -----
    This is useful for analytic SNORM tests such as Scordelis-Lo, where the
    cylindrical nodal normal is known exactly.  For general meshes use
    compute_averaged_nodal_normals(elements).
    """
    for elem in elements:
        avg = np.zeros((4, 3), dtype=float)
        for i, node in enumerate(elem.nodes):
            if node in normal_by_node:
                n = np.asarray(normal_by_node[node], dtype=float)
            elif id(node) in normal_by_node:
                n = np.asarray(normal_by_node[id(node)], dtype=float)
            elif hasattr(node, "nid") and node.nid in normal_by_node:
                n = np.asarray(normal_by_node[node.nid], dtype=float)
            else:
                raise KeyError("normal_by_node missing one or more element nodes")
            nm = np.linalg.norm(n)
            avg[i] = n / nm if nm > 1e-15 else np.array([0.0, 0.0, 1.0])
        elem._averaged_normals = avg
        elem.recompute_stiffness()


__all__ = [
    "DKMQ24",
    "DKMQ24Katili",
    "shape_N",
    "shape_dN",
    "shape_P",
    "shape_dP",
    "GAUSS_POINTS",
    "compute_averaged_nodal_normals",
    "inject_nodal_normals",
]

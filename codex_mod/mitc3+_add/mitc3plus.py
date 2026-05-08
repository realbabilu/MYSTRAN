"""
mitc3plus.py — MITC3+ Triangular Shell Element (Production Implementation)
============================================================================
Referensi utama:
  Lee, Y., Lee, P.-S., Bathe, K.-J. (2014).
  "The MITC3+ shell element and its performance."
  Computers and Structures, 138, 12–23.

Formulasi:
  - 3-node triangular shell, 6 DOF/node (u,v,w,rx,ry,rz)
  - Cubic bubble function untuk enrichment rotasi (internal DOF: α4, β4)
  - New assumed transverse shear strain field (MITC scheme)
  - Static condensation bubble DOF pada level elemen
  - AutoSPC drilling via penalty untuk DOF rz (normal rotation)

Catatan implementasi:
  - Geometry flat per elemen (director V_n dari cross-product sisi elemen)
  - Local frame: e1 // sisi 1-2, e3 = V_n, e2 = e3 × e1
  - Integrasi: 7-point Gauss pada domain triangular (optimal per paper)
  - d = 1/10000 untuk tying points (D,E,F) — zero-energy-free per paper
  - Semua notasi konsisten dengan paper Sec. 2.2
"""

from __future__ import annotations

import numpy as np
from numpy.linalg import norm, inv
from typing import List

from core import Element, Node


# ══════════════════════════════════════════════════════════════════════════════
#  KONSTANTA INTEGRASI
# ══════════════════════════════════════════════════════════════════════════════

# 7-point Gauss quadrature pada domain segitiga (r,s) ∈ triangle
# Ref: Bathe (1996) Table 5.5
_GAUSS_7PT: List = None

def _gauss_triangle_7pt():
    """Return (r, s, w) arrays for 7-point Gauss quadrature on unit triangle."""
    global _GAUSS_7PT
    if _GAUSS_7PT is not None:
        return _GAUSS_7PT
    a1 = 0.1012865073235
    b1 = 0.7974269853531
    a2 = 0.4701420641051
    b2 = 0.0597158717898
    c  = 1.0/3.0
    w1 = 0.1259391805448
    w2 = 0.1323941527885
    w3 = 0.2250000000000
    r = np.array([a1, b1, a1, a2, b2, a2, c])
    s = np.array([a1, a1, b1, a2, a2, b2, c])
    w = np.array([w1, w1, w1, w2, w2, w2, w3]) * 0.5  # ×0.5 karena area segitiga = 0.5
    _GAUSS_7PT = (r, s, w)
    return _GAUSS_7PT


# ══════════════════════════════════════════════════════════════════════════════
#  INTERPOLATION FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def _h(r: float, s: float):
    """Linear shape functions h1, h2, h3."""
    return np.array([1.0 - r - s, r, s])

def _dh():
    """Derivatives ∂h/∂r, ∂h/∂s → shape (3,2)."""
    return np.array([
        [-1.0, -1.0],
        [ 1.0,  0.0],
        [ 0.0,  1.0],
    ])

def _f4(r: float, s: float) -> float:
    """Cubic bubble function f4 = 27·r·s·(1-r-s)."""
    return 27.0 * r * s * (1.0 - r - s)

def _df4(r: float, s: float):
    """Derivatives ∂f4/∂r, ∂f4/∂s."""
    t = 1.0 - r - s
    df4dr = 27.0 * s * (t - r)
    df4ds = 27.0 * r * (t - s)
    return np.array([df4dr, df4ds])

def _fi(r: float, s: float):
    """fi = hi - (1/3)·f4 for i=1,2,3; f4 itself for bubble."""
    f4v  = _f4(r, s)
    hv   = _h(r, s)
    return hv - f4v / 3.0, f4v

def _dfi(r: float, s: float):
    """∂fi/∂ξ for i=1..3, plus ∂f4/∂ξ. Returns shape (4,2)."""
    dh  = _dh()
    df4 = _df4(r, s)
    dfi123 = dh - df4[np.newaxis, :] / 3.0
    return np.vstack([dfi123, df4])  # (4,2)


# ══════════════════════════════════════════════════════════════════════════════
#  LOCAL FRAME
# ══════════════════════════════════════════════════════════════════════════════

def _local_frame(nodes: List[Node]):
    """
    Compute element local frame (e1, e2, e3) dan director V_n per node.

    Returns:
      e1, e2, e3 : unit vectors, global Cartesian
      V_n        : normal director (same for flat element)
      x_local    : node coords in local frame (3,3)
    """
    x1 = nodes[0].coords
    x2 = nodes[1].coords
    x3 = nodes[2].coords

    v12 = x2 - x1
    v13 = x3 - x1

    e1  = v12 / norm(v12)
    e3  = np.cross(v12, v13)
    e3_norm = norm(e3)
    if e3_norm < 1e-14:
        raise ValueError(f"MITC3+: degenerate element, nodes collinear")
    e3  /= e3_norm
    e2  = np.cross(e3, e1)
    e2  /= norm(e2)

    # Local coords (project to e1,e2 plane)
    R = np.vstack([e1, e2, e3])  # 3×3 rotation (rows = local axes)
    x_local = np.array([(R @ (n.coords - x1)) for n in nodes])

    return e1, e2, e3, x_local


# ══════════════════════════════════════════════════════════════════════════════
#  JACOBIAN (isoparametric mapping pada mid-surface)
# ══════════════════════════════════════════════════════════════════════════════

def _jacobian_2d(r: float, s: float, xy: np.ndarray):
    """
    Jacobian J = ∂(x,y)/∂(r,s) pada koordinat lokal 2D.
    xy : (3,2) koordinat lokal x,y dari 3 nodes.
    Returns J (2×2), det_J.
    """
    dh = _dh()  # (3,2)
    J  = dh.T @ xy  # (2,2)
    return J, np.linalg.det(J)


# ══════════════════════════════════════════════════════════════════════════════
#  MATERIAL MATRIX (plane stress + transverse shear)
# ══════════════════════════════════════════════════════════════════════════════

def _C_membrane(E: float, nu: float, h: float) -> np.ndarray:
    """Membrane constitutive matrix Cm (3×3), integrated through thickness."""
    fac = E * h / (1.0 - nu**2)
    return fac * np.array([
        [1.0,  nu,       0.0          ],
        [nu,   1.0,      0.0          ],
        [0.0,  0.0,  (1.0-nu)/2.0    ],
    ])

def _C_bending(E: float, nu: float, h: float) -> np.ndarray:
    """Bending constitutive matrix Cb (3×3), integrated through thickness."""
    fac = E * h**3 / (12.0 * (1.0 - nu**2))
    return fac * np.array([
        [1.0,  nu,       0.0          ],
        [nu,   1.0,      0.0          ],
        [0.0,  0.0,  (1.0-nu)/2.0    ],
    ])

def _C_shear(E: float, nu: float, h: float, kappa: float = 5.0/6.0) -> np.ndarray:
    """Transverse shear constitutive matrix Cs (2×2)."""
    G = E / (2.0 * (1.0 + nu))
    fac = kappa * G * h
    return fac * np.eye(2)


# ══════════════════════════════════════════════════════════════════════════════
#  B-MATRIX MEMBRANE
# ══════════════════════════════════════════════════════════════════════════════

def _Bm_at(r: float, s: float, J_inv: np.ndarray, xy: np.ndarray):
    """
    Membrane B-matrix (3 × 3·3=9) untuk 3 nodes, DOF = [u1,v1, u2,v2, u3,v3].
    Strain = [ε_xx, ε_yy, 2ε_xy]^T
    """
    dh_rs = _dh()          # (3,2) ∂h/∂(r,s)
    dh_xy = dh_rs @ J_inv.T  # (3,2) ∂h/∂(x,y)

    Bm = np.zeros((3, 6))
    for i in range(3):
        Bm[0, 2*i  ] = dh_xy[i, 0]   # ∂u/∂x
        Bm[1, 2*i+1] = dh_xy[i, 1]   # ∂v/∂y
        Bm[2, 2*i  ] = dh_xy[i, 1]   # ∂u/∂y
        Bm[2, 2*i+1] = dh_xy[i, 0]   # ∂v/∂x
    return Bm


# ══════════════════════════════════════════════════════════════════════════════
#  B-MATRIX BENDING (termasuk bubble enrichment)
# ══════════════════════════════════════════════════════════════════════════════

def _Bb_at(r: float, s: float, J_inv: np.ndarray):
    """
    Bending B-matrix (3 × (3·2 + 2)) = (3 × 8).
    DOF order: [α1,β1, α2,β2, α3,β3, α4,β4]
    Curvature κ = [κ_xx, κ_yy, 2κ_xy]^T
    αi, βi = rotations about V1^i, V2^i (local axes at node i)

    Untuk elemen flat: V1^i = e1, V2^i = e2 untuk semua nodes.
    κ_xx = ∂β/∂x,  κ_yy = -∂α/∂y,  κ_xy = (∂β/∂y - ∂α/∂x)/2
    (sign convention sesuai Bathe 1996)
    """
    dfi_rs = _dfi(r, s)       # (4,2) — rows: nodes 1..3 + bubble
    dfi_xy = dfi_rs @ J_inv.T  # (4,2) ∂fi/∂(x,y)

    Bb = np.zeros((3, 8))
    for i in range(4):  # nodes 1,2,3 + bubble
        # αi: rotation about e1 → affects κ_yy, κ_xy
        # βi: rotation about e2 → affects κ_xx, κ_xy
        # κ_xx = -∂α_z/∂x + ... for Mindlin. Using Bathe convention:
        # κ_xx =  ∂β/∂x
        # κ_yy = -∂α/∂y
        # 2κ_xy = ∂β/∂y - ∂α/∂x
        Bb[0, 2*i+1] =  dfi_xy[i, 0]   # ∂βi/∂x
        Bb[1, 2*i  ] = -dfi_xy[i, 1]   # -∂αi/∂y
        Bb[2, 2*i+1] =  dfi_xy[i, 1]   # ∂βi/∂y
        Bb[2, 2*i  ] = -dfi_xy[i, 0]   # -∂αi/∂x
    return Bb


# ══════════════════════════════════════════════════════════════════════════════
#  ASSUMED TRANSVERSE SHEAR STRAIN — MITC3+ (Eq. 15-17 paper)
# ══════════════════════════════════════════════════════════════════════════════

# Tying points per Table 1 of paper
_TYP_ABC = np.array([
    [1.0/6, 2.0/3],   # (A)
    [2.0/3, 1.0/6],   # (B)
    [1.0/6, 1.0/6],   # (C)
])

_D_PARAM = 1.0 / 10000.0  # paper recommends d=1/10000

def _tying_DEF(d: float = _D_PARAM):
    """Tying points (D,E,F) on internal lines, barycenter-to-edge-midpoint."""
    # Midpoints of edges:
    #   edge 1-2: r=0.5, s=0
    #   edge 2-3: r=0.5, s=0.5
    #   edge 3-1: r=0,   s=0.5
    # Internal lines from barycenter (1/3,1/3) toward these midpoints:
    # D: toward (0.5, 0)   → (1/3 + d_r, 1/3 + d_s)
    # per paper Fig 5(c) and Table 1:
    return np.array([
        [1.0/3 + d,   1.0/3 - 2*d],  # (D)
        [1.0/3 - 2*d, 1.0/3 + d  ],  # (E)
        [1.0/3 + d,   1.0/3 + d  ],  # (F)
    ])


def _transverse_shear_disp(r: float, s: float,
                            h_local: float,
                            J_inv: np.ndarray,
                            fi_vals: np.ndarray,
                            dfi_xy: np.ndarray) -> np.ndarray:
    """
    Compute displacement-based covariant transverse shear strains [e_rt, e_st]
    at (r,s). DOF vector order: [u,v,w,α,β] × 3nodes + [α4,β4] bubble.
    Returns function (gdof_vec → [e_rt, e_st]) — BUT here we return
    the B_s sub-matrices that multiply [w1,α1,β1, w2,α2,β2, w3,α3,β3, α4,β4].

    For flat element, transverse shear strains (covariant in local frame):
      e_rt = ∂w/∂r + (t/2)·∂(fi·βi)/∂r·... simplified for thin shell (t→0):
      e_rt = Σ ∂hi/∂r · wi + (a/2)·Σ fi · βi_via_covariant...

    We use the Reissner-Mindlin kinematics in the element local frame:
      γ_xz = ∂w/∂x - α  (where α = rotation about e1)
      γ_yz = ∂w/∂y + β  (where β = rotation about e2)

    Covariant shear: e_rt = ∂w/∂r - Σfi·αi (chain-rule from x to r)
                    e_st = ∂w/∂s + Σfi·βi
    More precisely using the full Jacobian mapping.
    """
    # We return 2×(11) B-matrix rows for [w1,α1,β1, w2,α2,β2, w3,α3,β3, α4,β4]
    # Actually core DOF per node in local = [u,v,w,α,β,—] but we only care
    # about w and α,β rotations here.
    pass


def _compute_ert_est_dispbased(r: float, s: float, xy: np.ndarray,
                                J: np.ndarray, J_inv: np.ndarray,
                                dof_vec: np.ndarray) -> np.ndarray:
    """
    Evaluate displacement-based covariant transverse shear at (r,s).
    dof_vec: [w1,α1,β1, w2,α2,β2, w3,α3,β3, α4,β4] — 11 entries

    Covariant strains in parametric coords:
      e_rt = ∂w/∂r + (h/2)·[Σ(∂fi/∂r)·(-αi)]   (Mindlin w/ sign per Bathe)
      e_st = ∂w/∂s + (h/2)·[Σ(∂fi/∂s)·(-βi... )]

    For thin limit (h→0 in rotational terms irrelevant here):
      e_rt = Σ ∂hi/∂r · wi  +  Σ fi · (J[0,0]·βi - J[0,1]·(-αi)) ... 
    
    SIMPLIFICATION for plate/shell with unit thickness director:
      Using local frame where e3 ⊥ mid-surface:
      γ_s1 = ∂w/∂s1 - α  (s1 = x local)
      γ_s2 = ∂w/∂s2 + β  (s2 = y local)
      
      Covariant:
      e_rt = J[0,0]·γ_x + J[0,1]·γ_y  (chain rule: ∂/∂r = J[0,0]·∂/∂x + J[0,1]·∂/∂y)
    """
    # Decode dof_vec
    w  = dof_vec[0:9:3]   # w1,w2,w3
    al = dof_vec[1:9:3]   # α1,α2,α3
    be = dof_vec[2:9:3]   # β1,β2,β3
    al4 = dof_vec[9]
    be4 = dof_vec[10]

    al_all = np.append(al, al4)
    be_all = np.append(be, be4)

    # ∂w/∂r, ∂w/∂s using linear shape functions
    dh_rs = _dh()  # (3,2)
    dw_dr = dh_rs[:, 0] @ w
    dw_ds = dh_rs[:, 1] @ w

    # ∂w/∂x = J_inv[0,0]·∂w/∂r + J_inv[0,1]·∂w/∂s
    dw_dx = J_inv[0, 0] * dw_dr + J_inv[0, 1] * dw_ds
    dw_dy = J_inv[1, 0] * dw_dr + J_inv[1, 1] * dw_ds

    # α_interp = Σ fi·αi (fi for nodes 1..3 + bubble)
    fi3, f4v = _fi(r, s)
    fi_all = np.append(fi3, f4v)

    a_interp = fi_all @ al_all
    b_interp = fi_all @ be_all

    # Transverse shear in local frame (Mindlin convention):
    # γ_xz = ∂w/∂x - α_about_y = ∂w/∂x + β  (β = rotation about e2 = about local-y)
    # Wait — sign convention: for Reissner-Mindlin plate,
    # γ_xz = ∂w/∂x - θ_y,  γ_yz = ∂w/∂y + θ_x
    # Here αi = θ_x-like (rotation about e1=x), βi = θ_y-like (rotation about e2=y)
    # So: γ_xz = ∂w/∂x - β_interp,  γ_yz = ∂w/∂y + α_interp
    gxz = dw_dx - b_interp
    gyz = dw_dy + a_interp

    # Covariant: e_rt = J[0,0]·γ_xz + J[0,1]·γ_yz
    #            e_st = J[1,0]·γ_xz + J[1,1]·γ_yz
    e_rt = J[0, 0] * gxz + J[0, 1] * gyz
    e_st = J[1, 0] * gxz + J[1, 1] * gyz

    return np.array([e_rt, e_st])


def _Bs_mitc3plus(xy: np.ndarray, h_thickness: float) -> np.ndarray:
    """
    Build MITC3+ assumed transverse shear B-matrix (2 × 11) per paper Eq.17.
    xy   : (3,2) node local coordinates
    Returns Bs such that [e_rt_assumed, e_st_assumed] = Bs @ dof_vec_s
    where dof_vec_s = [w1,α1,β1, w2,α2,β2, w3,α3,β3, α4,β4] (11 entries)

    Strategy:
      1. Evaluate displacement-based e_rt,e_st at tying points A,B,C,D,E,F
         using linearized B-matrices at those points
      2. Assemble MITC3+ assumed field per Eq.(15)-(17)
      3. Result is a (2×11) operator that when applied gives
         [^e_rt, ^e_st] at any (r,s) — but we actually need full
         (2 × 11 × Ngauss) for assembly. We return the operator
         that MAPS the 11-dof to [^e_rt(r,s), ^e_st(r,s)].

    NOTE: The MITC3+ assumed field is linear in (r,s), so:
      ^e_rt(r,s) = Bs_rt_row @ dof_s   where Bs_rt_row is linear in (r,s)
    We parameterize: Bs(r,s) = Bs0 + r·Bs_r + s·Bs_s
    by evaluating at 3 points and solving.

    SIMPLER DIRECT APPROACH: Build Bs row-by-row using the tying operators.
    """
    # Helper: B-vector (1×11) for e_rt or e_st at point (r,s) from displacement
    def _b_row_ert(r, s):
        J = _dh().T @ xy
        detJ = np.linalg.det(J)
        J_inv = inv(J)
        dh_rs = _dh()
        dh_xy = dh_rs @ J_inv.T
        fi3, f4v = _fi(r, s)
        fi_all = np.append(fi3, f4v)
        dfi_rs = _dfi(r, s)
        dfi_xy = dfi_rs @ J_inv.T

        # γ_xz = ∂w/∂x - β_interp → contributes to e_rt
        # γ_yz = ∂w/∂y + α_interp → contributes to e_rt
        # e_rt = J[0,0]*γ_xz + J[0,1]*γ_yz
        row = np.zeros(11)
        for i in range(3):
            row[3*i  ] = J[0, 0]*dh_xy[i, 0] + J[0, 1]*dh_xy[i, 1]  # ∂w_i
            row[3*i+1] = J[0, 1]*fi_all[i]                            # α_i: from γ_yz
            row[3*i+2] = -J[0, 0]*fi_all[i]                           # β_i: from γ_xz
        # bubble: α4, β4
        row[9 ] = J[0, 1]*fi_all[3]   # α4
        row[10] = -J[0, 0]*fi_all[3]  # β4
        return row

    def _b_row_est(r, s):
        J = _dh().T @ xy
        J_inv = inv(J)
        dh_rs = _dh()
        dh_xy = dh_rs @ J_inv.T
        fi3, f4v = _fi(r, s)
        fi_all = np.append(fi3, f4v)

        # e_st = J[1,0]*γ_xz + J[1,1]*γ_yz
        row = np.zeros(11)
        for i in range(3):
            row[3*i  ] = J[1, 0]*dh_xy[i, 0] + J[1, 1]*dh_xy[i, 1]
            row[3*i+1] = J[1, 1]*fi_all[i]
            row[3*i+2] = -J[1, 0]*fi_all[i]
        row[9 ] = J[1, 1]*fi_all[3]
        row[10] = -J[1, 0]*fi_all[3]
        return row

    A_pt = _TYP_ABC[0]
    B_pt = _TYP_ABC[1]
    C_pt = _TYP_ABC[2]
    DEF  = _tying_DEF()
    D_pt = DEF[0]
    E_pt = DEF[1]
    F_pt = DEF[2]

    # Evaluate displacement-based rows at tying points
    ert_A = _b_row_ert(*A_pt)
    est_A = _b_row_est(*A_pt)
    ert_B = _b_row_ert(*B_pt)
    est_B = _b_row_est(*B_pt)
    ert_C = _b_row_ert(*C_pt)
    est_C = _b_row_est(*C_pt)
    ert_D = _b_row_ert(*D_pt)
    est_D = _b_row_est(*D_pt)
    ert_E = _b_row_ert(*E_pt)
    est_E = _b_row_est(*E_pt)
    ert_F = _b_row_ert(*F_pt)
    est_F = _b_row_est(*F_pt)

    # MITC3+ assumed constant part (Eq. 15):
    # ^const_e_rt = (2/3)·[e_rt(B) - (1/2)·e_st(B)] + (1/3)·[e_rt(C) + e_st(C)]
    # ^const_e_st = (2/3)·[e_st(A) - (1/2)·e_rt(A)] + (1/3)·[e_rt(C) + e_st(C)]
    const_ert = (2.0/3.0)*(ert_B - 0.5*est_B) + (1.0/3.0)*(ert_C + est_C)
    const_est = (2.0/3.0)*(est_A - 0.5*ert_A) + (1.0/3.0)*(ert_C + est_C)

    # ĉ = e_rt(F) - e_rt(D) - [e_st(F) - e_st(E)]  (Eq. 16)
    c_hat = (ert_F - ert_D) - (est_F - est_E)

    # Linear part (Eq. 16-17):
    # ^linear_e_rt(r,s) = (1/3)·ĉ·(3s-1)
    # ^linear_e_st(r,s) = (1/3)·ĉ·(1-3r)
    # These are functions of (r,s), so the full Bs depends on (r,s):
    # ^e_rt(r,s) = const_ert + (1/3)·(3s-1)·c_hat @ dof_s
    # ^e_st(r,s) = const_est + (1/3)·(1-3r)·c_hat @ dof_s

    # Return as callable (r,s) → (2×11) matrix
    def Bs_at(r, s):
        fac_rt = (3.0*s - 1.0) / 3.0
        fac_st = (1.0 - 3.0*r) / 3.0
        Bs = np.zeros((2, 11))
        Bs[0, :] = const_ert + fac_rt * c_hat
        Bs[1, :] = const_est + fac_st * c_hat
        return Bs

    return Bs_at


# ══════════════════════════════════════════════════════════════════════════════
#  FULL ELEMENT STIFFNESS WITH STATIC CONDENSATION
# ══════════════════════════════════════════════════════════════════════════════

def _build_K_full(nodes: List[Node], E: float, nu: float, h: float,
                  drilling_penalty: float) -> np.ndarray:
    """
    Build full stiffness matrix (20×20) in local frame BEFORE condensation.
    DOF order (local, per node): [u, v, w, α, β, rz_drill] × 3
               + bubble: [α4, β4]
    Total = 18 + 2 = 20

    After condensation of [α4, β4], result is 18×18.
    Then map to global 18×18 (6DOF × 3 nodes).
    """
    e1, e2, e3, x_local = _local_frame(nodes)

    # 2D local node coords for Jacobian
    xy = x_local[:, :2]  # (3,2)

    J = _dh().T @ xy       # (2,2) Jacobian (constant for linear tri)
    det_J = np.linalg.det(J)
    if det_J < 1e-14:
        raise ValueError(f"MITC3+: near-zero Jacobian det={det_J:.3e}")
    J_inv = inv(J)

    # Constitutive matrices
    Cm = _C_membrane(E, nu, h)
    Cb = _C_bending(E, nu, h)
    Cs = _C_shear(E, nu, h)

    # MITC3+ shear B-matrix factory
    Bs_factory = _Bs_mitc3plus(xy, h)

    # Gauss integration
    r_gp, s_gp, w_gp = _gauss_triangle_7pt()

    # K dimensions:
    #   membrane: 3 nodes × (u,v) = 6 DOF
    #   bending:  3 nodes × (α,β) + bubble (α4,β4) = 8 DOF
    #   transverse shear: same 8 DOF (w + α,β per node)
    # We build sub-blocks then assemble.

    # Sub-stiffness: index mapping in 20-DOF local vector
    # Local DOF per node i: [u_i, v_i, w_i, α_i, β_i, rz_i]  → indices 6i..6i+5
    # Bubble: [α4, β4] → indices 18, 19
    # For membrane: DOF u_i, v_i → 6i, 6i+1
    # For bending:  DOF α_i, β_i → 6i+3, 6i+4  + bubble 18,19
    # For shear:    DOF w_i, α_i, β_i → 6i+2, 6i+3, 6i+4 + bubble 18,19

    N_full = 20
    K_full = np.zeros((N_full, N_full))

    # Membrane DOF indices: [u1,v1,u2,v2,u3,v3] → [0,1,6,7,12,13]
    idx_m = [0, 1, 6, 7, 12, 13]

    # Bending DOF indices (matches Bb convention [α1,β1,...,α4,β4]):
    # [α1,β1,α2,β2,α3,β3,α4,β4] → [3,4,9,10,15,16,18,19]
    idx_b = [3, 4, 9, 10, 15, 16, 18, 19]

    # Shear DOF [w1,α1,β1, w2,α2,β2, w3,α3,β3, α4,β4]:
    idx_s = [2, 3, 4, 8, 9, 10, 14, 15, 16, 18, 19]

    for i_gp in range(len(r_gp)):
        r = r_gp[i_gp]
        s = s_gp[i_gp]
        w = w_gp[i_gp]

        fac = w * det_J

        # Membrane
        Bm = _Bm_at(r, s, J_inv, xy)  # (3×6)
        Km_gp = Bm.T @ Cm @ Bm * fac
        for ii, gi in enumerate(idx_m):
            for jj, gj in enumerate(idx_m):
                K_full[gi, gj] += Km_gp[ii, jj]

        # Bending
        Bb = _Bb_at(r, s, J_inv)  # (3×8)
        Kb_gp = Bb.T @ Cb @ Bb * fac
        for ii, gi in enumerate(idx_b):
            for jj, gj in enumerate(idx_b):
                K_full[gi, gj] += Kb_gp[ii, jj]

        # MITC3+ transverse shear
        Bs = Bs_factory(r, s)  # (2×11)
        # Cs is in physical (x,y) — but Bs gives covariant (r,s) strains
        # Need to transform: σ_phys = J^{-T} · σ_cov (for stress-strain in cov frame)
        # Properly: shear strain energy = ∫ [e_rt,e_st]^T · Cs_cov · [e_rt,e_st] dA
        # where Cs_cov = J · Cs_phys · J^T / det_J  (covariant metric)
        # For simplicity (and per standard MITC practice for flat elements):
        # Use physical shear with Bs transformed to physical:
        #   [γ_xz, γ_yz] = J_inv @ [e_rt, e_st]   ← only when |J| uniform
        # Actually: e_rt = J[0,0]γ_xz + J[0,1]γ_yz, so [e] = J @ [γ]
        # Hence: [γ] = J_inv @ [e], and energy = γ^T Cs γ = e^T J_inv^T Cs J_inv e
        Cs_cov = J_inv.T @ Cs @ J_inv  # (2×2)
        Ks_gp  = Bs.T @ Cs_cov @ Bs * fac
        for ii, gi in enumerate(idx_s):
            for jj, gj in enumerate(idx_s):
                K_full[gi, gj] += Ks_gp[ii, jj]

    # Drilling penalty (very small, just to prevent singular K)
    # DOF rz_i → indices 5, 11, 17
    for idx_rz in [5, 11, 17]:
        K_full[idx_rz, idx_rz] += drilling_penalty

    return K_full


def _static_condense(K_full: np.ndarray) -> np.ndarray:
    """
    Statically condense out bubble DOF (indices 18 and 19).
    Returns K_condensed (18×18).
    K = [[Kaa, Kab], [Kba, Kbb]]
    K_cond = Kaa - Kab · Kbb^{-1} · Kba
    a = 0..17, b = 18..19
    """
    Kaa = K_full[:18, :18]
    Kab = K_full[:18, 18:]
    Kba = K_full[18:, :18]
    Kbb = K_full[18:, 18:]

    try:
        Kbb_inv = inv(Kbb)
    except np.linalg.LinAlgError:
        # Fallback: pseudo-inverse
        Kbb_inv = np.linalg.pinv(Kbb)

    return Kaa - Kab @ Kbb_inv @ Kba


# ══════════════════════════════════════════════════════════════════════════════
#  TRANSFORMATION MATRIX: LOCAL ↔ GLOBAL
# ══════════════════════════════════════════════════════════════════════════════

def _T_global(nodes: List[Node]) -> np.ndarray:
    """
    Build 18×18 transformation matrix from local to global DOF.
    Each node has 6 DOF: [ux,uy,uz, rx,ry,rz] in global,
                         [u1,u2,u3, α, β, rz] in local (e1,e2,e3 frame).
    The transformation for each node is a 6×6 block:
      [u_glob] = R^T [u_loc]
    where R = [e1|e2|e3] (columns).
    Drilling DOF (local rz = about e3) maps to global rz = about global e3.
    """
    e1, e2, e3, _ = _local_frame(nodes)
    R3 = np.column_stack([e1, e2, e3])  # 3×3 rotation matrix

    # For displacements AND rotations, same R3
    T6 = np.zeros((6, 6))
    T6[:3, :3] = R3
    T6[3:, 3:] = R3

    # 18×18 block diagonal
    T = np.zeros((18, 18))
    for i in range(3):
        T[6*i:6*i+6, 6*i:6*i+6] = T6
    return T


# ══════════════════════════════════════════════════════════════════════════════
#  MITC3+ ELEMENT CLASS
# ══════════════════════════════════════════════════════════════════════════════

class MITC3Plus(Element):
    """
    MITC3+ Triangular Shell Element (Lee, Lee, Bathe 2014).

    Parameters
    ----------
    eid          : element ID
    nodes        : list of 3 Node objects
    E            : Young's modulus
    nu           : Poisson's ratio
    h            : shell thickness
    drilling_penalty : penalty stiffness for drilling DOF (default: auto)

    Notes
    -----
    - 6 DOF per node: [ux, uy, uz, rx, ry, rz]
    - Bubble DOF (α4, β4) condensed at element level
    - MITC3+ assumed transverse shear strain, d=1/10000
    - 7-point Gauss integration
    """

    ndof_per_node: int = 6

    def __init__(self, eid: int, nodes: List[Node],
                 E: float, nu: float, h: float,
                 drilling_penalty: float = None):
        if len(nodes) != 3:
            raise ValueError("MITC3+ requires exactly 3 nodes")
        self.eid   = eid
        self.nodes = nodes
        self.E     = float(E)
        self.nu    = float(nu)
        self.h     = float(h)

        # Auto drilling penalty: fraction of membrane stiffness
        if drilling_penalty is None:
            G = E / (2.0 * (1.0 + nu))
            self.drilling_penalty = 1e-5 * G * h   # calibrated: 1e-5 optimal for curved shells
        else:
            self.drilling_penalty = float(drilling_penalty)

        self._Ke_cache: np.ndarray = None

    def k_local(self) -> np.ndarray:
        """
        Compute 18×18 element stiffness in LOCAL frame (after static condensation).
        """
        if self._Ke_cache is not None:
            return self._Ke_cache

        K_full = _build_K_full(
            self.nodes, self.E, self.nu, self.h, self.drilling_penalty
        )
        Ke_loc = _static_condense(K_full)  # 18×18

        self._Ke_cache = Ke_loc
        return Ke_loc

    def T_matrix(self) -> np.ndarray:
        """18×18 transformation from local to global."""
        return _T_global(self.nodes)

    def k_global(self) -> np.ndarray:
        """18×18 element stiffness in GLOBAL frame."""
        T  = self.T_matrix()   # local→global: u_global = T @ u_local
        Kl = self.k_local()
        # K_global = T @ K_local @ T^T  (since u_local = T^T @ u_global)
        return T @ Kl @ T.T

    def global_dof_indices(self) -> list:
        dofs = []
        for node in self.nodes:
            dofs.extend(node.dofs[:6])
        return dofs

    def area(self) -> float:
        """Element area."""
        _, _, _, x_loc = _local_frame(self.nodes)
        xy = x_loc[:, :2]
        J  = _dh().T @ xy
        return 0.5 * abs(np.linalg.det(J))

    def aspect_ratio(self) -> float:
        """Max edge / min edge."""
        c = [n.coords for n in self.nodes]
        edges = [norm(c[1]-c[0]), norm(c[2]-c[1]), norm(c[0]-c[2])]
        return max(edges) / (min(edges) + 1e-30)

    def __repr__(self):
        nids = [n.nid for n in self.nodes]
        return f"MITC3Plus(eid={self.eid}, nodes={nids}, h={self.h})"

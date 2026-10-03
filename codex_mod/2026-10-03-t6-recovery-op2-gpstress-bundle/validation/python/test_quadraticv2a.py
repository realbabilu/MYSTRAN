"""
test_quadraticv2a.py  (revisi dari test_quadraticv2.py)
=======================
PERUBAHAN v2a:
  * T2 (patch bending) dan T5 (patch twist) untuk Q8: BC rotasi bisa tangan-kanan
    (thx=+dw/dy, thy=-dw/dx) per kelas lewat ROT_FLIP. Versi lama selalu tangan-kiri
    (thx=-dw/dy, thy=+dw/dx) sehingga elemen tangan-kanan (Simo v1p8/V10) tampak
    gagal 200% padahal elemennya benar. Kelas yang tidak tercantum di ROT_FLIP
    tetap memakai BC lama (perilaku tidak berubah). T6 (segitiga) tidak diubah.
  * Ditambah Simo_V10 (Simo1993_Q8_ShellElement_v10_standalone) dan Simo_V11 (..._v11_standalone) di Q8_ELEMS.
Test suite gabungan: Q8 Shell Elements + T6 Triangular Shell Element

Elemen Q8:
  1. HBQ8       :Darilmaz & Kumbasar 2006
  2. ANS8_BDG6  :Jung & Han 2013

Elemen T6:
  3. T6_Shell   :Simo1993 T6 v1.8 (FIXED)

Tests:
  T1.  Zero eigenvalue         (target = 6)
  T2.  Patch bending           (full-prescribed, tol < 1%)
  T3.  Patch membran           (full-prescribed, tol < 1%)
  T4.  Patch transverse shear  (1 elemen, tol < 5%)
  T5.  Patch twist             (1 elemen, tol < 5%)
  T6.  Kantilever tip momen    (w/w_ref konvergensi)
  T7.  Kantilever tip beban    (w/w_ref konvergensi)
  T8.  Pelat SSSS              (w_center/w_ref)
  T9.  Pinched cylinder        (w_A/w_ref, ref=1.8248e-5)
  T11. Distorted panel         (self-convergence)
"""

import sys, types, re
import numpy as np
sys.path.insert(0, '/mnt/user-data/uploads')
sys.path.insert(0, '/mnt/user-data/outputs')
from core import Model
from Simo1993_Tri6_ShellElement_v2 import Simo1993_Tri6_ShellElement_v2
from MacNeal_MH6T_Tri_v3 import MacNeal_MH6T_Tri_v3
from MITC6_Tri_v4 import MITC6_Tri_v4
from Rezaiee2017_Tri6_v3 import Rezaiee2017_Tri6_ShellElement_v3
from HBQ8_ShellElement      import HBQ8_ShellElement
from ANS8_BDG6_ShellElement import ANS8_BDG6_ShellElement
from Simo1993_Q8_ShellElement_v1p8_standalone import Simo1993_Q8_ShellElement_v1p8
from Simo1993_Q8_ShellElement_v10_standalone import Simo1993_Q8_ShellElement_v10
from Simo1993_Q8_ShellElement_v11_standalone import Simo1993_Q8_ShellElement_v11
from MITC8_ShellElement_v1p2 import MITC8_ShellElement_v1p2
from  MITC6_Tri_v1 import  MITC6_Tri_v1
from  MITC6_Tri_v0 import  MITC6_Tri_v0

# #########Daftar elemen #########
Q8_ELEMS = {'HBQ8': HBQ8_ShellElement, 'ANS8_BDG6': ANS8_BDG6_ShellElement,
            'MITC8_v1p2': MITC8_ShellElement_v1p2,
            'SimoQ8': Simo1993_Q8_ShellElement_v1p8,
            'Simo_V10': Simo1993_Q8_ShellElement_v10,
            'Simo_V11': Simo1993_Q8_ShellElement_v11}

# Konvensi BC rotasi untuk T2/T5 (Q8): +1 = BC lama (tangan-kiri), -1 = tangan-kanan.
# Hanya kelas yang DIKETAHUI tangan-kanan dicantumkan; sisanya default +1 (tidak berubah).
ROT_FLIP = {Simo1993_Q8_ShellElement_v1p8: -1, Simo1993_Q8_ShellElement_v10: -1, Simo1993_Q8_ShellElement_v11: -1}
def _rf(cls):
    return ROT_FLIP.get(cls, +1)
T6_ELEMS = {
    "SIMOT6": Simo1993_Tri6_ShellElement_v2,
    "MITC6": MITC6_Tri_v4,
    "MH6T": MacNeal_MH6T_Tri_v3,
    "REZAIEE": Rezaiee2017_Tri6_ShellElement_v3,
}

# #########Helpers #########
def Cmat(E, nu):
    f = E / (1 - nu ** 2)
    return f * np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]])

SEP = '\u2500' * 86
SEP2 = '\u2550' * 86

def mark(b):
    return '\u2713 PASS' if b else '\u2717 FAIL'

def hdr(s):
    print(f'\n{SEP2}\n  {s}\n{SEP2}')

def trow(name, res, Ns, w=10):
    v = {r[0]: r[-1] for r in res}
    ps = [f"{v.get(N, np.nan):>{w}.4f}" for N in Ns]
    fin = v.get(Ns[-1], np.nan)
    ok = '\u2713' if (not np.isnan(fin)) and 0.9 < fin < 1.1 else '\u2717'
    print(f"  {name:18s}  {'  '.join(ps)}  {ok}")

#  ######### Q8 mesh builder #########
def q8_rect(Lx, Ly, Nx, Ny, E, nu, h, cls, bc_fn=None, **kw):
    fem = Model(ndof_per_node=6, autospc=True, autospc_tol=1e-10)
    dx = Lx / Nx
    dy = Ly / Ny
    corner = {}
    nid = 1
    for j in range(Ny + 1):
        for i in range(Nx + 1):
            fem.add_node(nid, i * dx, j * dy, 0.)
            corner[(i, j)] = nid
            nid += 1
    hmid = {}
    for j in range(Ny + 1):
        for i in range(Nx):
            fem.add_node(nid, (i + .5) * dx, j * dy, 0.)
            hmid[(i, j)] = nid
            nid += 1
    vmid = {}
    for j in range(Ny):
        for i in range(Nx + 1):
            fem.add_node(nid, i * dx, (j + .5) * dy, 0.)
            vmid[(i, j)] = nid
            nid += 1
    elems = []
    eid = 1
    for j in range(Ny):
        for i in range(Nx):
            ns = [corner[(i, j)], corner[(i + 1, j)], corner[(i + 1, j + 1)],
                  corner[(i, j + 1)], hmid[(i, j)], vmid[(i + 1, j)],
                  hmid[(i, j + 1)], vmid[(i, j)]]
            e = cls(eid, [fem.nodes[n] for n in ns], E, nu, h, **kw)
            fem.add_element(e)
            elems.append(e)
            eid += 1
    if bc_fn:
        bc_fn(fem, corner, hmid, vmid, Nx, Ny)
    fem.solve_static()
    return fem, elems, corner, hmid, vmid

#  ######### T6 mesh builder #########
def t6_rect(Lx, Ly, Nx, Ny, E, nu, h, cls, bc_fn=None, **kw):
    fem = Model(ndof_per_node=6, autospc=True, autospc_tol=1e-10)
    dx = Lx / Nx
    dy = Ly / Ny
    corner = {}
    nid = 1
    for j in range(Ny + 1):
        for i in range(Nx + 1):
            fem.add_node(nid, i * dx, j * dy, 0.)
            corner[(i, j)] = nid
            nid += 1
    edge_mid = {}

    def get_mid(a, b):
        nonlocal nid
        key = tuple(sorted([a, b]))
        if key not in edge_mid:
            xa, ya, _ = fem.nodes[a].x, fem.nodes[a].y, 0.
            xb, yb, _ = fem.nodes[b].x, fem.nodes[b].y, 0.
            fem.add_node(nid, (xa + xb) / 2, (ya + yb) / 2, 0.)
            edge_mid[key] = nid
            nid += 1
        return edge_mid[key]

    elems = []
    eid = 1
    for j in range(Ny):
        for i in range(Nx):
            n1 = corner[(i, j)]
            n2 = corner[(i + 1, j)]
            n3 = corner[(i + 1, j + 1)]
            n4 = corner[(i, j + 1)]
            m12 = get_mid(n1, n2)
            m24 = get_mid(n2, n4)
            m41 = get_mid(n4, n1)
            e1 = cls(eid, [fem.nodes[n] for n in [n1, n2, n4, m12, m24, m41]],
                     E, nu, h, **kw)
            fem.add_element(e1)
            elems.append(e1)
            eid += 1
            m23 = get_mid(n2, n3)
            m34 = get_mid(n3, n4)
            e2 = cls(eid, [fem.nodes[n] for n in [n2, n3, n4, m23, m34, m24]],
                     E, nu, h, **kw)
            fem.add_element(e2)
            elems.append(e2)
            eid += 1

    if bc_fn:
        bc_fn(fem, corner, edge_mid, Nx, Ny)
    fem.solve_static()
    return fem, elems, corner, edge_mid

#  ######### Ekstrak Bm/Bb #########
def get_Bm(e):
    if hasattr(e, '_compute_Bm'):
        xi0 = 1. / 3.
        try:
            return e._compute_Bm(xi0, xi0)
        except Exception:
            return e._compute_Bm(0., 0.)
    return None

def get_Bb(e):
    if hasattr(e, '_compute_Bb'):
        xi0 = 1. / 3.
        try:
            return e._compute_Bb(xi0, xi0)
        except Exception:
            return e._compute_Bb(0., 0.)
    return None


# ----------------------------------------------------------------------
#  GLOBAL elem_stress (rotates local → global)
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
#  FULL 3D STRESS ROTATION HELPER
# ----------------------------------------------------------------------
def _get_element_basis(e):
    """Return (e1, e2) orthonormal basis for the element.
    Uses element.E1/.E2 if available; otherwise defaults to global X,Y
    (valid for all flat patch tests in the XY plane)."""
    try:
        return e.E1, e.E2
    except AttributeError:
        # Fallback for element classes that don't store local basis
        return np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])


def _rotate_stress_to_global(s_local, e1, e2):
    """
    Rotate a 3-component local stress vector [sxx, syy, sxy]
    to global coordinates using the element's local orthonormal
    basis vectors e1, e2 (e3 = e1 x e2).
    """
    e1 = np.asarray(e1, dtype=float)
    e2 = np.asarray(e2, dtype=float)
    e3 = np.cross(e1, e2)
    Q = np.column_stack([e1, e2, e3])           # 3x3 rotation matrix
    S_loc = np.array([[s_local[0], s_local[2], 0.0],
                      [s_local[2], s_local[1], 0.0],
                      [0.0,        0.0,        0.0]])
    S_glob = Q @ S_loc @ Q.T
    # Return global components in order [Sxx, Syy, Sxy]
    return np.array([S_glob[0, 0], S_glob[1, 1], S_glob[0, 1]])


# ----------------------------------------------------------------------
#  GLOBAL elem_stress (rotates local → global)
# ----------------------------------------------------------------------
def elem_stress(e, u):
    """Compute element stress in GLOBAL coordinates (rotated from local)."""
    idx = e.global_dof_indices()
    u_loc = e.T_matrix() @ u[idx]
    C = Cmat(e.E, e.nu)
    h_val = float(np.mean(e.h))

    Bm = get_Bm(e)
    Bb = get_Bb(e)

    sig_local = C @ (Bm @ u_loc) if Bm is not None else np.zeros(3)
    mom_local = (C @ (Bb @ u_loc)) * (h_val ** 3 / 12.0) if Bb is not None else np.zeros(3)

    e1, e2 = _get_element_basis(e)
    sig_global = _rotate_stress_to_global(sig_local, e1, e2)
    mom_global = _rotate_stress_to_global(mom_local, e1, e2)

    return sig_global, mom_global


#------------- before
def elem_stress_1(e, u):
    idx = e.global_dof_indices()
    u_loc = e.T_matrix() @ u[idx]
    C = Cmat(e.E, e.nu)
    h_val = float(np.mean(e.h))
    Bm = get_Bm(e)
    Bb = get_Bb(e)
    sig = C @ (Bm @ u_loc) if Bm is not None else np.zeros(3)
    mom = C @ (Bb @ u_loc) * h_val ** 3 / 12. if Bb is not None else np.zeros(3)
    return sig, mom
#------------- before old
def elem_stress2(e, u):
    idx = e.global_dof_indices()
    u_loc = e.T_matrix() @ u[idx]
    C = Cmat(e.E, e.nu)
    h = e.h
    Bm = get_Bm(e)
    Bb = get_Bb(e)
    sig = C @ (Bm @ u_loc) if Bm is not None else np.zeros(3)
    mom = C @ (Bb @ u_loc) * h ** 3 / 12. if Bb is not None else np.zeros(3)
    return sig, mom

# #########
# T1: ZERO EIGENVALUE
# #########
def t1_q8(cls, **kw):
    fem = Model(ndof_per_node=6, autospc=False)
    for i, (x, y, z) in enumerate([(-1, -1, 0), (1, -1, 0), (1, 1, 0), (-1, 1, 0),
                                     (0, -1, 0), (1, 0, 0), (0, 1, 0), (-1, 0, 0)], 1):
        fem.add_node(i, x, y, z)
    e = cls(1, [fem.nodes[i] for i in range(1, 9)], 2.1e5, 0.3, 0.01, **kw)
    K = e.k_local()
    v = np.linalg.eigvalsh(K)
    pos = v[v > 1.]
    thr = 1e-5 * pos.min() if len(pos) > 0 else 1e-6
    nz = int(np.sum(np.abs(v) < thr))
    nn = int(np.sum(v < -thr))
    return nz == 6 and nn == 0, nz, nn, v

def t1_t6(cls, **kw):
    fem = Model(ndof_per_node=6, autospc=False)
    for i, (x, y, z) in enumerate([(0, 0, 0), (1, 0, 0), (0, 1, 0),
                                     (0.5, 0, 0), (0.5, 0.5, 0), (0, 0.5, 0)], 1):
        fem.add_node(i, x, y, z)
    e = cls(1, [fem.nodes[i] for i in range(1, 7)], 2.1e5, 0.3, 0.01, **kw)
    K = e.k_local()
    v = np.linalg.eigvalsh(K)
    pos = v[v > 1e-3]
    thr = 1e-5 * pos.min() if len(pos) > 0 else 1e-8
    nz = int(np.sum(np.abs(v) < thr))
    nn = int(np.sum(v < -thr))
    return nz == 6 and nn == 0, nz, nn, v

# #########
# T2: PATCH BENDING (full-prescribed)
# #########
def _patch_bend_q8(cls, Ns, **kw):
    E = 1e6
    nu = 0.3
    h = 0.01
    a = b = 1e-4
    D = E * h ** 3 / (12 * (1 - nu ** 2))
    ref_mxx = D * (a + b * nu)
    ref_myy = D * (b + a * nu)
    res = {}
    for N, lbl in zip(Ns, ['1x1', '2x2', '3x3']):
        Lx = Ly = 0.12
        dx = Lx / N
        dy = Ly / N
        fem = Model(ndof_per_node=6, autospc=True, autospc_tol=1e-10)
        corner = {}
        nid = 1
        for j in range(N + 1):
            for i in range(N + 1):
                fem.add_node(nid, i * dx, j * dy, 0.)
                corner[(i, j)] = nid
                nid += 1
        hmid = {}
        for j in range(N + 1):
            for i in range(N):
                fem.add_node(nid, (i + .5) * dx, j * dy, 0.)
                hmid[(i, j)] = nid
                nid += 1
        vmid = {}
        for j in range(N):
            for i in range(N + 1):
                fem.add_node(nid, i * dx, (j + .5) * dy, 0.)
                vmid[(i, j)] = nid
                nid += 1
        elems = []
        eid = 1
        for j in range(N):
            for i in range(N):
                ns = [corner[(i, j)], corner[(i + 1, j)], corner[(i + 1, j + 1)],
                      corner[(i, j + 1)], hmid[(i, j)], vmid[(i + 1, j)],
                      hmid[(i, j + 1)], vmid[(i, j)]]
                e = cls(eid, [fem.nodes[n] for n in ns], E, nu, h, **kw)
                fem.add_element(e)
                elems.append(e)
                eid += 1
        for nid_, nd in fem.nodes.items():
            x, y = nd.x, nd.y
            fem.add_bc(nid_, 0, 0.)
            fem.add_bc(nid_, 1, 0.)
            fem.add_bc(nid_, 2, a * x ** 2 / 2 + b * y ** 2 / 2)
            fem.add_bc(nid_, 3, _rf(cls) * (-b * y))
            fem.add_bc(nid_, 4, _rf(cls) * (a * x))
            fem.add_bc(nid_, 5, 0.)
        fem.solve_static()
        errs = []
        for e in elems:
            _, mom = elem_stress(e, fem.u)
            errs.append(abs(mom[0] - ref_mxx) / (ref_mxx + 1e-30))
            errs.append(abs(mom[1] - ref_myy) / (ref_myy + 1e-30))
        res[lbl] = float(np.max(errs)) if errs else np.nan
    return res

def _patch_bend_t6(cls, Ns, **kw):
    E = 1e6
    nu = 0.3
    h = 0.01
    a = b = 1e-4
    D = E * h ** 3 / (12 * (1 - nu ** 2))
    ref_mxx = D * (a + b * nu)
    ref_myy = D * (b + a * nu)
    res = {}
    for N, lbl in zip(Ns, ['1x1', '2x2', '3x3']):
        Lx = Ly = 0.12
        dx = Lx / N
        dy = Ly / N
        fem = Model(ndof_per_node=6, autospc=True, autospc_tol=1e-10)
        corner = {}
        nid = 1
        for j in range(N + 1):
            for i in range(N + 1):
                fem.add_node(nid, i * dx, j * dy, 0.)
                corner[(i, j)] = nid
                nid += 1
        edge_mid = {}

        def get_mid(a_n, b_n):
            nonlocal nid
            key = tuple(sorted([a_n, b_n]))
            if key not in edge_mid:
                xa, ya = fem.nodes[a_n].x, fem.nodes[a_n].y
                xb, yb = fem.nodes[b_n].x, fem.nodes[b_n].y
                fem.add_node(nid, (xa + xb) / 2, (ya + yb) / 2, 0.)
                edge_mid[key] = nid
                nid += 1
            return edge_mid[key]

        elems = []
        eid = 1
        for j in range(N):
            for i in range(N):
                n1 = corner[(i, j)]
                n2 = corner[(i + 1, j)]
                n3 = corner[(i + 1, j + 1)]
                n4 = corner[(i, j + 1)]
                m12 = get_mid(n1, n2)
                m24 = get_mid(n2, n4)
                m41 = get_mid(n4, n1)
                m23 = get_mid(n2, n3)
                m34 = get_mid(n3, n4)
                e1 = cls(eid, [fem.nodes[n] for n in [n1, n2, n4, m12, m24, m41]],
                         E, nu, h, **kw)
                e2 = cls(eid + 1, [fem.nodes[n] for n in [n2, n3, n4, m23, m34, m24]],
                         E, nu, h, **kw)
                fem.add_element(e1)
                fem.add_element(e2)
                elems += [e1, e2]
                eid += 2
        for nid_, nd in fem.nodes.items():
            x, y = nd.x, nd.y
            fem.add_bc(nid_, 0, 0.)
            fem.add_bc(nid_, 1, 0.)
            fem.add_bc(nid_, 2, a * x ** 2 / 2 + b * y ** 2 / 2)
            fem.add_bc(nid_, 3, -b * y)
            fem.add_bc(nid_, 4, a * x)
            fem.add_bc(nid_, 5, 0.)
        fem.solve_static()
        errs = []
        for e in elems:
            _, mom = elem_stress(e, fem.u)
            errs.append(abs(mom[0] - ref_mxx) / (ref_mxx + 1e-30))
            errs.append(abs(mom[1] - ref_myy) / (ref_myy + 1e-30))
        res[lbl] = float(np.max(errs)) if errs else np.nan
    return res

# #########
# T3: PATCH MEMBRAN (full-prescribed)
# #########
def _patch_memb(cls, is_t6, Ns, **kw):
    E = 1e6
    nu = 0.3
    h = 0.01
    a_, b_, c_, d_ = 1e-5, 2e-5, 3e-5, 4e-5
    C = Cmat(E, nu)
    sig_ref = C @ np.array([a_, d_, b_ + c_])
    res = {}
    for N, lbl in zip(Ns, ['1x1', '2x2', '3x3']):
        Lx = Ly = 0.12
        dx = Lx / N
        dy = Ly / N
        fem = Model(ndof_per_node=6, autospc=True, autospc_tol=1e-10)
        corner = {}
        nid = 1
        for j in range(N + 1):
            for i in range(N + 1):
                fem.add_node(nid, i * dx, j * dy, 0.)
                corner[(i, j)] = nid
                nid += 1
        edge_mid = {} if is_t6 else None

        def get_mid_t6(a_n, b_n):
            nonlocal nid
            key = tuple(sorted([a_n, b_n]))
            if key not in edge_mid:
                xa, ya = fem.nodes[a_n].x, fem.nodes[a_n].y
                xb, yb = fem.nodes[b_n].x, fem.nodes[b_n].y
                fem.add_node(nid, (xa + xb) / 2, (ya + yb) / 2, 0.)
                edge_mid[key] = nid
                nid += 1
            return edge_mid[key]

        hmid = {}
        vmid = {}

        def get_mid_q8_h(i_, j_):
            nonlocal nid
            if (i_, j_) not in hmid:
                fem.add_node(nid, (i_ + .5) * dx, j_ * dy, 0.)
                hmid[(i_, j_)] = nid
                nid += 1
            return hmid[(i_, j_)]

        def get_mid_q8_v(i_, j_):
            nonlocal nid
            if (i_, j_) not in vmid:
                fem.add_node(nid, i_ * dx, (j_ + .5) * dy, 0.)
                vmid[(i_, j_)] = nid
                nid += 1
            return vmid[(i_, j_)]

        elems = []
        eid = 1
        for j in range(N):
            for i in range(N):
                n1 = corner[(i, j)]
                n2 = corner[(i + 1, j)]
                n3 = corner[(i + 1, j + 1)]
                n4 = corner[(i, j + 1)]
                if is_t6:
                    m12 = get_mid_t6(n1, n2)
                    m24 = get_mid_t6(n2, n4)
                    m41 = get_mid_t6(n4, n1)
                    m23 = get_mid_t6(n2, n3)
                    m34 = get_mid_t6(n3, n4)
                    e1 = cls(eid, [fem.nodes[n] for n in [n1, n2, n4, m12, m24, m41]],
                             E, nu, h, **kw)
                    e2 = cls(eid + 1, [fem.nodes[n] for n in [n2, n3, n4, m23, m34, m24]],
                             E, nu, h, **kw)
                    fem.add_element(e1)
                    fem.add_element(e2)
                    elems += [e1, e2]
                    eid += 2
                else:
                    mh = get_mid_q8_h(i, j)
                    mv = get_mid_q8_v(i + 1, j)
                    mht = get_mid_q8_h(i, j + 1)
                    mvl = get_mid_q8_v(i, j)
                    ns = [n1, n2, n3, n4, mh, mv, mht, mvl]
                    e = cls(eid, [fem.nodes[n] for n in ns], E, nu, h, **kw)
                    fem.add_element(e)
                    elems.append(e)
                    eid += 1
        for nid_, nd in fem.nodes.items():
            x, y = nd.x, nd.y
            fem.add_bc(nid_, 0, a_ * x + b_ * y)
            fem.add_bc(nid_, 1, c_ * x + d_ * y)
            fem.add_bc(nid_, 2, 0.)
            fem.add_bc(nid_, 3, 0.)
            fem.add_bc(nid_, 4, 0.)
            fem.add_bc(nid_, 5, 0.)
        fem.solve_static()
        errs = []
        for e in elems:
            sig, _ = elem_stress(e, fem.u)
            for k in range(3):
                errs.append(abs(sig[k] - sig_ref[k]) / (abs(sig_ref[k]) + 1e-30))
        res[lbl] = float(np.max(errs)) if errs else np.nan
    return res

# #########
# T4: PATCH TRANSVERSE SHEAR (1 elemen)
# #########
def t4_shear_q8(cls, **kw):
    E = 1e6
    nu = 0.3
    h = 0.1
    a = 1e-4
    G = E / (2 * (1 + nu))
    ks = 5. / 6.
    Qx_ref = ks * G * h * a
    fem = Model(ndof_per_node=6, autospc=True, autospc_tol=1e-10)
    for i, (x, y, z) in enumerate([(-1, -1, 0), (1, -1, 0), (1, 1, 0), (-1, 1, 0),
                                     (0, -1, 0), (1, 0, 0), (0, 1, 0), (-1, 0, 0)], 1):
        fem.add_node(i, x, y, z)
    e = cls(1, [fem.nodes[i] for i in range(1, 9)], E, nu, h, **kw)
    fem.add_element(e)
    for nid in range(1, 9):
        nd = fem.nodes[nid]
        x = nd.x
        fem.add_bc(nid, 0, 0.)
        fem.add_bc(nid, 1, 0.)
        fem.add_bc(nid, 2, a * x)
        fem.add_bc(nid, 3, 0.)
        fem.add_bc(nid, 4, 0.)
        fem.add_bc(nid, 5, 0.)
    fem.solve_static()
    idx = e.global_dof_indices()
    u_loc = e.T_matrix() @ fem.u[idx]
    if hasattr(e, '_compute_Bs'):
        Bs = e._compute_Bs(0., 0.)
    else:
        N8, dN = e._shape_q8(0., 0.) if hasattr(e, '_shape_q8') else (None, None)
        if N8 is None:
            return False, np.nan
        J = e._jacobian(0., 0.)
        Jinv = np.linalg.inv(J)
        dNxy = Jinv @ dN
        Bs = np.zeros((2, 48))
        for i in range(8):
            col = 6 * i
            Bs[0, col + 2] = dNxy[0, i]
            Bs[0, col + 4] = N8[i]
            Bs[1, col + 2] = dNxy[1, i]
            Bs[1, col + 3] = -N8[i]
    g = Bs @ u_loc
    Q = ks * G * np.eye(2) @ g * h
    err = abs(Q[0] - Qx_ref) / (Qx_ref + 1e-30)
    return err < 0.05, err

def t4_shear_t6(cls, **kw):
    E = 1e6
    nu = 0.3
    h = 0.1
    a = 1e-4
    G = E / (2 * (1 + nu))
    ks = 5. / 6.
    Qx_ref = ks * G * h * a
    fem = Model(ndof_per_node=6, autospc=True, autospc_tol=1e-10)
    for i, (x, y, z) in enumerate([(0, 0, 0), (1, 0, 0), (0, 1, 0),
                                     (0.5, 0, 0), (0.5, 0.5, 0), (0, 0.5, 0)], 1):
        fem.add_node(i, x, y, z)
    e = cls(1, [fem.nodes[i] for i in range(1, 7)], E, nu, h, **kw)
    fem.add_element(e)
    for nid in range(1, 7):
        nd = fem.nodes[nid]
        x = nd.x
        fem.add_bc(nid, 0, 0.)
        fem.add_bc(nid, 1, 0.)
        fem.add_bc(nid, 2, a * x)
        fem.add_bc(nid, 3, 0.)
        fem.add_bc(nid, 4, 0.)
        fem.add_bc(nid, 5, 0.)
    fem.solve_static()
    idx = e.global_dof_indices()
    u_loc = e.T_matrix() @ fem.u[idx]
    Bs = e._compute_Bs(1. / 3., 1. / 3.)
    g = Bs @ u_loc
    Q = ks * G * np.eye(2) @ g * h
    err = abs(Q[0] - Qx_ref) / (Qx_ref + 1e-30)
    return err < 0.05, err

# #########
# T5: PATCH TWIST (1 elemen)
# #########
def t5_twist(cls, is_t6, **kw):
    E = 1e6
    nu = 0.3
    h = 0.01
    a = 1e-4
    D = E * h ** 3 / (12 * (1 - nu ** 2))
    ref_mxy = D * (1 - nu) * a
    fem = Model(ndof_per_node=6, autospc=True, autospc_tol=1e-10)
    if is_t6:
        for i, (x, y, z) in enumerate([(0, 0, 0), (1, 0, 0), (0, 1, 0),
                                         (0.5, 0, 0), (0.5, 0.5, 0), (0, 0.5, 0)], 1):
            fem.add_node(i, x, y, z)
        e = cls(1, [fem.nodes[i] for i in range(1, 7)], E, nu, h, **kw)
        fem.add_element(e)
        for nid in range(1, 7):
            nd = fem.nodes[nid]
            x, y = nd.x, nd.y
            fem.add_bc(nid, 0, 0.)
            fem.add_bc(nid, 1, 0.)
            fem.add_bc(nid, 2, a * x * y)
            fem.add_bc(nid, 3, -a * x)
            fem.add_bc(nid, 4, a * y)
            fem.add_bc(nid, 5, 0.)
    else:
        for i, (x, y, z) in enumerate([(-1, -1, 0), (1, -1, 0), (1, 1, 0), (-1, 1, 0),
                                         (0, -1, 0), (1, 0, 0), (0, 1, 0), (-1, 0, 0)], 1):
            fem.add_node(i, x, y, z)
        e = cls(1, [fem.nodes[i] for i in range(1, 9)], E, nu, h, **kw)
        fem.add_element(e)
        for nid in range(1, 9):
            nd = fem.nodes[nid]
            x, y = nd.x, nd.y
            fem.add_bc(nid, 0, 0.)
            fem.add_bc(nid, 1, 0.)
            fem.add_bc(nid, 2, a * x * y)
            fem.add_bc(nid, 3, _rf(cls) * (-a * x))
            fem.add_bc(nid, 4, _rf(cls) * (a * y))
            fem.add_bc(nid, 5, 0.)
    fem.solve_static()
    idx = e.global_dof_indices()
    u_loc = e.T_matrix() @ fem.u[idx]
    Bb = get_Bb(e)
    C = Cmat(E, nu)
    mom = C @ (Bb @ u_loc) * h ** 3 / 12.
    err = abs(mom[2] - ref_mxy) / (abs(ref_mxy) + 1e-30)
    return err < 0.05, err

# #########
# T6-T8: KONVERGENSI (shared logic)
# #########
def t6_moment(cls, is_t6, Ns, **kw):
    L = 10.0
    b = 1.0
    h = 0.1
    E = 1e6
    nu = 0.3
    M = 1.0                     # total applied moment
    I = b * h**3 / 12.0
    w_ref = M * L**2 / (2.0 * E * I)

    results = []
    for N in Ns:
        Nx = N
        Ny = max(1, N // 5)
        try:
            if is_t6:
                # T6: trapezoidal lumping – correct total moment, but not
                # exactly consistent (no midside nodes on edge are loaded).
                # Still gives ~1.0 as seen in your output.
                def bc(fem, corner, emid, Nx_, Ny_):
                    for j in range(Ny_ + 1):
                        fem.fix_node(corner[(0, j)])
                    dy_ = b / Ny_
                    for j in range(Ny_):
                        fem.add_load(corner[(Nx_, j)],     4, -M * dy_ / 2.0)
                        fem.add_load(corner[(Nx_, j + 1)], 4, -M * dy_ / 2.0)
                fem, _, corner, _ = t6_rect(L, b, Nx, Ny, E, nu, h, cls, bc_fn=bc, **kw)
            else:
                # Q8: consistent quadratic loads for a constant distributed moment
                def bc(fem, corner, hmid, vmid, Nx_, Ny_):
                    for j in range(Ny_ + 1):
                        fem.fix_node(corner[(0, j)])
                    dy_ = b / Ny_
                    for j in range(Ny_):
                        fem.add_load(corner[(Nx_, j)],     4, -M * dy_ / 6.0)
                        fem.add_load(vmid[(Nx_, j)],       4, -M * dy_ * 2.0 / 3.0)
                        fem.add_load(corner[(Nx_, j + 1)], 4, -M * dy_ / 6.0)
                fem, _, corner, _, _ = q8_rect(L, b, Nx, Ny, E, nu, h, cls, bc_fn=bc, **kw)

            w = abs(fem.u[fem.nodes[corner[(Nx, 0)]].dofs[2]])
            results.append((N, w, w / w_ref))
        except Exception:
            results.append((N, np.nan, np.nan))

    fin = results[-1][2]
    return (not np.isnan(fin)) and 0.9 < fin < 1.1, results


def t7_load(cls, is_t6, Ns, **kw):
    L = 10.0
    b = 1.0
    h = 0.1
    E = 1e6
    nu = 0.0
    P = 1.0                     # total tip load
    I = b * h**3 / 12.0
    w_ref = P * L**3 / (3.0 * E * I)

    results = []
    for N in Ns:
        Nx = N
        Ny = max(1, N // 10)
        try:
            if is_t6:
                def bc(fem, corner, emid, Nx_, Ny_):
                    for j in range(Ny_ + 1):
                        fem.fix_node(corner[(0, j)])
                    for j in range(Ny_ + 1):
                        fem.add_load(corner[(Nx_, j)], 2, P / (Ny_ + 1))
                fem, _, corner, _ = t6_rect(L, b, Nx, Ny, E, nu, h, cls, bc_fn=bc, **kw)
            else:
                # Q8: equal division to corner nodes only (simplified).
                # For exact consistency you could use quadratic distribution,
                # but this works well enough for convergence.
                def bc(fem, corner, hmid, vmid, Nx_, Ny_):
                    for j in range(Ny_ + 1):
                        fem.fix_node(corner[(0, j)])
                    for j in range(Ny_ + 1):
                        fem.add_load(corner[(Nx_, j)], 2, P / (Ny_ + 1))
                fem, _, corner, _, _ = q8_rect(L, b, Nx, Ny, E, nu, h, cls, bc_fn=bc, **kw)

            w = abs(fem.u[fem.nodes[corner[(Nx, 0)]].dofs[2]])
            results.append((N, w, w / w_ref))
        except Exception:
            results.append((N, np.nan, np.nan))

    fin = results[-1][2]
    return (not np.isnan(fin)) and 0.9 < fin < 1.1, results

def t6_moment_old(cls, is_t6, Ns, **kw):
    L = 10.
    b = 1.
    h = 0.1
    E = 1e6
    nu = 0.3
    M = 1.
    I = b * h ** 3 / 12.
    w_ref = M * L ** 2 / (2. * E * I)
    results = []
    for N in Ns:
        Nx = N
        Ny = max(1, N // 5)
        try:
            if is_t6:
                def bc(fem, corner, emid, Nx_, Ny_):
                    for j in range(Ny_ + 1):
                        fem.fix_node(corner[(0, j)])
                    dy_ = b / Ny_
                    for j in range(Ny_ + 1):
                        wm = dy_ / 2 if (j == 0 or j == Ny_) else dy_
                        fem.add_load(corner[(Nx_, j)], 4, -M * wm)

                fem, _, corner, emid = t6_rect(L, b, Nx, Ny, E, nu, h, cls,
                                                  bc_fn=bc, **kw)
            else:
                def bc(fem, corner, hmid, vmid, Nx_, Ny_):
                    for j in range(Ny_ + 1):
                        fem.fix_node(corner[(0, j)])
                    dy_ = b / Ny_
                    for j in range(Ny_ + 1):
                        wm = dy_ / 2 if (j == 0 or j == Ny_) else dy_
                        fem.add_load(corner[(Nx_, j)], 4, -M * wm)
                    for j in range(Ny_):
                        fem.add_load(vmid[(Nx_, j)], 4, -M * dy_ * (2. / 3.))
                    for j in range(Ny_ + 1):
                        wm = dy_ / 6 if (j == 0 or j == Ny_) else dy_ / 3.
                        fem.add_load(corner[(Nx_, j)], 4, -M * wm)

                fem, _, corner, hmid, vmid = q8_rect(L, b, Nx, Ny, E, nu, h, cls,
                                                      bc_fn=bc, **kw)
            w = abs(fem.u[fem.nodes[corner[(Nx, 0)]].dofs[2]])
            results.append((N, w, w / w_ref))
        except Exception:
            results.append((N, np.nan, np.nan))
    fin = results[-1][2]
    return (not np.isnan(fin)) and 0.9 < fin < 1.1, results

def t7_load_old(cls, is_t6, Ns, **kw):
    L = 10.
    b = 1.
    h = 0.1
    E = 1e6
    nu = 0.0
    P = 1.
    I = b * h ** 3 / 12.
    w_ref = P * L ** 3 / (3. * E * I)
    results = []
    for N in Ns:
        Nx = N
        Ny = max(1, N // 10)
        try:
            if is_t6:
                def bc(fem, corner, emid, Nx_, Ny_):
                    for j in range(Ny_ + 1):
                        fem.fix_node(corner[(0, j)])
                    for j in range(Ny_ + 1):
                        fem.add_load(corner[(Nx_, j)], 2, P / (Ny_ + 1))

                fem, _, corner, emid = t6_rect(L, b, Nx, Ny, E, nu, h, cls,
                                                  bc_fn=bc, **kw)
            else:
                def bc(fem, corner, hmid, vmid, Nx_, Ny_):
                    for j in range(Ny_ + 1):
                        fem.fix_node(corner[(0, j)])
                    for j in range(Ny_ + 1):
                        fem.add_load(corner[(Nx_, j)], 2, P / (Ny_ + 1))

                fem, _, corner, hmid, vmid = q8_rect(L, b, Nx, Ny, E, nu, h, cls,
                                                      bc_fn=bc, **kw)
            w = abs(fem.u[fem.nodes[corner[(Nx, 0)]].dofs[2]])
            results.append((N, w, w / w_ref))
        except Exception:
            results.append((N, np.nan, np.nan))
    fin = results[-1][2]
    return (not np.isnan(fin)) and 0.9 < fin < 1.1, results

def t8_ssplate(cls, is_t6, Ns, **kw):
    a_ = 1.
    q = 1.
    E = 1e6
    nu = 0.3
    h = 0.01
    D = E * h ** 3 / (12 * (1 - nu ** 2))
    w_ref = 0.00406 * q * a_ ** 4 / D
    results = []
    for N in Ns:
        dx = a_ / N
        dy = a_ / N
        try:
            if is_t6:
                def bc(fem, corner, emid, Nx_, Ny_):
                    for j in range(Ny_ + 1):
                        fem.add_bc(corner[(0, j)], 2, 0.)
                        fem.add_bc(corner[(Nx_, j)], 2, 0.)
                    for i in range(Nx_ + 1):
                        fem.add_bc(corner[(i, 0)], 2, 0.)
                        fem.add_bc(corner[(i, Ny_)], 2, 0.)
                    for j in range(Ny_ + 1):
                        for i in range(Nx_ + 1):
                            fx = 0.5 if (i == 0 or i == Nx_) else 1.
                            fy = 0.5 if (j == 0 or j == Ny_) else 1.
                            fem.add_load(corner[(i, j)], 2, q * dx * dy * fx * fy)

                fem, _, corner, emid = t6_rect(a_, a_, N, N, E, nu, h, cls,
                                                bc_fn=bc, **kw)
            else:
                def bc(fem, corner, hmid, vmid, Nx_, Ny_):
                    for j in range(Ny_ + 1):
                        fem.add_bc(corner[(0, j)], 2, 0.)
                        fem.add_bc(corner[(Nx_, j)], 2, 0.)
                    for i in range(Nx_ + 1):
                        fem.add_bc(corner[(i, 0)], 2, 0.)
                        fem.add_bc(corner[(i, Ny_)], 2, 0.)
                    for j in range(Ny_ + 1):
                        for i in range(Nx_ + 1):
                            fx = 0.5 if (i == 0 or i == Nx_) else 1.
                            fy = 0.5 if (j == 0 or j == Ny_) else 1.
                            fem.add_load(corner[(i, j)], 2, q * dx * dy * fx * fy)

                fem, _, corner, hmid, vmid = q8_rect(a_, a_, N, N, E, nu, h, cls,
                                                      bc_fn=bc, **kw)
            mid = N // 2
            w = abs(fem.u[fem.nodes[corner[(mid, mid)]].dofs[2]])
            results.append((N, w, w / w_ref))
        except Exception:
            results.append((N, np.nan, np.nan))
    fin = results[-1][2]
    return (not np.isnan(fin)) and 0.9 < fin < 1.1, results

# #########
# T9: PINCHED CYLINDER (shared, pakai Q8 atau T6)
# #########
def t9_cyl(cls, is_t6, Ns, **kw):
    R = 300.
    L = 600.
    h = 3.
    E = 3e6
    nu = 0.3
    P = 1.
    w_ref = 1.8248e-5
    results = []
    for N in Ns:
        Nth = N
        Nz = N
        dth = (np.pi / 2) / Nth
        dz = (L / 2) / Nz
        try:
            fem = Model(ndof_per_node=6, autospc=True, autospc_tol=1e-10)
            nid = 1
            corner = {}
            hmid = {}
            vmid = {}
            for j in range(Nz + 1):
                for i in range(Nth + 1):
                    th = i * dth
                    z = j * dz
                    fem.add_node(nid, R * np.cos(th), R * np.sin(th), z)
                    corner[(i, j)] = nid
                    nid += 1
            if not is_t6:
                for j in range(Nz + 1):
                    for i in range(Nth):
                        th = (i + .5) * dth
                        z = j * dz
                        fem.add_node(nid, R * np.cos(th), R * np.sin(th), z)
                        hmid[(i, j)] = nid
                        nid += 1
                for j in range(Nz):
                    for i in range(Nth + 1):
                        th = i * dth
                        z = (j + .5) * dz
                        fem.add_node(nid, R * np.cos(th), R * np.sin(th), z)
                        vmid[(i, j)] = nid
                        nid += 1
            else:
                edge_mid = {}

                def get_mid_c(a_n, b_n):
                    nonlocal nid
                    key = tuple(sorted([a_n, b_n]))
                    if key not in edge_mid:
                        xa, ya, za = fem.nodes[a_n].x, fem.nodes[a_n].y, fem.nodes[a_n].z
                        xb, yb, zb = fem.nodes[b_n].x, fem.nodes[b_n].y, fem.nodes[b_n].z
                        xm, ym, zm = (xa + xb) / 2, (ya + yb) / 2, (za + zb) / 2
                        r = np.sqrt(xm ** 2 + ym ** 2)
                        fac = R / r if r > 1e-10 else 1.
                        fem.add_node(nid, xm * fac, ym * fac, zm)
                        edge_mid[key] = nid
                        nid += 1
                    return edge_mid[key]

            eid = 1
            for j in range(Nz):
                for i in range(Nth):
                    n1 = corner[(i, j)]
                    n2 = corner[(i + 1, j)]
                    n3 = corner[(i + 1, j + 1)]
                    n4 = corner[(i, j + 1)]
                    if is_t6:
                        m12 = get_mid_c(n1, n2)
                        m24 = get_mid_c(n2, n4)
                        m41 = get_mid_c(n4, n1)
                        m23 = get_mid_c(n2, n3)
                        m34 = get_mid_c(n3, n4)
                        e1 = cls(eid, [fem.nodes[n] for n in [n1, n2, n4, m12, m24, m41]],
                                 E, nu, h, **kw)
                        e2 = cls(eid + 1, [fem.nodes[n] for n in [n2, n3, n4, m23, m34, m24]],
                                 E, nu, h, **kw)
                        fem.add_element(e1)
                        fem.add_element(e2)
                        eid += 2
                    else:
                        ns = [n1, n2, n3, n4, hmid[(i, j)], vmid[(i + 1, j)],
                              hmid[(i, j + 1)], vmid[(i, j)]]
                        fem.add_element(cls(eid, [fem.nodes[n] for n in ns],
                                            E, nu, h, **kw))
                        eid += 1
            for j in range(Nz + 1):
                fem.add_bc(corner[(0, j)], 1, 0.)
                fem.add_bc(corner[(0, j)], 5, 0.)
                fem.add_bc(corner[(Nth, j)], 0, 0.)
                fem.add_bc(corner[(Nth, j)], 5, 0.)
            for i in range(Nth + 1):
                fem.add_bc(corner[(i, 0)], 2, 0.)
                fem.add_bc(corner[(i, 0)], 3, 0.)
                fem.add_bc(corner[(i, 0)], 4, 0.)
                fem.add_bc(corner[(i, Nz)], 2, 0.)
                th = i * dth
                if np.cos(th) > 0.7:
                    fem.add_bc(corner[(i, Nz)], 1, 0.)
                if np.sin(th) > 0.7:
                    fem.add_bc(corner[(i, Nz)], 0, 0.)
            fem.add_load(corner[(0, 0)], 0, -P / 4.)
            fem.solve_static()
            w = abs(fem.u[fem.nodes[corner[(0, 0)]].dofs[0]])
            results.append((N, w, w / w_ref))
        except Exception as ex:
            results.append((N, np.nan, np.nan))
    fin = results[-1][2]
    return (not np.isnan(fin)) and 0.8 < fin < 1.2, results

# #########
# T11: DISTORTED PANEL
# #########
def t11_panel(cls, is_t6, Ns, **kw):
    L = 10.
    b = 1.
    h = 0.1
    E = 1e6
    nu = 0.3
    P = 1.
    res_abs = []
    for N in Ns:
        Nx = N
        Ny = max(1, N // 5)
        try:
            if is_t6:
                def bc(fem, corner, emid, Nx_, Ny_):
                    for j in range(Ny_ + 1):
                        fem.fix_node(corner[(0, j)])
                    for j in range(Ny_ + 1):
                        fem.add_load(corner[(Nx_, j)], 1, P / (Ny_ + 1))

                fem, _, corner, emid = t6_rect(L, b, Nx, Ny, E, nu, h, cls,
                                                  bc_fn=bc, **kw)
            else:
                def bc(fem, corner, hmid, vmid, Nx_, Ny_):
                    for j in range(Ny_ + 1):
                        fem.fix_node(corner[(0, j)])
                    for j in range(Ny_ + 1):
                        fem.add_load(corner[(Nx_, j)], 1, P / (Ny_ + 1))

                fem, _, corner, hmid, vmid = q8_rect(L, b, Nx, Ny, E, nu, h, cls,
                                                      bc_fn=bc, **kw)
            u = fem.u[fem.nodes[corner[(Nx, 0)]].dofs[1]]
            res_abs.append((N, u))
        except Exception:
            res_abs.append((N, np.nan))
    ref = res_abs[-1][1]
    return [(N, v, v / ref if not np.isnan(ref) and ref != 0 else np.nan)
            for N, v in res_abs]

# #########
# MAIN
# #########
if __name__ == '__main__':
    print(SEP2)
    print('  TEST SUITE GABUNGAN: Q8 + T6 SHELL ELEMENTS')
    print(SEP2)

    ALL = {**{n: (c, False) for n, c in Q8_ELEMS.items()},
           **{n: (c, True) for n, c in T6_ELEMS.items()}}

    # ######### T1 #########
    hdr('T1: ZERO EIGENVALUE  (target = 6 zero modes, 0 negative)')
    print(f"  {'Elemen':18s}  {'nzero':>6}  {'nneg':>5}  {'Status':>8}  EV[6:9]")
    print(f"  {'-'*18}  {'-'*6}  {'-'*5}  {'-'*8}  {'-'*30}")
    for name, (cls, is_t6) in ALL.items():
        fn = t1_t6 if is_t6 else t1_q8
        ok, nz, nn, v = fn(cls)
        ev = [f'{x:.2e}' for x in v[6:9]] if isinstance(v, np.ndarray) else [str(v)]
        print(f"  {name:18s}  {nz:>6}  {nn:>5}  {mark(ok):>8}  {ev}")

    # ######### T2 #########
    Ns_p = [1, 2, 3]
    lbls = ['1x1', '2x2', '3x3']
    hdr('T2: PATCH BENDING  (full-prescribed semua node, tol < 1%)')
    print(f"  {'Elemen':18s}  {'1x1':>9}  {'2x2':>9}  {'3x3':>9}  Status")
    print(f"  {'-'*18}  {'-'*9}  {'-'*9}  {'-'*9}  ------")
    for name, (cls, is_t6) in ALL.items():
        try:
            res = _patch_bend_t6(cls, Ns_p) if is_t6 else _patch_bend_q8(cls, Ns_p)
            vals = [res.get(l, np.nan) for l in lbls]
            ok = all(v < 0.01 for v in vals if not np.isnan(v))
            cells = [f"{v*100:8.3f}%" if not np.isnan(v) else "    N/A  " for v in vals]
            print(f"  {name:18s}  {'  '.join(cells)}  {mark(ok)}")
        except Exception as ex:
            print(f"  {name:18s}  ERROR: {ex}")

    # ######### T3 #########
    hdr('T3: PATCH MEMBRAN  (full-prescribed semua node, tol < 1%)')
    print(f"  {'Elemen':18s}  {'1x1':>9}  {'2x2':>9}  {'3x3':>9}  Status")
    print(f"  {'-'*18}  {'-'*9}  {'-'*9}  {'-'*9}  ------")
    for name, (cls, is_t6) in ALL.items():
        try:
            res = _patch_memb(cls, is_t6, Ns_p)
            vals = [res.get(l, np.nan) for l in lbls]
            ok = all(v < 0.01 for v in vals if not np.isnan(v))
            cells = [f"{v*100:8.3f}%" if not np.isnan(v) else "    N/A  " for v in vals]
            print(f"  {name:18s}  {'  '.join(cells)}  {mark(ok)}")
        except Exception as ex:
            print(f"  {name:18s}  ERROR: {ex}")

    # ######### T4 #########
    hdr('T4: PATCH TRANSVERSE SHEAR  (1 elemen, w=ax theta_y=0, tol < 5%)')
    print(f"  {'Elemen':18s}  {'Error':>10}  Status")
    print(f"  {'-'*18}  {'-'*10}  ------")
    for name, (cls, is_t6) in ALL.items():
        fn = t4_shear_t6 if is_t6 else t4_shear_q8
        ok, err = fn(cls)
        cell = f"{err*100:9.3f}%" if not np.isnan(err) else "     N/A "
        print(f"  {name:18s}  {cell}  {mark(ok)}")

    # ######### T5 #########
    hdr('T5: PATCH TWIST  (1 elemen, kappa_xy konstan, tol < 5%)')
    print(f"  {'Elemen':18s}  {'Error':>10}  Status")
    print(f"  {'-'*18}  {'-'*10}  ------")
    for name, (cls, is_t6) in ALL.items():
        ok, err = t5_twist(cls, is_t6)
        cell = f"{err*100:9.3f}%" if not np.isnan(err) else "     N/A "
        print(f"  {name:18s}  {cell}  {mark(ok)}")

    # ######### T6 #########
    Ns = [2, 4, 8, 16]
    hdr('T6: KANTILEVER TIP MOMEN  (|w|/w_ref, ref=ML^2/2EI)')
    print(f"  {'Elemen':18s}  {'N=2':>10}  {'N=4':>10}  {'N=8':>10}  {'N=16':>10}  OK?")
    print(f"  {'-'*18}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  ---")
    for name, (cls, is_t6) in ALL.items():
        _, res = t6_moment(cls, is_t6, Ns)
        trow(name, res, Ns)

    # ######### T7 #########
    hdr('T7: KANTILEVER TIP BEBAN  (|w|/w_ref, ref=PL^3/3EI, nu=0)')
    print(f"  {'Elemen':18s}  {'N=2':>10}  {'N=4':>10}  {'N=8':>10}  {'N=16':>10}  OK?")
    print(f"  {'-'*18}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  ---")
    for name, (cls, is_t6) in ALL.items():
        _, res = t7_load(cls, is_t6, Ns)
        trow(name, res, Ns)

    # ######### T8 #########
    hdr('T8: PELAT SSSS  (w_center/w_ref, ref=0.00406qa^4/D)')
    print(f"  {'Elemen':18s}  {'N=2':>10}  {'N=4':>10}  {'N=8':>10}  {'N=16':>10}  OK?")
    print(f"  {'-'*18}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  ---")
    for name, (cls, is_t6) in ALL.items():
        _, res = t8_ssplate(cls, is_t6, Ns)
        trow(name, res, Ns)

    # ######### T9 #########
    Ns9 = [4, 8, 16]
    hdr('T9: PINCHED CYLINDER  (w/w_ref, ref=1.8248e-5, Belytschko 1985)')
    print(f"  {'Elemen':18s}  {'N=4':>10}  {'N=8':>10}  {'N=16':>10}  OK?")
    print(f"  {'-'*18}  {'-'*10}  {'-'*10}  {'-'*10}  ---")
    for name, (cls, is_t6) in ALL.items():
        _, res = t9_cyl(cls, is_t6, Ns9)
        trow(name, res, Ns9)

    # ######### T11 #########
    hdr('T11: DISTORTED PANEL  (self-convergence, norm ke N=8)')
    print(f"  {'Elemen':18s}  {'N=2':>10}  {'N=4':>10}  {'N=8':>10}")
    print(f"  {'-'*18}  {'-'*10}  {'-'*10}  {'-'*10}")
    for name, (cls, is_t6) in ALL.items():
        res = t11_panel(cls, is_t6, [2, 4, 8])
        v = {r[0]: r[2] for r in res}
        ps = [f"{v.get(N, np.nan):>10.4f}" for N in [2, 4, 8]]
        print(f"  {name:18s}  {'  '.join(ps)}")

    # ######### RINGKASAN #########
    print(f'\n{SEP2}')
    print('  INTERPRETASI:')
    print('  T1     nzero=6, nneg=0      -> tidak ada spurious mode / mode negatif')
    print('  T2/T3  error < 1%           -> patch test lulus (konsistensi variasional)')
    print('  T4/T5  error < 5%           -> shear/twist reproducible dengan benar')
    print('  T6-T8  ratio 0.90-1.10      -> konvergensi baik di mesh sedang')
    print('  T9     ratio N=16 ~ 1.0     -> lulus obstacle course (shell tipis)')
    print('  T11    monoton -> 1.0        -> konvergensi stabil pada geometri distorted')
    print(f'  Q8 DOF/elem=48, T6 DOF/elem=36; T6 ~25% lebih hemat per elemen')
    print(SEP2)

def _rotate_stress_to_global(s_local, e1, e2):
    """
    Rotate a 3-component local stress vector [sxx, syy, sxy]
    to global coordinates, using the element's local orthonormal basis
    vectors e1, e2 (and e3 = e1 x e2).
    This is identical to core._rotate_inplane_tensor_to_global.
    """
    e1 = np.asarray(e1, dtype=float)
    e2 = np.asarray(e2, dtype=float)
    e3 = np.cross(e1, e2)
    Q = np.column_stack([e1, e2, e3])           # 3x3 rotation matrix
    # Build local stress tensor in the (e1,e2) plane
    S_loc = np.array([[s_local[0], s_local[2], 0.0],
                      [s_local[2], s_local[1], 0.0],
                      [0.0,        0.0,        0.0]])
    S_glob = Q @ S_loc @ Q.T
    # Extract global components in standard order [sxx, syy, sxy]
    return np.array([S_glob[0, 0], S_glob[1, 1], S_glob[0, 1]])



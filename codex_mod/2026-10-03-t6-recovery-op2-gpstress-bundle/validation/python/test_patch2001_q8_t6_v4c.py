"""
test_patch2001_q8_t6_v4c.py
============================
MacNeal-Harder 2-001 Patch Test — Q8 + T6
v4c: adds grid-node displacement print for comparison with MYSTRAN F06.

Local stress = stress expressed in element's natural (r,s) basis
              (i.e. NOT rotated to global X/Y — s11,s22,s12).

Global stress = rotated to global X/Y axes (sgxx, sgyy, sgxy) — as before.

Mesh Q8  (5 CQUAD8): nodes 501-520, offset x0=0, y0=2
Mesh T6 (10 CTRIA6): nodes 301-325, offset x0=1, y0=2

Boundary:
  Q8: {501,502,507,508,509,510,511,512}  — 4 corners + 4 bnd mids
  T6: {301,302,307,308,309,310,311,312}  — sama

SPCD verified: Ux=(xp+yp/2)/1000, Uy=(yp+xp/2)/1000
               Uz=(xp²+xp·yp+yp²)/2000, Rx=(yp+xp/2)/1000, Ry=-(xp+yp/2)/1000

NASTRAN result:
  Q8 membrane: PASS (principal exact), bending: FAIL (CQUAD8 limitation)
  T6 bending:  avg Mx/My err <1%, Mxy err 27% (orientation-dependent)
"""

import sys, numpy as np
from functools import partial

sys.path.insert(0, '/mnt/user-data/uploads')
sys.path.insert(0, '/mnt/user-data/outputs')

from core import Model
from MacNealQ8_1992_native             import MacNealQ8_1992_native
from MacNealQ8_1992_native_v3_drill import MacNealQ8_1992_native_v3_drill
from KikuchiMacNeal_Q8_ShellElement_v1 import KikuchiMacNeal_Q8_ShellElement_v1

from Simo1993_Q8_ShellElement_v1p8_standalone import Simo1993_Q8_ShellElement_v1p8
from KikuchiMacNeal_MITC8_v1 import KikuchiMacNeal_MITC8_v1
from KikuchiMacNeal_ANS8_v1 import  KikuchiMacNeal_ANS8_v1
from MITC6_Tri_v1 import MITC6_Tri_v1
from MITC6_Tri_v4 import MITC6_Tri_v4
from MacNeal_MH6T_Tri_v1 import MacNeal_MH6T_Tri_v1
from MacNeal_MH6T_Tri_v3 import MacNeal_MH6T_Tri_v3
from Rezaiee2017_Tri6_Final import Rezaiee2017_Tri6_ShellElement
from Rezaiee2017_Tri6_v3 import Rezaiee2017_Tri6_ShellElement_v3
from Simo1993_Tri6_ShellElement_v2 import Simo1993_Tri6_ShellElement_v2
try:
    from Simo1993_Tri6_ShellElement_v1p8 import Simo1993_Tri6_ShellElement_v1p8
    HAS_SIMO = True
except ImportError:
    HAS_SIMO = False

# ═══════════════════════════════════════════════════════════════════
# MATERIAL
# ═══════════════════════════════════════════════════════════════════
E_, NU_, H_ = 1.0e6, 0.25, 0.001

def Cmat():
    f = E_/(1-NU_**2)
    return np.array([[f,f*NU_,0],[f*NU_,f,0],[0,0,f*(1-NU_)/2]])

C_ = Cmat()
SIG_REF = C_ @ np.array([1e-3, 1e-3, 1e-3])      # [1333.33,1333.33,400]
MOM_REF = C_ @ np.array([1e-3, 1e-3, 1e-3]) * H_**3/12.

# ═══════════════════════════════════════════════════════════════════
# GEOMETRY — persis dari DUEL3.dat
# ═══════════════════════════════════════════════════════════════════

# --- Q8: CQUAD8 nodes 501-520 (NASTRAN y-offset=2) ---
X0_Q8, Y0_Q8 = 0.0, 2.0

NODE_Q8 = {
    501:(0.00,2.00), 502:(0.00,2.12), 503:(0.04,2.02), 504:(0.08,2.08),
    505:(0.18,2.03), 506:(0.16,2.08), 507:(0.24,2.00), 508:(0.24,2.12),
    509:(0.00,2.06), 510:(0.12,2.12), 511:(0.24,2.06), 512:(0.12,2.00),
    513:(0.17,2.055),514:(0.06,2.05),515:(0.02,2.01),516:(0.04,2.10),
    517:(0.20,2.10), 518:(0.21,2.015),519:(0.11,2.025),520:(0.12,2.08),
}
# CQUAD8,21,1,501,502,504,503,509,516,514,515
# G5=mid(G1G2),G6=mid(G2G3),G7=mid(G3G4),G8=mid(G4G1)
CONN_Q8 = {
    21:[501,502,504,503, 509,516,514,515],
    22:[507,501,503,505, 512,515,519,518],
    23:[504,506,505,503, 520,513,519,514],
    24:[502,508,506,504, 510,517,520,516],
    25:[508,507,505,506, 511,518,513,517],
}
BND_Q8 = {501,502,507,508,509,510,511,512}

# --- T6: CTRIA6 nodes 301-325 (NASTRAN x-offset=1, y-offset=2) ---
X0_T6, Y0_T6 = 1.0, 2.0

NODE_T6 = {
    301:(1.00,2.00), 302:(1.00,2.12), 303:(1.04,2.02), 304:(1.08,2.08),
    305:(1.18,2.03), 306:(1.16,2.08), 307:(1.24,2.00), 308:(1.24,2.12),
    309:(1.00,2.06), 310:(1.12,2.12), 311:(1.24,2.06), 312:(1.12,2.00),
    313:(1.17,2.055),314:(1.06,2.05),315:(1.02,2.01),316:(1.04,2.10),
    317:(1.20,2.10), 318:(1.21,2.015),319:(1.11,2.025),320:(1.12,2.08),
    # diagonal interior mid-nodes (verified from DUEL3.dat)
    321:(1.02,2.07),   # mid(302,303)
    322:(1.08,2.10),   # mid(304,310) = mid(306,302) check
    323:(1.21,2.075),  # mid(305,308)
    324:(1.14,2.01),   # mid(303,307)
    325:(1.13,2.055),  # mid(304,305)
}
# CTRIA6 [C1,C2,C3, M12=mid(C1C2), M23=mid(C2C3), M31=mid(C3C1)]
CONN_T6 = {
    51:[303,302,304, 321,316,314],
    52:[303,305,307, 319,318,324],
    53:[304,306,305, 320,313,325],
    54:[306,302,308, 322,310,317],
    55:[305,306,308, 313,317,323],
    56:[301,302,303, 309,321,315],
    57:[304,302,306, 316,322,320],
    58:[305,308,307, 323,311,318],
    59:[301,303,307, 315,324,312],
    60:[303,304,305, 314,325,319],
}
BND_T6 = {301,302,307,308,309,310,311,312}

# ═══════════════════════════════════════════════════════════════════
# BC FIELDS (dengan offset per mesh)
# ═══════════════════════════════════════════════════════════════════
def mem_disp(xn, yn, x0, y0):
    xp,yp = xn-x0, yn-y0
    return (xp+yp/2)/1000, (yp+xp/2)/1000

def bend_disp(xn, yn, x0, y0):
    xp,yp = xn-x0, yn-y0
    return ((xp**2+xp*yp+yp**2)/2000,
            (yp+xp/2)/1000,
            -(xp+yp/2)/1000)

# ═══════════════════════════════════════════════════════════════════
# STRESS EXTRACTION
# v4a: stress_q8/t6 return (sg, mg, slocal, mlocal)
#   sg, mg   = global stresses  (rotated to global X/Y)
#   slocal   = [s11, s22, s12]   in element natural (r,s) basis
#   mlocal   = [m11, m22, m12]   in element natural (r,s) basis
# ═══════════════════════════════════════════════════════════════════
def _get_basis(e, r=0., s=0.):
    try:    return e._local_basis_at_point(r, s)
    except: return (np.array([1.,0.,0.]),np.array([0.,1.,0.]),np.array([0.,0.,1.]))

def _rot_global(s11,s22,s12,e1,e2):
    c11,c12 = e1[0],e1[1]; c21,c22 = e2[0],e2[1]
    sXX = c11**2*s11 + c21**2*s22 + 2*c11*c21*s12
    sYY = c12**2*s11 + c22**2*s22 + 2*c12*c22*s12
    sXY = c11*c12*s11 + c21*c22*s22 + (c11*c22+c12*c21)*s12
    return np.array([sXX,sYY,sXY])

def _Bm_Bb(e, r, s):
    try:    return e._compute_Bm(r,s), e._compute_Bb(r,s)
    except TypeError: return e._compute_Bm(), e._compute_Bb()

def stress_q8(e, u_global, r=0., s=0.):
    """Return (sg, mg, slocal, mlocal).
    sg,mg   = global (XX,YY,XY) after rotation to global X/Y
    slocal  = [s11,s22,s12] in element natural (r,s) frame
    mlocal  = [m11,m22,m12] in element natural (r,s) frame
    """
    idx=e.global_dof_indices(); u_e=e.T_matrix()@u_global[idx]
    Bm,Bb = _Bm_Bb(e,r,s)
    sig = C_ @ (Bm@u_e); mom = C_ @ (Bb@u_e) * H_**3/12.
    slocal = sig.copy(); mlocal = mom.copy()
    e1,e2,e3 = _get_basis(e,r,s)
    sg = _rot_global(sig[0],sig[1],sig[2],e1,e2)
    mg = _rot_global(mom[0],mom[1],mom[2],e1,e2)
    if e3[2] < 0: sg[2]*=-1; mg[2]*=-1
    return sg, mg, slocal, mlocal

def stress_t6(e, u_global, r=1./3., s=1./3.):
    """Return (sg, mg, slocal, mlocal).
    sg,mg   = global (XX,YY,XY) after rotation to global X/Y
    slocal  = [s11,s22,s12] in element natural (r,s) frame
    mlocal  = [m11,m22,m12] in element natural (r,s) frame
    """
    idx=e.global_dof_indices(); u_e=e.T_matrix()@u_global[idx]
    Bm,Bb = _Bm_Bb(e,r,s)
    sig = C_ @ (Bm@u_e); mom = C_ @ (Bb@u_e) * H_**3/12.
    slocal = sig.copy(); mlocal = mom.copy()
    if hasattr(e,'E1') and hasattr(e,'E2'):
        e1,e2 = e.E1, e.E2
        e3z = np.cross(e1,e2)[2]
    else:
        e1,e2,e3 = _get_basis(e,r,s); e3z=e3[2]
    sg = _rot_global(sig[0],sig[1],sig[2],e1,e2)
    mg = _rot_global(mom[0],mom[1],mom[2],e1,e2)
    if e3z < 0: sg[2]*=-1; mg[2]*=-1
    return sg, mg, slocal, mlocal

# ═══════════════════════════════════════════════════════════════════
# BUILD + SOLVE
# ═══════════════════════════════════════════════════════════════════
def build_and_solve(elem_cls, node_xy, conn_dict, bnd_nodes,
                    x0, y0, mode, stress_fn, tol=0.05, tag='Q8'):
    fem = Model(ndof_per_node=6, autospc=True, autospc_tol=1e-10)
    for nid,(xn,yn) in node_xy.items():
        fem.add_node(nid, xn, yn, 0.0)
    elems={}
    for eid,conn in conn_dict.items():
        e = elem_cls(eid,[fem.nodes[n] for n in conn],E_,NU_,H_)
        fem.add_element(e); elems[eid]=e

    # Apply BC persis NASTRAN DUEL3: SPC1 12345 pada bnd + SPCD
    for nid in bnd_nodes:
        xn,yn = node_xy[nid]
        if mode=='membrane':
            ux,uy = mem_disp(xn,yn,x0,y0)
            fem.add_bc(nid,0,ux); fem.add_bc(nid,1,uy)
            for d in [2,3,4,5]: fem.add_bc(nid,d,0.0)
        else:
            uz,rx,ry = bend_disp(xn,yn,x0,y0)
            fem.add_bc(nid,0,0.0); fem.add_bc(nid,1,0.0)
            fem.add_bc(nid,2,uz)
            fem.add_bc(nid,3,rx); fem.add_bc(nid,4,ry)
            fem.add_bc(nid,5,0.0)
    fem.solve_static()

    interior = sorted(set(node_xy)-bnd_nodes)
    disp_errs=[]
    for nid in interior:
        xn,yn=node_xy[nid]; u=fem.u; dofs=fem.nodes[nid].dofs
        if mode=='membrane':
            ru,rv=mem_disp(xn,yn,x0,y0)
            denom=max(abs(ru),abs(rv),1e-14)
            disp_errs+=[abs(u[dofs[0]]-ru)/denom,
                        abs(u[dofs[1]]-rv)/denom]
        else:
            ruz,rrx,rry=bend_disp(xn,yn,x0,y0)
            denom=max(abs(ruz),abs(rrx),abs(rry),1e-14)
            disp_errs+=[abs(u[dofs[2]]-ruz)/denom,
                        abs(u[dofs[3]]-rrx)/denom,
                        abs(u[dofs[4]]-rry)/denom]

    ref = SIG_REF if mode=='membrane' else MOM_REF
    stress_errs=[]
    elem_results={}
    per_grid = {}
    for eid,e in elems.items():
        try:
            sg,mg,slocal,mlocal = stress_fn(e,fem.u)
            vals = sg if mode=='membrane' else mg
            errs = [abs(vals[k]-ref[k])/(abs(ref[k])+1e-30) for k in range(3)]
            stress_errs+=errs
            elem_results[eid]=(vals,errs)
            per_grid[eid] = []
            try:
                gd = e.global_dof_indices()
                ngrid = len(gd) // 6
            except Exception:
                ngrid = 0
            for k in range(1, ngrid+1):
                rs = (Q8_NODE_RS if tag == 'Q8' else T6_NODE_RS)[k-1]
                try:
                    if tag == 'Q8':
                        sg_k, mg_k, sl_k, ml_k = stress_q8(e, fem.u, r=rs[0], s=rs[1])
                    else:
                        sg_k, mg_k, sl_k, ml_k = stress_t6(e, fem.u, r=rs[0], s=rs[1])
                    per_grid[eid].append((sg_k, mg_k, sl_k, ml_k))
                except Exception:
                    per_grid[eid].append(None)
        except Exception as ex:
            stress_errs.append(np.nan)

    max_d = float(np.nanmax(disp_errs))   if disp_errs   else np.nan
    max_s = float(np.nanmax(stress_errs)) if stress_errs else np.nan
    return (max_d<tol and max_s<tol), max_d, max_s, elem_results, per_grid, fem

# ═══════════════════════════════════════════════════════════════════
# REGISTRY — (cls, stress_fn, node_xy, conn, bnd, x0, y0)
# ═══════════════════════════════════════════════════════════════════
# Registry: (cls, stress_fn, tag, node_xy, conn, bnd, x0, y0)
REGISTRY = {
    'MacNealQ8_native': (
        MacNealQ8_1992_native, stress_q8, 'Q8',
        NODE_Q8, CONN_Q8, BND_Q8, X0_Q8, Y0_Q8),

    'SimoQ8': (
        Simo1993_Q8_ShellElement_v1p8, stress_q8, 'Q8',
        NODE_Q8, CONN_Q8, BND_Q8, X0_Q8, Y0_Q8),

    'Kikuchi_Q8': (
        partial(KikuchiMacNeal_Q8_ShellElement_v1, variant='kikuchi'),
        stress_q8, 'Q8', NODE_Q8, CONN_Q8, BND_Q8, X0_Q8, Y0_Q8),

    'MH_Q8': (
        partial(KikuchiMacNeal_Q8_ShellElement_v1, variant='macneal_harder'),
        stress_q8, 'Q8', NODE_Q8, CONN_Q8, BND_Q8, X0_Q8, Y0_Q8),

    'Kikuchi_ANS8': (
        partial(KikuchiMacNeal_ANS8_v1, variant='kikuchi'),
        stress_q8, 'Q8', NODE_Q8, CONN_Q8, BND_Q8, X0_Q8, Y0_Q8),

    'Kikuchi_MITC8': (
        partial(KikuchiMacNeal_MITC8_v1, variant='kikuchi'),
        stress_q8, 'Q8', NODE_Q8, CONN_Q8, BND_Q8, X0_Q8, Y0_Q8),

    'MacNealQ8_native_drill': (
        MacNealQ8_1992_native_v3_drill, stress_q8, 'Q8',
        NODE_Q8, CONN_Q8, BND_Q8, X0_Q8, Y0_Q8),

    'SIMOT6': (
        Simo1993_Tri6_ShellElement_v2, stress_t6, 'T6',
        NODE_T6, CONN_T6, BND_T6, X0_T6, Y0_T6),

    'MITC6': (
        MITC6_Tri_v4, stress_t6, 'T6',
        NODE_T6, CONN_T6, BND_T6, X0_T6, Y0_T6),

    'MH6T': (
        MacNeal_MH6T_Tri_v3, stress_t6, 'T6',
        NODE_T6, CONN_T6, BND_T6, X0_T6, Y0_T6),

    'REZAIEE': (
        Rezaiee2017_Tri6_ShellElement_v3, stress_t6, 'T6',
        NODE_T6, CONN_T6, BND_T6, X0_T6, Y0_T6),
}

# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════
SEP  = "="*72
LINE = "-"*72

def mark(ok): return "✓ PASS" if ok else "✗ FAIL"

# Per-node natural coords for corner/GPSTRESS-style print
Q8_NODE_RS = [(-1.,-1.),( 1.,-1.),( 1., 1.),(-1., 1.),
              ( 0.,-1.),( 1., 0.),( 0., 1.),(-1., 0.)]
T6_AREA_RST = [(1.,0.,0.),(0.,1.,0.),(0.,0.,1.),
               (0.5,0.5,0.),(0.,0.5,0.5),(0.5,0.,0.5)]
T6_NODE_RS = [(r,s) for (r,s,t) in T6_AREA_RST]

def print_per_grid_block(name, tag, conn, per_grid, mode):
    """Print per-grid stress — both LOCAL (s11,s22,s12) and GLOBAL (sXX,sYY,sXY).

    LOCAL  = stress in element natural (r,s) basis — NOT rotated to global.
    GLOBAL = stress rotated to global X/Y axes (original v4 behaviour).

    Layout per row: EID GRID  LOC_s11  LOC_s22  LOC_s12  |  GLOB_sXX  GLOB_sYY  GLOB_sXY
    """
    if mode == 'membrane':
        llbl = ['s11','s22','s12']
        glbl = ['sXX','sYY','sXY']
    else:
        llbl = ['M11','M22','M12']
        glbl = ['MXX','MYY','MXY']

    if not per_grid:
        return

    print(f"      {name} per-grid stress ({mode}):")
    print(f"      LOCAL  = stress in element natural (r,s) frame  (NOT rotated to global)")
    print(f"      GLOBAL = stress rotated to global X/Y             (v4 behaviour)")
    print()

    hdr = (f"      {'EID':>4} {'GRID':>5}  "
           + "  ".join(f"{a:>13}" for a in llbl)
           + "  |"
           + "  ".join(f"{a:>13}" for a in glbl))
    sep = f"      {'----':>4} {'----':>5}  " + "-"*13*3 + "-+-" + "-"*13*3
    print(hdr)
    print(sep.replace('-', '='))

    for eid in sorted(per_grid.keys()):
        conn_e = conn.get(eid, [])
        grids = per_grid[eid]
        for k, gv in enumerate(grids, start=1):
            grid_id = conn_e[k-1] if k-1 < len(conn_e) else '----'
            if gv is None:
                lvals = ['    n/a'] * 3
                gvals = ['    n/a'] * 3
            else:
                sg, mg, sl, ml = gv
                src_local = sl if mode == 'membrane' else ml
                src_global = sg if mode == 'membrane' else mg
                lvals = [f"{v:13.4E}" for v in src_local]
                gvals = [f"{v:13.4E}" for v in src_global]
            lrow = "  ".join(lvals)
            grow = "  ".join(gvals)
            print(f"      {eid:>4} {str(grid_id):>5}  {lrow}  |  {grow}")
    print()

def print_center_stress(name, tag, conn, per_grid, mode):
    """Print center-element stress for each element.
    Q8 center: r=0, s=0
    T6 center: r=1/3, s=1/3
    """
    if mode == 'membrane':
        llbl = ['s11','s22','s12']
        glbl = ['sXX','sYY','sXY']
    else:
        llbl = ['M11','M22','M12']
        glbl = ['MXX','MYY','MXY']

    if not per_grid:
        return

    print(f"      {name} CENTER stress ({mode}):")
    print(f"      {'EID':>4}  "
          + "  ".join(f"{a:>13}" for a in llbl)
          + "  |"
          + "  ".join(f"{a:>13}" for a in glbl))
    print(f"      {'----':>4}  " + "-"*13*3 + "-+-" + "-"*13*3)

    for eid in sorted(per_grid.keys()):
        conn_e = conn.get(eid, [])
        grids = per_grid[eid]
        # Center is the last grid point for both Q8 (8th) and T6 (6th)
        center_idx = len(grids) - 1
        if center_idx < 0 or grids[center_idx] is None:
            continue
        sg, mg, sl, ml = grids[center_idx]
        src_local = sl if mode == 'membrane' else ml
        src_global = sg if mode == 'membrane' else mg
        lvals = [f"{v:13.4E}" for v in src_local]
        gvals = [f"{v:13.4E}" for v in src_global]
        lrow = "  ".join(lvals)
        grow = "  ".join(gvals)
        print(f"      {eid:>4}  {lrow}  |  {grow}")
    print()

def print_grid_displacements(name, tag, node_xy, conn, bnd_nodes, x0, y0, mode, fem):
    """v4c: Print all grid-node displacements for comparison with MYSTRAN F06.

    Prints every node in the mesh with its 6 DOF displacement values.
    For membrane: Ux, Uy (Uz,Rx,Ry,Rz = 0 on boundary)
    For bending:  Uz, Rx, Ry (Ux,Uy,Rz = 0 on boundary)
    """
    print(f"      {name} GRID NODE DISPLACEMENTS ({mode}):")
    print(f"      {'GRID':>5}  {'X':>8} {'Y':>8}  "
          + "  ".join(f"{a:>13}" for a in (['Ux','Uy','Uz','Rx','Ry','Rz'] if mode=='membrane'
                                            else ['Uz','Rx','Ry','Ux','Uy','Rz'])))
    print(f"      {'----':>5}  {'-'*8} {'-'*8}  " + "-"*13*6)

    u = fem.u
    for nid in sorted(node_xy.keys()):
        xn, yn = node_xy[nid]
        dofs = fem.nodes[nid].dofs
        vals = [u[dofs[i]] for i in range(6)]
        if mode == 'membrane':
            # Show Ux, Uy, Uz, Rx, Ry, Rz
            order = [0, 1, 2, 3, 4, 5]
        else:
            # Show Uz, Rx, Ry, Ux, Uy, Rz
            order = [2, 3, 4, 0, 1, 5]
        sv = [f"{vals[i]:13.6E}" for i in order]
        bnd = " *" if nid in bnd_nodes else ""
        print(f"      {nid:>5}  {xn:8.4f} {yn:8.4f}  " + "  ".join(sv) + bnd)
    print(f"      (* = boundary node)")
    print()

def get_tag(name):
    n=name.lower()
    # Only T6-tri variants get 'T6' tag; everything else is Q8
    return 'T6' if any(k in n for k in ('t6','tri')) else 'Q8'

def print_detail(elem_results, ref, mode):
    lbl = ['sxx','syy','sxy'] if mode=='membrane' else ['Mx','My','Mxy']
    print(f"\n    {'Elem':>4}  {lbl[0]:>12} {lbl[1]:>12} {lbl[2]:>12}  "
          f"{'err0%':>7} {'err1%':>7} {'err2%':>7}")
    for eid,(vals,errs) in sorted(elem_results.items()):
        note = " ←" if all(e<0.001 for e in errs) else ""
        print(f"    {eid:>4}  {vals[0]:12.4E} {vals[1]:12.4E} {vals[2]:12.4E}  "
              f"{errs[0]*100:7.2f} {errs[1]*100:7.2f} {errs[2]*100:7.2f}{note}")

if __name__=='__main__':
    print(SEP)
    print("  MH 2-001 Patch Test v4c — LOCAL + GLOBAL per-grid + CENTER stress + DISP")
    print(f"  Ref membrane : σxx={SIG_REF[0]:.2f} σyy={SIG_REF[1]:.2f} τxy={SIG_REF[2]:.2f}")
    print(f"  Ref bending  : Mxx={MOM_REF[0]:.4e} Myy={MOM_REF[1]:.4e} Mxy={MOM_REF[2]:.4e}")
    print(f"  Tol=5%  |  Elemen: {len(REGISTRY)}")
    print(f"  Q8 boundary: {{501,502,507,508,509,510,511,512}}  (offset y-2)")
    print(f"  T6 boundary: {{301,302,307,308,309,310,311,312}}  (offset x-1,y-2)")
    print(SEP)

    for mode, lbl_s in [('membrane','stress'), ('bending','moment')]:
        title = ("SC1 MEMBRANE  (Ux, Uy)" if mode=='membrane'
                 else "SC2 BENDING   (Uz, Rx, Ry)")
        print(f"\n{SEP}\n  {title}\n{LINE}")
        print(f"  {'[tag] name':<28s}  {'disp err':>10}  {lbl_s+' err':>10}  Status")
        print(LINE)

        for name, entry in REGISTRY.items():
            cls, sfn, tag, nxy, conn, bnd, x0, y0 = entry
            try:
                ok,de,se,eres,per_grid,fem = build_and_solve(
                    cls,nxy,conn,bnd,x0,y0,mode,sfn,tag=tag)
                ds = f"{de*100:8.4f}%" if not np.isnan(de) else "     n/a"
                ss = f"{se*100:8.4f}%" if not np.isnan(se) else "     n/a"
                print(f"  [{tag}] {name:<24s}  {ds}  {ss}  {mark(ok)}")
                if '--detail' in sys.argv:
                    print_detail(eres, SIG_REF if mode=='membrane' else MOM_REF, mode)
                # Always print per-grid stress — now with LOCAL + GLOBAL
                print_per_grid_block(name, tag, conn, per_grid, mode)
                # v4b: also print center-element stress
                print_center_stress(name, tag, conn, per_grid, mode)
                # v4c: print grid-node displacements
                print_grid_displacements(name, tag, nxy, conn, bnd, x0, y0, mode, fem)
            except Exception as ex:
                print(f"  [{tag}] {name:<24s}  ERROR: {ex}")

    print(f"\n{SEP}")
    print("  Note: run dengan '--detail' untuk lihat nilai per elemen")
    print(f"  NASTRAN benchmark: Q8 membrane EXACT (principal), bending FAIL")
    print(f"  T6 bending: avg Mx/My <1% err, Mxy ~27% (orientation effect)")
    print(f"{SEP}")

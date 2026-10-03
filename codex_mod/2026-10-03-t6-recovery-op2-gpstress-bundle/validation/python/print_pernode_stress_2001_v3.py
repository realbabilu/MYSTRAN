"""
print_pernode_stress_2001_v3.py
=================================
Berbasis test_patch2001_q8_t6_v4.py / print_pernode_stress_2001_v2.py
(mesh & BC identik). Untuk pembandingan dengan MYSTRAN / NASTRAN.

Output per elemen (CTR + tiap node corner & mid):
  1. STRESS LOKAL   : dalam frame elemen (satu frame per elemen, dievaluasi di pusat)
  2. STRESS GLOBAL  : dirotasi ke XY global, per elemen per node
  3. GP-AVERAGE     : rata-rata sederhana stress global di tiap grid dari semua
                      elemen yang berbagi grid (pendekatan "surface GPSTRESS")
  4. DISPLACEMENT   : global (ux uy uz rx ry rz) dan lokal-elemen

Untuk mode bending juga dicetak stress serat Z1=-h/2 dan Z2=+h/2
(sigma = sigma_membran -/+ 6M/h^2). Tanda tergantung konvensi Bb elemen.

Usage:
    python print_pernode_stress_2001_v3.py                 # semua elemen
    python print_pernode_stress_2001_v3.py --detail        # + error vs medan eksak
    python print_pernode_stress_2001_v3.py --only MITC6    # satu elemen saja
    python print_pernode_stress_2001_v3.py --noflip        # matikan flip shear global (e3z<0)
    python print_pernode_stress_2001_v3.py --mode bending  # satu mode saja
"""
import sys
import inspect
import numpy as np
from functools import partial

sys.path.insert(0, '.')

from core import Model
from MacNealQ8_1992_native             import MacNealQ8_1992_native
from MacNealQ8_1992_native_v3_drill import MacNealQ8_1992_native_v3_drill
from KikuchiMacNeal_Q8_ShellElement_v1 import KikuchiMacNeal_Q8_ShellElement_v1
from Simo1993_Q8_ShellElement_v1p8_standalone import Simo1993_Q8_ShellElement_v1p8
from KikuchiMacNeal_MITC8_v1 import KikuchiMacNeal_MITC8_v1
from KikuchiMacNeal_ANS8_v1 import KikuchiMacNeal_ANS8_v1
from MITC6_Tri_v1 import MITC6_Tri_v1
from MITC6_Tri_v4 import MITC6_Tri_v4
from Rezaiee2017_Tri6_Final import Rezaiee2017_Tri6_ShellElement
from Rezaiee2017_Tri6_v3 import Rezaiee2017_Tri6_ShellElement_v3
from  MacNeal_MH6T_Tri_v3 import MacNeal_MH6T_Tri_v3
from Simo1993_Tri6_ShellElement_v2 import Simo1993_Tri6_ShellElement_v2
try:
    from Simo1993_Tri6_ShellElement_v1p8 import Simo1993_Tri6_ShellElement_v1p8
    HAS_SIMO_T6 = True
except ImportError:
    HAS_SIMO_T6 = False

# ═══════════════════════════════════════════════════════════════════
# MATERIAL
# ═══════════════════════════════════════════════════════════════════
E_, NU_, H_ = 1.0e6, 0.25, 0.001

def Cmat():
    f = E_/(1-NU_**2)
    return np.array([[f,f*NU_,0],[f*NU_,f,0],[0,0,f*(1-NU_)/2]])
C_ = Cmat()
SIG_REF = C_ @ np.array([1e-3, 1e-3, 1e-3])
MOM_REF = C_ @ np.array([1e-3, 1e-3, 1e-3]) * H_**3/12.

# ═══════════════════════════════════════════════════════════════════
# GEOMETRY — persis dari DUEL3.dat
# ═══════════════════════════════════════════════════════════════════
X0_Q8, Y0_Q8 = 0.0, 2.0
NODE_Q8 = {
    501:(0.00,2.00), 502:(0.00,2.12), 503:(0.04,2.02), 504:(0.08,2.08),
    505:(0.18,2.03), 506:(0.16,2.08), 507:(0.24,2.00), 508:(0.24,2.12),
    509:(0.00,2.06), 510:(0.12,2.12), 511:(0.24,2.06), 512:(0.12,2.00),
    513:(0.17,2.055),514:(0.06,2.05),515:(0.02,2.01),516:(0.04,2.10),
    517:(0.20,2.10), 518:(0.21,2.015),519:(0.11,2.025),520:(0.12,2.08),
}
CONN_Q8 = {
    21:[501,502,504,503, 509,516,514,515],
    22:[507,501,503,505, 512,515,519,518],
    23:[504,506,505,503, 520,513,519,514],
    24:[502,508,506,504, 510,517,520,516],
    25:[508,507,505,506, 511,518,513,517],
}
BND_Q8 = {501,502,507,508,509,510,511,512}

X0_T6, Y0_T6 = 1.0, 2.0
NODE_T6 = {
    301:(1.00,2.00), 302:(1.00,2.12), 303:(1.04,2.02), 304:(1.08,2.08),
    305:(1.18,2.03), 306:(1.16,2.08), 307:(1.24,2.00), 308:(1.24,2.12),
    309:(1.00,2.06), 310:(1.12,2.12), 311:(1.24,2.06), 312:(1.12,2.00),
    313:(1.17,2.055),314:(1.06,2.05),315:(1.02,2.01),316:(1.04,2.10),
    317:(1.20,2.10), 318:(1.21,2.015),319:(1.11,2.025),320:(1.12,2.08),
    321:(1.02,2.07),   322:(1.08,2.10), 323:(1.21,2.075),
    324:(1.14,2.01),   325:(1.13,2.055),
}
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
# BC / MEDAN EKSAK
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
# NATURAL COORDS
# ═══════════════════════════════════════════════════════════════════
Q8_NODE_RS = [(-1.,-1.),( 1.,-1.),( 1., 1.),(-1., 1.),
              ( 0.,-1.),( 1., 0.),( 0., 1.),(-1., 0.)]
Q8_NODE_LBL = ['C1','C2','C3','C4','M5','M6','M7','M8']
Q8_CENTER = (0., 0.)

T6_NODE_RS = [(1.,0.),(0.,1.),(0.,0.),
              (0.5,0.5),(0.,0.5),(0.5,0.)]
T6_NODE_LBL = ['C1','C2','C3','M12','M23','M31']
T6_CENTER = (1./3., 1./3.)          # centroid (BUKAN (0,0) = corner C3)

def q8_N(r, s):
    N = np.zeros(8)
    N[0] = 0.25*(1-r)*(1-s)*(-r-s-1)
    N[1] = 0.25*(1+r)*(1-s)*( r-s-1)
    N[2] = 0.25*(1+r)*(1+s)*( r+s-1)
    N[3] = 0.25*(1-r)*(1+s)*(-r+s-1)
    N[4] = 0.5*(1-r**2)*(1-s)
    N[5] = 0.5*(1+r)*(1-s**2)
    N[6] = 0.5*(1-r**2)*(1+s)
    N[7] = 0.5*(1-r)*(1-s**2)
    return N

def t6_N(r, s):
    t = 1.0 - r - s
    return np.array([r*(2*r-1), s*(2*s-1), t*(2*t-1), 4*r*s, 4*s*t, 4*t*r])

# ═══════════════════════════════════════════════════════════════════
# FRAME & ROTASI
# ═══════════════════════════════════════════════════════════════════
_WARNED = set()

# Flip tanda shear saat e3z<0 (warisan v2/v4). HANYA berlaku untuk output GLOBAL.
# Output LOKAL tidak pernah di-flip. Matikan dgn --noflip bila frame elemen (E1,E2)
# sudah kanan-tangan terhadap normalnya (E2 = E3 x E1).
FLIP_SHEAR_IF_NEG_E3 = True

def get_basis(e, r, s):
    """(e1,e2,e3) frame lokal elemen di (r,s). Fallback ke sumbu global
    hanya jika elemen memang tidak menyediakan basis (dengan peringatan)."""
    if hasattr(e, 'E1') and hasattr(e, 'E2'):
        e1 = np.asarray(e.E1, float); e2 = np.asarray(e.E2, float)
        return e1, e2, np.cross(e1, e2)
    if hasattr(e, '_local_basis_at_point'):
        return e._local_basis_at_point(r, s)
    key = type(e).__name__
    if key not in _WARNED:
        _WARNED.add(key)
        print(f"  [WARN] {key}: tidak ada E1/E2 maupun _local_basis_at_point; "
              f"memakai sumbu global sebagai frame lokal.")
    return np.array([1.,0.,0.]), np.array([0.,1.,0.]), np.array([0.,0.,1.])

def rot_to_global(s11, s22, s12, e1, e2):
    c11,c12 = e1[0],e1[1]; c21,c22 = e2[0],e2[1]
    return np.array([
        c11**2*s11 + c21**2*s22 + 2*c11*c21*s12,
        c12**2*s11 + c22**2*s22 + 2*c12*c22*s12,
        c11*c12*s11 + c21*c22*s22 + (c11*c22+c12*c21)*s12])

def rot_to_local(sXX, sYY, sXY, e1, e2):
    c11,c12 = e1[0],e1[1]; c21,c22 = e2[0],e2[1]
    return np.array([
        c11**2*sXX + c12**2*sYY + 2*c11*c12*sXY,
        c21**2*sXX + c22**2*sYY + 2*c21*c22*sXY,
        c11*c21*sXX + c12*c22*sYY + (c11*c22+c12*c21)*sXY])

def get_Bm_Bb(e, r, s):
    """Pilih signature berdasarkan inspect (tidak menelan TypeError internal)."""
    def call(fn):
        n = len([p for p in inspect.signature(fn).parameters.values()
                 if p.default is p.empty])
        return fn(r, s) if n >= 2 else fn()
    return call(e._compute_Bm), call(e._compute_Bb)

# ═══════════════════════════════════════════════════════════════════
# EVALUASI STRESS DI SATU TITIK
# ═══════════════════════════════════════════════════════════════════
def eval_point(e, u_global, r, s, ctr):
    """Return dict: sig_loc, mom_loc (frame elemen di pusat),
    sig_glb, mom_glb, frame angle (deg)."""
    idx = e.global_dof_indices()
    u_e = e.T_matrix() @ u_global[idx]
    Bm, Bb = get_Bm_Bb(e, r, s)
    sig = C_ @ (Bm @ u_e)
    mom = C_ @ (Bb @ u_e) * H_**3/12.

    e1, e2, e3 = get_basis(e, r, s)
    sg0 = rot_to_global(*sig, e1, e2)           # global, tanpa flip
    mg0 = rot_to_global(*mom, e1, e2)

    # frame elemen tunggal (pusat). LOKAL = nilai di frame ini, TANPA flip tanda apa pun.
    c1, c2, c3 = get_basis(e, ctr[0], ctr[1])
    sl = rot_to_local(*sg0, c1, c2)
    ml = rot_to_local(*mg0, c1, c2)

    sg = sg0.copy(); mg = mg0.copy()
    if FLIP_SHEAR_IF_NEG_E3 and e3[2] < 0:      # hanya global
        sg[2] *= -1; mg[2] *= -1

    ang = np.degrees(np.arctan2(c1[1], c1[0]))
    return dict(sig_loc=sl, mom_loc=ml, sig_glb=sg, mom_glb=mg, ang=ang,
                sig_pt=sig, mom_pt=mom)

def fibers(sig, mom):
    """Stress serat: Z1=-h/2, Z2=+h/2."""
    k = 6.0/H_**2
    return sig - k*mom, sig + k*mom

# ═══════════════════════════════════════════════════════════════════
# BUILD + SOLVE
# ═══════════════════════════════════════════════════════════════════
def build_and_solve(elem_cls, node_xy, conn_dict, bnd_nodes, x0, y0, mode):
    fem = Model(ndof_per_node=6, autospc=True, autospc_tol=1e-10)
    for nid,(xn,yn) in node_xy.items():
        fem.add_node(nid, xn, yn, 0.0)
    elems = {}
    for eid,conn in conn_dict.items():
        e = elem_cls(eid, [fem.nodes[n] for n in conn], E_, NU_, H_)
        fem.add_element(e); elems[eid] = e

    for nid in bnd_nodes:
        xn,yn = node_xy[nid]
        if mode == 'membrane':
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
    return fem, elems

# ═══════════════════════════════════════════════════════════════════
# PRINT HELPERS
# ═══════════════════════════════════════════════════════════════════
def fmt(vals, w=12, p=4, sci=True):
    f = f"{{:{w}.{p}E}}" if sci else f"{{:{w}.{p}f}}"
    return "  ".join(f.format(v) for v in vals)

def collect_element(e, fem, conn, mode, tag):
    """Hitung semua titik (CTR + nodes) sekali; dipakai print & averaging."""
    is_q8 = (tag == 'Q8')
    node_rs  = Q8_NODE_RS  if is_q8 else T6_NODE_RS
    node_lbl = Q8_NODE_LBL if is_q8 else T6_NODE_LBL
    ctr      = Q8_CENTER   if is_q8 else T6_CENTER
    Nfun     = q8_N        if is_q8 else t6_N

    u_global = fem.u
    idx = e.global_dof_indices()
    u_loc = e.T_matrix() @ u_global[idx]
    n = len(node_rs)
    u_loc_n = u_loc.reshape(n, 6)
    u_glb_n = np.array([u_global[fem.nodes[nid].dofs][:6] for nid in conn])

    pts = []
    for lbl, (r, s), k in [('CTR', ctr, None)] + \
            [(l, rs, i) for i, (l, rs) in enumerate(zip(node_lbl, node_rs))]:
        res = eval_point(e, u_global, r, s, ctr)
        if k is None:
            N = Nfun(r, s)
            ug = N @ u_glb_n; ul = N @ u_loc_n; nid = None
        else:
            ug = u_glb_n[k]; ul = u_loc_n[k]; nid = conn[k]
        res.update(lbl=lbl, nid=nid, u_glb=ug, u_loc=ul)
        pts.append(res)
    return pts

def print_element(name, eid, pts, mode):
    if mode == 'membrane':
        ls = ['s11','s22','s12']; gs = ['sXX','sYY','sXY']
        key_l, key_g = 'sig_loc', 'sig_glb'
    else:
        ls = ['M11','M22','M12']; gs = ['MXX','MYY','MXY']
        key_l, key_g = 'mom_loc', 'mom_glb'
    dl = ['ux','uy','uz','rx','ry','rz']

    print(f"\n    --- {name}  elemen {eid}   (sudut e1 vs X = {pts[0]['ang']:.4f} deg) ---")

    # 1. LOKAL
    print(f"    [LOCAL {'stress' if mode=='membrane' else 'moment'}, frame elemen]")
    hdr = f"    {'pt':>4} {'grid':>5}  " + "  ".join(f"{a:>12}" for a in ls)
    if mode == 'bending':
        hdr += "  |  " + "  ".join(f"{a:>12}" for a in ['Z1:s11','Z1:s22','Z1:s12',
                                                      'Z2:s11','Z2:s22','Z2:s12'])
    print(hdr)
    for p in pts:
        row = f"    {p['lbl']:>4} {str(p['nid'] or '---'):>5}  " + fmt(p[key_l])
        if mode == 'bending':
            z1, z2 = fibers(p['sig_loc'], p['mom_loc'])
            row += "  |  " + fmt(z1) + "  " + fmt(z2)
        print(row)

    # 2. GLOBAL
    print(f"    [GLOBAL {'stress' if mode=='membrane' else 'moment'}, per-elemen]")
    hdr = f"    {'pt':>4} {'grid':>5}  " + "  ".join(f"{a:>12}" for a in gs)
    if mode == 'bending':
        hdr += "  |  " + "  ".join(f"{a:>12}" for a in ['Z1:sXX','Z1:sYY','Z1:sXY',
                                                      'Z2:sXX','Z2:sYY','Z2:sXY'])
    print(hdr)
    for p in pts:
        row = f"    {p['lbl']:>4} {str(p['nid'] or '---'):>5}  " + fmt(p[key_g])
        if mode == 'bending':
            z1, z2 = fibers(p['sig_glb'], p['mom_glb'])
            row += "  |  " + fmt(z1) + "  " + fmt(z2)
        print(row)

    # 3. DISPLACEMENT
    print(f"    [DISPLACEMENT]  global | lokal-elemen")
    print(f"    {'pt':>4} {'grid':>5}  " + "  ".join(f"{a:>11}" for a in dl)
          + "  |  " + "  ".join(f"{a:>11}" for a in dl))
    for p in pts:
        print(f"    {p['lbl']:>4} {str(p['nid'] or '---'):>5}  "
              + fmt(p['u_glb'], 11, 5) + "  |  " + fmt(p['u_loc'], 11, 5))

def print_gp_average(name, acc, mode, ref):
    """Rata-rata stress GLOBAL per grid dari semua elemen yang berbagi grid."""
    gs = ['sXX','syy','sXY'] if mode == 'membrane' else ['MXX','MYY','MXY']
    gs = [g.upper() if g[0] in 'sM' else g for g in gs]
    print(f"\n    [GP-AVERAGE global {mode}]  {name}  "
          f"(rata-rata sederhana antar elemen; BUKAN algoritma patch MYSTRAN)")
    print(f"    {'grid':>5} {'n':>2}  " + "  ".join(f"{a:>12}" for a in gs)
          + "  |  " + "  ".join(f"{'err%'+a[-2:]:>9}" for a in gs))
    for nid in sorted(acc):
        v = np.mean(acc[nid], axis=0)
        err = [abs(v[k]-ref[k])/(abs(ref[k])+1e-30)*100 for k in range(3)]
        print(f"    {nid:>5} {len(acc[nid]):>2}  " + fmt(v) + "  |  "
              + "  ".join(f"{x:9.3f}" for x in err))

def print_detail(name, fem, node_xy, bnd, x0, y0, mode):
    """Error displacement tiap node vs medan eksak (interior saja)."""
    print(f"\n    [DETAIL displacement error vs eksak, node interior]  {name}")
    errs = []
    for nid in sorted(set(node_xy) - bnd):
        xn, yn = node_xy[nid]; d = fem.nodes[nid].dofs; u = fem.u
        if mode == 'membrane':
            ex = mem_disp(xn, yn, x0, y0); got = (u[d[0]], u[d[1]])
        else:
            ex = bend_disp(xn, yn, x0, y0); got = (u[d[2]], u[d[3]], u[d[4]])
        den = max(max(abs(x) for x in ex), 1e-14)
        errs.append(max(abs(g-x) for g, x in zip(got, ex))/den)
    print(f"      max err = {max(errs)*100:.6f} %   (jumlah node interior = {len(errs)})")

# ═══════════════════════════════════════════════════════════════════
# REGISTRY
# ═══════════════════════════════════════════════════════════════════
def _q8(cls): return (cls, NODE_Q8, CONN_Q8, BND_Q8, X0_Q8, Y0_Q8, 'Q8')
def _t6(cls): return (cls, NODE_T6, CONN_T6, BND_T6, X0_T6, Y0_T6, 'T6')

REGISTRY = {
    'MacNealQ8_native':       _q8(MacNealQ8_1992_native),
    'SimoQ8':                 _q8(Simo1993_Q8_ShellElement_v1p8),
    'Kikuchi_Q8':             _q8(partial(KikuchiMacNeal_Q8_ShellElement_v1, variant='kikuchi')),
    'MH_Q8':                  _q8(partial(KikuchiMacNeal_Q8_ShellElement_v1, variant='macneal_harder')),
    'Kikuchi_ANS8':           _q8(partial(KikuchiMacNeal_ANS8_v1, variant='kikuchi')),
    'Kikuchi_MITC8':          _q8(partial(KikuchiMacNeal_MITC8_v1, variant='kikuchi')),
    'MacNealQ8_native_drill': _q8(MacNealQ8_1992_native_v3_drill),
    'MITC6':                  _t6(MITC6_Tri_v1),
    'MITC6v4':                  _t6(MITC6_Tri_v4),
    'Rezaiee_v1':                  _t6(Rezaiee2017_Tri6_ShellElement),
    'Rezaiee_v3':                  _t6(Rezaiee2017_Tri6_ShellElement_v3),
     'MHT6v3':                  _t6(MacNeal_MH6T_Tri_v3),
     'SimoT6v2':                 _t6(Simo1993_Tri6_ShellElement_v2),
}
if HAS_SIMO_T6:
    REGISTRY['SimoT6'] = _t6(Simo1993_Tri6_ShellElement_v1p8)

# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════
SEP  = "="*100
LINE = "-"*100

def _arg(flag):
    if flag in sys.argv:
        i = sys.argv.index(flag)
        if i+1 < len(sys.argv): return sys.argv[i+1]
    return None

if __name__ == '__main__':
    detail = '--detail' in sys.argv
    if '--noflip' in sys.argv:
        FLIP_SHEAR_IF_NEG_E3 = False
    only   = _arg('--only')
    only_mode = _arg('--mode')

    print(SEP)
    print("  PATCH TEST 2-001 v3 — Local + Global Stress, GP-average, Displacement")
    print(f"  Ref membrane : sigma_xx={SIG_REF[0]:.2f}  sigma_yy={SIG_REF[1]:.2f}  tau_xy={SIG_REF[2]:.2f}")
    print(f"  Ref bending  : Mxx={MOM_REF[0]:.4e}  Myy={MOM_REF[1]:.4e}  Mxy={MOM_REF[2]:.4e}")
    print(f"  Q8: {len(CONN_Q8)} elemen (nodes 501-520) | T6: {len(CONN_T6)} elemen (nodes 301-325)")
    print(f"  Frame lokal = frame elemen di pusat; CTR Q8=(0,0), CTR T6=(1/3,1/3); lokal tidak di-flip")
    print(f"  Flip shear global bila e3z<0: {'AKTIF' if FLIP_SHEAR_IF_NEG_E3 else 'MATI (--noflip)'}")
    print(SEP)

    for mode in ('membrane', 'bending'):
        if only_mode and only_mode != mode:
            continue
        title = ("SC1 MEMBRANE  (Ux, Uy)" if mode == 'membrane'
                 else "SC2 BENDING   (Uz, Rx, Ry)")
        ref = SIG_REF if mode == 'membrane' else MOM_REF
        print(f"\n{SEP}\n  {title}\n{LINE}")

        for name, (cls, nxy, conn, bnd, x0, y0, tag) in REGISTRY.items():
            if only and only != name:
                continue
            try:
                fem, elems = build_and_solve(cls, nxy, conn, bnd, x0, y0, mode)
                print(f"\n  [{tag}] {name}")
                acc = {}
                for eid in sorted(elems):
                    pts = collect_element(elems[eid], fem, conn[eid], mode, tag)
                    print_element(name, eid, pts, mode)
                    key = 'sig_glb' if mode == 'membrane' else 'mom_glb'
                    for p in pts:
                        if p['nid'] is not None:
                            acc.setdefault(p['nid'], []).append(p[key])
                print_gp_average(name, acc, mode, ref)
                if detail:
                    print_detail(name, fem, nxy, bnd, x0, y0, mode)
            except Exception as ex:
                print(f"\n  [{tag}] {name}  ERROR: {ex}")
                import traceback; traceback.print_exc()

    print(f"\n{SEP}\n  Done.\n{SEP}")

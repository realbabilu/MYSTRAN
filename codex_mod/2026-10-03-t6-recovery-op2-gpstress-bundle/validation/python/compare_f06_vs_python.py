"""
compare_f06_vs_python.py
========================
Parse F06 MYSTRAN / NASTRAN (T6, patch test 2-001) dan bandingkan dengan:
  (a) medan EKSAK patch test        -> selalu jalan (tidak perlu modul Python elemen)
  (b) hasil Python v3 (opsional)    -> --python NAMA  (mis. MITC6, SimoT6)

Usage:
  python compare_f06_vs_python.py duel3a_mitc6_1.F06 duel3a_simo_1.F06 duel3a_mh6t_1.F06 \
         --nastran duel3a_t6_nastran.f06
  python compare_f06_vs_python.py duel3a_mitc6_1.F06 --python MITC6      # + Python v3
  python compare_f06_vs_python.py ... > report.txt

Urutan pemeriksaan (hulu -> hilir):
  A. stress lokal elemen di CENTER
  B. stress lokal elemen di GRD (corner)   <- baris "GRD" pada F06 MYSTRAN termodifikasi
  C. GP stress (surface, global, sudah dirata-rata)
Opsi: --normal auto|cw|ccw  (default auto = normal dari urutan node, konvensi NASTRAN/MYSTRAN)

Besaran yang dibandingkan:
  * displacement tiap node vs eksak
  * stress elemen (frame lokal, CENTER): trace(s11+s22) dan principal
    (invarian terhadap rotasi frame -> tidak tergantung konvensi e1/e2)
    + komponen lokal eksak untuk frame n1->n2 (informasi)
  * GP stress (surface, sistem global) vs eksak
  * Uji hipotesis bug: strain RECOV249 direproduksi dari geometri GLOBAL +
    displacement LOKAL (frame tidak konsisten)
"""
import re
import sys
import numpy as np

# ─────────────── material / mesh (sama dgn v3) ───────────────
E_, NU_, H_ = 1.0e6, 0.25, 0.001
_f = E_/(1-NU_**2)
C_ = np.array([[_f,_f*NU_,0],[_f*NU_,_f,0],[0,0,_f*(1-NU_)/2]])
X0, Y0 = 1.0, 2.0
CONN = {51:[303,302,304,321,316,314], 52:[303,305,307,319,318,324],
        53:[304,306,305,320,313,325], 54:[306,302,308,322,310,317],
        55:[305,306,308,313,317,323], 56:[301,302,303,309,321,315],
        57:[304,302,306,316,322,320], 58:[305,308,307,323,311,318],
        59:[301,303,307,315,324,312], 60:[303,304,305,314,325,319]}
XY = {301:(1.00,2.00),302:(1.00,2.12),303:(1.04,2.02),304:(1.08,2.08),
      305:(1.18,2.03),306:(1.16,2.08),307:(1.24,2.00),308:(1.24,2.12),
      309:(1.00,2.06),310:(1.12,2.12),311:(1.24,2.06),312:(1.12,2.00),
      313:(1.17,2.055),314:(1.06,2.05),315:(1.02,2.01),316:(1.04,2.10),
      317:(1.20,2.10),318:(1.21,2.015),319:(1.11,2.025),320:(1.12,2.08),
      321:(1.02,2.07),322:(1.08,2.10),323:(1.21,2.075),324:(1.14,2.01),
      325:(1.13,2.055)}
CORNER_GRIDS = [301,302,303,304,305,306,307,308]
STRAIN_REF = np.array([1e-3,1e-3,1e-3])          # exx, eyy, gxy (uniform)
NORMAL_MODE = 'auto'                              # 'auto' | 'cw' | 'ccw' (override konvensi normal elemen)

def exact_disp(nid, mode):
    xp, yp = XY[nid][0]-X0, XY[nid][1]-Y0
    if mode == 'membrane':
        return np.array([(xp+yp/2)/1e3, (yp+xp/2)/1e3])
    return np.array([(xp**2+xp*yp+yp**2)/2e3, (yp+xp/2)/1e3, -(xp+yp/2)/1e3])

def rot_local(v, e1, e2):
    """[sxx,syy,sxy] global -> local (e1,e2)."""
    c11,c12,c21,c22 = e1[0],e1[1],e2[0],e2[1]
    return np.array([c11**2*v[0]+c12**2*v[1]+2*c11*c12*v[2],
                     c21**2*v[0]+c22**2*v[1]+2*c21*c22*v[2],
                     c11*c21*v[0]+c12*c22*v[1]+(c11*c22+c12*c21)*v[2]])

def elem_normal_z(eid):
    """+1 jika node 1-2-3 berlawanan jarum jam (normal +Z), -1 jika searah jarum jam."""
    if NORMAL_MODE == 'cw':  return -1.0
    if NORMAL_MODE == 'ccw': return 1.0
    P = [np.array(XY[n]) for n in CONN[eid][:3]]
    a = (P[1][0]-P[0][0])*(P[2][1]-P[0][1]) - (P[1][1]-P[0][1])*(P[2][0]-P[0][0])
    return 1.0 if a > 0 else -1.0

def frame_n1n2(eid):
    """e1 = n1->n2, e3 = normal dari urutan node, e2 = e3 x e1 (konvensi NASTRAN/MYSTRAN)."""
    P = [np.array(XY[n]) for n in CONN[eid][:3]]
    e1 = P[1]-P[0]; e1 /= np.linalg.norm(e1)
    return e1, elem_normal_z(eid)*np.array([-e1[1], e1[0]])

def exact_stress_global(mode):
    """membran: sigma; bending: fibre Z2 (+h/2)."""
    s = C_ @ STRAIN_REF
    return s if mode == 'membrane' else s*H_/2.0

def invariants(v):
    a, b, c = v
    m = 0.5*(a+b); r = np.hypot(0.5*(a-b), c)
    return a+b, m+r, m-r

# ─────────────── parser MYSTRAN ───────────────
def _floats(tokens):
    out = []
    for t in tokens:
        try: out.append(float(t))
        except ValueError: return None
    return out

def parse_mystran(path):
    """Mendukung format lama (hanya CENTER) dan baru (CENTER + baris GRD, engineering forces)."""
    L = open(path, errors='ignore').read().replace('\r', '').split('\n')
    D = dict(disp={}, elstr={}, elgrd={}, elfrc={}, gp={}, gpf={}, recov={})
    sc = 1; mode = None; gid = None; cur_eid = None
    recov_count = {}
    for i, l in enumerate(L):
        t = l.split()
        m = re.search(r'OUTPUT FOR SUBCASE\s+(\d+)', l)
        if m: sc = int(m.group(1)); mode = None
        m = re.match(r'^0\s+SUBCASE\s+(\d+)', l)
        if m: sc = int(m.group(1))
        if   'S T R E S S E S   A T   G R I D' in l:                          mode = 'gp'
        elif 'E L E M E N T   S T R E S S E S   I N   L O C A L' in l:        mode = 'elstr'
        elif 'E N G I N E E R I N G   F O R C E S' in l:                      mode = 'elfrc'
        elif 'F O R C E S   A T   G R I D' in l:                              mode = 'gpf'
        elif 'D I S P L A C E M E N T S' in l:                                mode = 'disp'
        if not t: continue

        if mode in (None, 'disp') and len(t) == 8 and t[0].isdigit() and t[1] == '0':
            v = _floats(t[2:])
            if v is not None: D['disp'].setdefault(sc, {})[int(t[0])] = np.array(v)

        elif mode == 'elstr' and t[0].isdigit() and len(t) >= 6 and t[1] in ('CENTER', 'GRD'):
            if t[1] == 'CENTER': z1 = _floats(t[3:6]); key = None
            else:                z1 = _floats(t[4:7]); key = int(t[2])
            z2 = _floats(L[i+1].split()[1:4])
            if z1 is None or z2 is None: continue
            val = (np.array(z1), np.array(z2)); eid = int(t[0])
            if key is None: D['elstr'].setdefault(sc, {})[eid] = val
            else:           D['elgrd'].setdefault(sc, {}).setdefault(eid, {})[key] = val

        elif mode == 'elfrc':
            if t[0].isdigit() and len(t) == 9 and _floats(t[1:]) is not None:
                cur_eid = int(t[0]); D['elfrc'].setdefault(sc, {}).setdefault(cur_eid, {})['CENTER'] = np.array(_floats(t[1:]))
            elif t[0] == 'GRD' and len(t) == 10 and cur_eid is not None and _floats(t[2:]) is not None:
                D['elfrc'].setdefault(sc, {}).setdefault(cur_eid, {})[int(t[1])] = np.array(_floats(t[2:]))

        elif mode == 'gpf' and len(t) == 10 and t[0].isdigit() and t[1].isdigit() and _floats(t[2:]) is not None:
            D['gpf'].setdefault(sc, {})[(int(t[0]), int(t[1]))] = np.array(_floats(t[2:]))

        elif mode == 'gp':
            if len(t) >= 11 and t[0].isdigit() and t[2] in ('Z1', 'Z2', 'MID'):
                gid = int(t[0]); D['gp'].setdefault(sc, {}).setdefault(gid, {})[t[2]] = np.array(_floats(t[3:6]))
            elif len(t) >= 4 and t[0] in ('Z1', 'Z2', 'MID') and gid is not None:
                D['gp'].setdefault(sc, {}).setdefault(gid, {})[t[0]] = np.array(_floats(t[1:4]))

        m = re.match(r'RECOV249 EID/PT\s+(\d+)\s+STR_PT_NUM\s+(\d+)', l)
        if m:
            key = (int(m.group(1)), int(m.group(2)))
            recov_count[key] = recov_count.get(key, 0) + 1
            rs = recov_count[key]; blk = {}
            j = i+1
            while j < len(L) and not L[j].startswith('RECOV249') and L[j].strip():
                tt = L[j].split()
                if tt and tt[0] in ('STRAIN1','STRAIN2','STRAIN3','STRESS1','STRESS2','STRESS3','EM','EB','ET','ZS'):
                    blk[tt[0]] = np.array(_floats(tt[1:]))
                j += 1
            D['recov'].setdefault(rs, {})[key] = blk
    return D

def parse_nastran_disp(path):
    L = open(path, errors='ignore').read().replace('\r', '').split('\n')
    out = {}; sc = 0
    for i, l in enumerate(L):
        if 'D I S P L A C E M E N T   V E C T O R' in l:
            sc += 1; out[sc] = {}
            for m in L[i+1:i+60]:
                if m.startswith('1 ') or 'SUBCASE' in m:      # page break / subcase berikutnya
                    break
                t = m.split()
                if len(t) == 8 and t[0].isdigit() and t[1] == 'G':
                    out[sc][int(t[0])] = np.array(_floats(t[2:]))
    return out

# ─────────────── T6 strain (uji hipotesis) ───────────────
def _dN(r, s):
    t = 1-r-s
    return (np.array([4*r-1,0,-(4*t-1),4*s,-4*s,4*(t-r)]),
            np.array([0,4*s-1,-(4*t-1),4*r,4*(t-s),-4*r]))

def t6_strain(X, U, r=1/3., s=1/3.):
    dr, ds = _dN(r, s)
    J = np.array([[dr@X[:,0], dr@X[:,1]],[ds@X[:,0], ds@X[:,1]]])
    d = np.linalg.solve(J, np.array([dr, ds]))
    return np.array([d[0]@U[:,0], d[1]@U[:,1], d[1]@U[:,0]+d[0]@U[:,1]])

def hypothesis_strain(eid):
    """geometri GLOBAL + displacement dirotasi ke frame (e1=n1->n2, e2=(e1y,-e1x))."""
    P = [np.array(XY[n]) for n in CONN[eid][:3]]
    e1 = P[1]-P[0]; e1 /= np.linalg.norm(e1); e2 = np.array([e1[1], -e1[0]])
    R = np.array([e1, e2])
    X = np.array([XY[n] for n in CONN[eid]])
    U = np.array([exact_disp(n, 'membrane') for n in CONN[eid]])
    return t6_strain(X, U @ R.T)

# ─────────────── laporan ───────────────
SEP = "="*112
def pct(a, b): return abs(a-b)/(abs(b)+1e-30)*100

def report_disp(name, disp):
    for sc, mode in ((1,'membrane'), (2,'bending')):
        if sc not in disp or not disp[sc]:
            print(f"  [{name}] SC{sc} {mode}: tidak ada displacement di F06"); continue
        mx = 0; worst = None
        for nid, u in disp[sc].items():
            ex = exact_disp(nid, mode)
            got = u[[0,1]] if mode == 'membrane' else u[[2,3,4]]
            den = max(np.abs(ex).max(), 1e-14)
            e = np.abs(got-ex).max()/den
            if e > mx: mx, worst = e, nid
        print(f"  [{name}] SC{sc} {mode:8s}: {len(disp[sc])} node, max rel err vs eksak = {mx:.3e} (node {worst})")

def report_elstress(name, D):
    for sc, mode in ((1,'membrane'), (2,'bending')):
        E = D['elstr'].get(sc)
        if not E: continue
        ex_g = exact_stress_global(mode)
        lab = 's11 s22 s12' if mode == 'membrane' else 'Z2:s11 s22 s12'
        print(f"\n  [{name}] SC{sc} {mode}  — stress elemen LOKAL di CENTER "
              f"({'Z1 row' if mode=='membrane' else 'Z2 row'})")
        print(f"  ex_trace={invariants(ex_g)[0]:.4f} ex_major={invariants(ex_g)[1]:.4f} ex_minor={invariants(ex_g)[2]:.4f}")
        print(f"  {'eid':>3} | {'MYSTRAN '+lab:^38} | {'exact local (e1=n1n2,e3=node-order)':^38} | "
              f"{'trace':>9} {'err%':>7} | {'major':>9} {'err%':>7} | {'minor':>9} {'err%':>7}")
        for eid in sorted(E):
            v = E[eid][0] if mode == 'membrane' else E[eid][1]
            e1, e2 = frame_n1n2(eid)
            exg_e = ex_g if mode == 'membrane' else ex_g*(-elem_normal_z(eid))
            exl = rot_local(exg_e, e1, e2)
            tr, ma, mi = invariants(v); etr, ema, emi = invariants(exg_e)
            print(f"  {eid:>3} | " + " ".join(f"{x:12.4f}" for x in v) + " | "
                  + " ".join(f"{x:12.4f}" for x in exl) + " | "
                  f"{tr:9.3f} {pct(tr,etr):7.2f} | {ma:9.3f} {pct(ma,ema):7.2f} | {mi:9.3f} {pct(mi,emi):7.2f}")

def _exact_local_row(eid, mode):
    """Stress lokal eksak (membran, atau fibre Z2 untuk bending) di frame elemen."""
    ex_g = exact_stress_global(mode)
    e1, e2 = frame_n1n2(eid)
    g = ex_g if mode == 'membrane' else ex_g*(-elem_normal_z(eid))
    return rot_local(g, e1, e2)

def _classify(v, exl, others, tol=1e-3):
    """OK | ZERO | SHIFT (= nilai eksak elemen lain) | OTHER."""
    scale = max(np.abs(exl).max(), 1e-30)
    if np.abs(v-exl).max() <= tol*scale: return 'OK', None
    if np.abs(v).max() <= 1e-6*scale:    return 'ZERO', None
    for e, c in others.items():
        if np.abs(c-v).max() <= tol*scale: return 'SHIFT', e
    return 'OTHER', None

def _others(mode):
    return {e: _exact_local_row(e, mode) for e in CONN}

def report_grid_local(name, D):
    """Tahap B: stress lokal elemen di baris GRD vs eksak. Regangan konstan -> harus sama dgn CENTER/eksak."""
    for sc, mode in ((1,'membrane'), (2,'bending')):
        G = D['elgrd'].get(sc)
        if not G: continue
        row = 0 if mode == 'membrane' else 1
        others = _others(mode)
        print(f"\n  [{name}] SC{sc} {mode} — TAHAP B: stress lokal elemen di GRD "
              f"({'Z1 row' if mode=='membrane' else 'Z2 row'}) vs eksak lokal; regangan konstan => semua grid = eksak elemen itu")
        print(f"  {'eid':>3} {'grid':>5} {'pos':>3} | {'MYSTRAN GRD':^38} | {'eksak lokal':^38} | max|d|/s | status")
        cnt = {'OK':0,'ZERO':0,'SHIFT':0,'OTHER':0}; bypos = {}
        for eid in sorted(G):
            exl = others[eid]
            for gid, vv in G[eid].items():
                v = vv[row]; tag, src = _classify(v, exl, {e: c for e, c in others.items() if e != eid})
                cnt[tag] += 1
                pos = CONN[eid].index(gid)+1 if gid in CONN[eid] else 0
                bypos.setdefault(pos, {'OK':0,'bad':0})['OK' if tag == 'OK' else 'bad'] += 1
                note = tag if tag != 'SHIFT' else f"berisi nilai elemen {src}"
                print(f"  {eid:>3} {gid:>5} {pos:>3} | " + " ".join(f"{x:12.5g}" for x in v) + " | "
                      + " ".join(f"{x:12.5g}" for x in exl)
                      + f" | {np.abs(v-exl).max()/max(np.abs(exl).max(),1e-30):.1e} | {note}")
        print(f"  --> ringkasan GRD: OK={cnt['OK']}  nol={cnt['ZERO']}  berisi elemen lain={cnt['SHIFT']}  "
              f"lain={cnt['OTHER']}  (total {sum(cnt.values())})")
        print("      per posisi grid di elemen: " + "  ".join(
            f"pos{k}: {v['OK']} OK / {v['bad']} salah" for k, v in sorted(bypos.items())))

def rot_global(v, e1, e2):
    """[s11,s22,s12] lokal (e1,e2) -> [sXX,sYY,sXY] global."""
    c11,c12,c21,c22 = e1[0],e1[1],e2[0],e2[1]
    return np.array([c11**2*v[0]+c21**2*v[1]+2*c11*c21*v[2],
                     c12**2*v[0]+c22**2*v[1]+2*c12*c22*v[2],
                     c11*c12*v[0]+c21*c22*v[1]+(c11*c22+c12*c21)*v[2]])

def elem_area(eid):
    P = [np.array(XY[n]) for n in CONN[eid][:3]]
    return 0.5*abs((P[1][0]-P[0][0])*(P[2][1]-P[0][1])-(P[1][1]-P[0][1])*(P[2][0]-P[0][0]))

def report_gp_reconstruct(name, D):
    """Tahap C': bangun ulang GP dari stress LOKAL GRD milik MYSTRAN sendiri (rotasi benar ke global,
    rata-rata berbobot luas). Jika lokal benar, hasilnya = eksak. Selisih terhadap GP MYSTRAN
    menunjukkan bug ada di rotasi/pembobotan GP, bukan di stress lokal."""
    for sc, mode in ((1,'membrane'), (2,'bending')):
        G = D['elgrd'].get(sc); P = D['gp'].get(sc)
        if not G or not P: continue
        row = 0 if mode == 'membrane' else 1; fib = 'Z1' if mode == 'membrane' else 'Z2'
        ex = exact_stress_global(mode)
        print(f"\n  [{name}] SC{sc} {mode} — TAHAP C': GP dibangun ulang dari stress LOKAL GRD MYSTRAN "
              f"(rotasi ke global benar, bobot luas elemen)")
        print(f"  {'grid':>5} {'elemen':<18} | {'rekonstruksi':^38} | {'GP MYSTRAN':^38} | max|MYSTRAN-rekon|/|eksak|")
        worst = 0.0
        for gid in CORNER_GRIDS:
            num = np.zeros(3); den = 0.0; mem = []
            for eid in sorted(G):
                if gid in G[eid]:
                    e1, e2 = frame_n1n2(eid); a = elem_area(eid)
                    num += a*rot_global(G[eid][gid][row], e1, e2); den += a; mem.append(eid)
            if not den or gid not in P or fib not in P[gid]: continue
            rec = num/den; my = P[gid][fib]; d = np.abs(my-rec).max()/np.abs(ex).max(); worst = max(worst, d)
            print(f"  {gid:>5} {str(mem):<18} | " + " ".join(f"{x:12.5g}" for x in rec) + " | "
                  + " ".join(f"{x:12.5g}" for x in my) + f" | {d:.2e}")
        print(f"  --> {'GP KONSISTEN dengan stress lokal' if worst < 1e-3 else 'GP TIDAK konsisten dengan stress lokal (bug di rotasi/bobot GP)'}"
              f" (maks {worst:.2e})")

def report_forces(name, D):
    """Engineering forces/moments lokal: M = Z2_eksak * h^2/6 (bending)."""
    F = D['elfrc'].get(2)
    if not F: return
    print(f"\n  [{name}] SC2 bending — element engineering moments (lokal) vs eksak (M = Z2*h^2/6)")
    print(f"  {'eid':>3} | {'MYSTRAN CENTER Mxx Myy Mxy':^42} | {'eksak':^42} | max|d|/|M|  | GRD rows")
    for eid in sorted(F):
        c = F[eid].get('CENTER')
        if c is None: continue
        exm = _exact_local_row(eid, 'bending')*H_**2/6.0
        err = np.abs(c[3:6]-exm).max()/np.abs(exm).max()
        g = [k for k in F[eid] if k != 'CENTER']
        gz = sum(1 for k in g if np.abs(F[eid][k]).max() < 1e-30)
        print(f"  {eid:>3} | " + " ".join(f"{x:13.5E}" for x in c[3:6]) + " | "
              + " ".join(f"{x:13.5E}" for x in exm) + f" | {err:10.2e} | {len(g)} baris, {gz} nol")

def report_summary(name, D):
    """Ringkasan tahap A -> B -> C."""
    print(f"\n  [{name}] RINGKASAN TAHAP (A: CENTER lokal -> B: GRD lokal -> C: GP global)")
    for sc, mode in ((1,'membrane'), (2,'bending')):
        C = D['elstr'].get(sc); G = D['elgrd'].get(sc); P = D['gp'].get(sc)
        row = 0 if mode == 'membrane' else 1
        a = b = c = 'n/a'
        if C:
            worst = max(np.abs(v[row] - _exact_local_row(e, mode)).max()/np.abs(_exact_local_row(e, mode)).max() for e, v in C.items())
            a = f"{'OK' if worst < 1e-3 else 'BEDA'} (maks rel {worst:.1e})"
        elif G:
            a = "baris CENTER tidak dicetak"
        if G:
            oth = _others(mode); bad = tot = 0
            for e in G:
                for gid, vv in G[e].items():
                    tot += 1; bad += (_classify(vv[row], oth[e], {k: x for k, x in oth.items() if k != e})[0] != 'OK')
            b = f"{'OK' if bad == 0 else 'BEDA'} ({tot-bad}/{tot} grid benar)"
        if P:
            fib = 'Z1' if mode == 'membrane' else 'Z2'; ex = exact_stress_global(mode)
            r = [(v[fib][0]+v[fib][1])/(ex[0]+ex[1]) for g, v in P.items() if fib in v and g in CORNER_GRIDS]
            dev = max(np.abs(v[fib]-ex).max()/np.abs(ex).max() for g, v in P.items() if fib in v and g in CORNER_GRIDS)
            c = (f"trace/eksak {min(r):.4f}..{max(r):.4f} ({'OK' if max(abs(x-1) for x in r) < 1e-3 else 'BEDA'}); "
                 f"komponen maks rel err {dev*100:.1f}% ({'OK' if dev < 1e-3 else 'BEDA'})")
        print(f"    SC{sc} {mode:8s}: A={a} | B={b} | C={c}")

def report_gp(name, D):
    for sc, mode in ((1,'membrane'), (2,'bending')):
        G = D['gp'].get(sc)
        if not G: continue
        ex = exact_stress_global(mode)
        fib = 'Z1' if mode == 'membrane' else 'Z2'
        print(f"\n  [{name}] SC{sc} {mode}  — GP stress (surface X-Y global), fiber {fib}; "
              f"eksak = {np.round(ex,4)}")
        print(f"  {'grid':>5} | {'sXX':>12} {'sYY':>12} {'sXY':>12} | {'err% XX':>8} {'err% YY':>8} {'err% XY':>8} | trace/eksak")
        for gid in CORNER_GRIDS:
            if gid not in G or fib not in G[gid]: continue
            v = G[gid][fib]
            print(f"  {gid:>5} | " + " ".join(f"{x:12.4f}" for x in v) + " | "
                  + " ".join(f"{pct(v[k],ex[k]):8.2f}" for k in range(3))
                  + f" | {(v[0]+v[1])/(ex[0]+ex[1]):9.4f}")

def exact_local_strain(eid):
    e1, e2 = frame_n1n2(eid)
    T = np.array([[STRAIN_REF[0], STRAIN_REF[2]/2],[STRAIN_REF[2]/2, STRAIN_REF[1]]])
    R = np.array([e1[:2], e2[:2]]); Tl = R @ T @ R.T
    return np.array([Tl[0,0], Tl[1,1], 2*Tl[0,1]])

def report_recov(name, D):
    R = D['recov'].get(1)
    if not R:
        print(f"  [{name}] tidak ada blok RECOV249 (DEBUG 249 tidak aktif)"); return
    print(f"\n  [{name}] RECOV249 STRAIN1 (blok ke-1, stress elemen SC1) vs regangan EKSAK lokal")
    print(f"  dan vs pola bug lama (geometri GLOBAL + displacement LOKAL); strain eksak global = {STRAIN_REF}")
    print(f"  {'eid':>3} | {'RECOV STRAIN1 (pt1)':^46} | {'eksak lokal':^46} | diff-eksak | diff-bug-lama")
    worst = 0; worst_old = 0
    for eid in sorted(CONN):
        blk = R.get((eid, 1))
        if not blk or 'STRAIN1' not in blk: continue
        s2 = blk['STRAIN1']; ex = exact_local_strain(eid); h = hypothesis_strain(eid)
        d = np.abs(s2-ex).max(); d_old = np.abs(s2-h).max(); worst = max(worst, d); worst_old = max(worst_old, d_old)
        print(f"  {eid:>3} | " + " ".join(f"{x:14.6E}" for x in s2) + " | "
              + " ".join(f"{x:14.6E}" for x in ex) + f" | {d:.2e}   | {d_old:.2e}")
    print(f"  --> maks selisih vs EKSAK: {worst:.2e} ({'COCOK' if worst < 1e-8 else 'TIDAK cocok'}); "
          f"vs bug lama: {worst_old:.2e} ({'masih pola bug lama' if worst_old < 1e-8 else 'bukan bug lama'})")

def report_python(pyname, D, name):
    try:
        sys.path.insert(0, '.')
        import print_pernode_stress_2001_v3 as v3
    except Exception as ex:
        print(f"\n  [python] gagal import v3: {ex}"); return
    from Simo1993_Tri6_ShellElement_v2 import Simo1993_Tri6_ShellElement_v2
    from MITC6_Tri_v4 import MITC6_Tri_v4
    from MacNeal_MH6T_Tri_v3 import MacNeal_MH6T_Tri_v3
    from Rezaiee2017_Tri6_v3 import Rezaiee2017_Tri6_ShellElement_v3
    final_t6 = {
        'SIMOT6': Simo1993_Tri6_ShellElement_v2,
        'MITC6': MITC6_Tri_v4,
        'MH6T': MacNeal_MH6T_Tri_v3,
        'REZAIEE': Rezaiee2017_Tri6_ShellElement_v3,
    }
    aliases = {'SimoT6': 'SIMOT6', 'SimoT6v2': 'SIMOT6', 'MITC6v4': 'MITC6',
               'MHT6v3': 'MH6T', 'Rezaiee_v3': 'REZAIEE'}
    pyname = aliases.get(pyname, pyname)
    v3.FLIP_SHEAR_IF_NEG_E3 = False
    if pyname in final_t6:
        from dump_python_debug import resolve_reference
        pyname, entry = resolve_reference(v3, pyname)
    elif pyname in v3.REGISTRY:
        entry = v3.REGISTRY[pyname]
    else:
        print(f"\n  [python] '{pyname}' tidak ada di REGISTRY: {list(v3.REGISTRY)}"); return
    cls, nxy, conn, bnd, x0, y0, tag = entry
    for sc, mode in ((1,'membrane'), (2,'bending')):
        Cm = D['elstr'].get(sc) or {}
        G = D['elgrd'].get(sc, {})
        if not Cm and not G: continue
        row = 0 if mode == 'membrane' else 1
        fem, elems = v3.build_and_solve(cls, nxy, conn, bnd, x0, y0, mode)
        print(f"\n  [PYTHON {pyname} vs MYSTRAN {name}] SC{sc} {mode} — TAHAP A+B: stress LOKAL per titik "
              f"({'membran' if mode=='membrane' else 'fibre Z2'}); frame eksak: e1=n1n2, e3 dari urutan node")
        print(f"  {'eid':>3} {'pt':>4} {'grid':>5} | {'python s11 s22 s12':^38} | {'eksak lokal':^38} | "
              f"{'MYSTRAN':^38} | d(py-ex) | d(py-my)")
        acc = {}; nbad_ex = nbad_my = nmy = 0; npts = 0
        for eid in sorted(elems):
            pts = v3.collect_element(elems[eid], fem, conn[eid], mode, tag)
            exl = _exact_local_row(eid, mode)
            for p in pts:
                py = p['sig_loc'] if mode == 'membrane' else v3.fibers(p['sig_loc'], p['mom_loc'])[1]
                gid = p['nid']
                if gid is None: my = Cm[eid][row] if eid in Cm else None
                else:           my = G.get(eid, {}).get(gid, (None, None))[row]
                scale = max(np.abs(exl).max(), 1e-30)
                d_ex = np.abs(py-exl).max()/scale
                d_my = (np.abs(py-my).max()/scale) if my is not None else None
                npts += 1; nbad_ex += d_ex > 1e-3
                if d_my is not None: nmy += 1; nbad_my += d_my > 1e-3
                mys = " ".join(f"{x:12.5g}" for x in my) if my is not None else f"{'-':^38}"
                print(f"  {eid:>3} {p['lbl']:>4} {str(gid or '---'):>5} | " + " ".join(f"{x:12.5g}" for x in py)
                      + " | " + " ".join(f"{x:12.5g}" for x in exl) + " | " + mys
                      + f" | {d_ex:8.1e} | " + (f"{d_my:8.1e}" if d_my is not None else "     -  "))
                if gid is not None and gid in CORNER_GRIDS:
                    key = 'sig_glb' if mode == 'membrane' else 'mom_glb'
                    g = p[key] if mode == 'membrane' else v3.fibers(p['sig_glb'], p['mom_glb'])[1]
                    acc.setdefault(gid, []).append(g)
        print(f"  --> Python vs eksak lokal: {npts-nbad_ex}/{npts} titik cocok (rel<1e-3); "
              f"Python vs MYSTRAN: {nmy-nbad_my}/{nmy} titik cocok")
        if mode == 'bending' and D['elfrc'].get(2):
            print(f"  Momen lokal CENTER Python vs MYSTRAN (engineering forces):")
            for eid in sorted(elems):
                p0 = v3.collect_element(elems[eid], fem, conn[eid], mode, tag)[0]
                f = D['elfrc'][2].get(eid, {}).get('CENTER')
                if f is not None:
                    print(f"  {eid:>3} | py " + " ".join(f"{x:13.5E}" for x in p0['mom_loc'])
                          + " | my " + " ".join(f"{x:13.5E}" for x in f[3:6]))
        Gp = D['gp'].get(sc, {}); fib = 'Z1' if mode == 'membrane' else 'Z2'
        print(f"  TAHAP C — GP-average Python (dari stress global per elemen-grid) vs GP MYSTRAN ({fib}) vs eksak:")
        ex = exact_stress_global(mode)
        for gid in CORNER_GRIDS:
            if gid in acc:
                pv = np.mean(acc[gid], axis=0)
                line = f"  {gid:>5} | py " + " ".join(f"{x:12.4f}" for x in pv)
                if gid in Gp and fib in Gp[gid]:
                    line += " | my " + " ".join(f"{x:12.4f}" for x in Gp[gid][fib])
                print(line + f" | trace py/eksak {(pv[0]+pv[1])/(ex[0]+ex[1]):.4f}")

def main(argv):
    global NORMAL_MODE
    files = []; nastran = None; pyname = None
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == '--nastran': nastran = argv[i+1]; i += 2
        elif a == '--python': pyname = argv[i+1]; i += 2
        elif a == '--normal': NORMAL_MODE = argv[i+1].lower(); i += 2
        else: files.append(a); i += 1
    if not files and not nastran:
        print(__doc__); return
    print(SEP); print("  F06 vs EKSAK  (patch test 2-001, T6)")
    ex = exact_stress_global('membrane'); exb = exact_stress_global('bending')
    print(f"  eksak membran  : sXX={ex[0]:.4f} sYY={ex[1]:.4f} sXY={ex[2]:.4f}   "
          f"trace={invariants(ex)[0]:.4f} major={invariants(ex)[1]:.4f} minor={invariants(ex)[2]:.4f}")
    print(f"  eksak bending Z2 (+h/2): {np.round(exb,5)}   (Z1 = tanda kebalikan)")
    print(SEP)
    if nastran:
        print("\n### NASTRAN displacement"); report_disp('NASTRAN', parse_nastran_disp(nastran))
    for f in files:
        name = re.sub(r'\.f06$', '', f.split('/')[-1], flags=re.I)
        D = parse_mystran(f)
        print(f"\n{SEP}\n### MYSTRAN: {f}\n{SEP}")
        report_disp(name, D['disp'])
        report_elstress(name, D)
        report_grid_local(name, D)
        report_forces(name, D)
        report_gp(name, D)
        report_gp_reconstruct(name, D)
        report_recov(name, D)
        if pyname: report_python(pyname, D, name)
        report_summary(name, D)
    print(f"\n{SEP}\n  Selesai.\n{SEP}")

if __name__ == '__main__':
    main(sys.argv[1:])

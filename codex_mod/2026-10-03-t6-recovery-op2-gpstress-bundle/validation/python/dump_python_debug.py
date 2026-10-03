"""
dump_python_debug.py
====================
Dump semua matriks/variabel pembentuk stress elemen Python ke file teks terpisah
(satu file per elemen-kelas per mode), dengan format baris seragam supaya mudah
di-diff dengan output debug Fortran MYSTRAN.

Format baris (satu nilai-array per baris):
    LABEL: v1 v2 v3 ...
Contoh Fortran yang menghasilkan format yang sama:
    WRITE(UNIT,'(A,":",100(1X,ES19.10E3))') TRIM(LABEL), (VAL(I), I=1,N)
(nilai dibaca sebagai float, jadi jumlah digit eksponen tidak masalah)

Label bertingkat:  E<eid>/<bagian>[/<titik>]
    E51/CONN                        id grid elemen
    E51/XY_GLOBAL/n1                x y z node 1
    E51/T_MATRIX/r01                baris 1 matriks transformasi global->lokal
    E51/U_GLOBAL                    displacement elemen (6 per node, urutan global)
    E51/U_LOCAL                     T @ U_GLOBAL
    E51/FRAME_CTR                   e1(3) e2(3) e3(3) di pusat
    E51/CTR/BM/r1 ...               baris ke-1 matriks Bm di titik CTR
    E51/CTR/STRAIN_M                Bm @ u_e
    E51/n1/STRESS_LOC               stress lokal (frame pusat), dst.

Usage:
    # Active T6: SIMOT6, MITC6, MH6T, REZAIEE (final Python references).
    python dump_python_debug.py MITC6                 # membran + bending -> debug_MITC6_*.txt
    python dump_python_debug.py MITC6 --mode membrane --elem 51 53
    python dump_python_debug.py MITC6 --attrs         # + semua array numpy atribut elemen
    python dump_python_debug.py --diff a.txt b.txt [--tol 1e-6]   # bandingkan dua dump
"""
import sys
import numpy as np


def resolve_reference(v3, name):
    """Resolve active T6 names to final references; retain Q8 registry entries."""
    from Simo1993_Tri6_ShellElement_v2 import Simo1993_Tri6_ShellElement_v2
    from MITC6_Tri_v4 import MITC6_Tri_v4
    from MacNeal_MH6T_Tri_v3 import MacNeal_MH6T_Tri_v3
    from Rezaiee2017_Tri6_v3 import Rezaiee2017_Tri6_ShellElement_v3
    references = {
        'SIMOT6': Simo1993_Tri6_ShellElement_v2,
        'MITC6': MITC6_Tri_v4,
        'MH6T': MacNeal_MH6T_Tri_v3,
        'REZAIEE': Rezaiee2017_Tri6_ShellElement_v3,
    }
    aliases = {'SimoT6': 'SIMOT6', 'SimoT6v2': 'SIMOT6',
               'MITC6v4': 'MITC6', 'MHT6v3': 'MH6T', 'MHT6': 'MH6T',
               'Rezaiee_v3': 'REZAIEE'}
    name = aliases.get(name, name.upper() if name.upper() in references else name)
    if name in references:
        # Final classes use L1=1-r-s, L2=r, L3=s (the legacy helper swaps nodes).
        v3.T6_NODE_RS = [(0., 0.), (1., 0.), (0., 1.),
                         (0.5, 0.), (0.5, 0.5), (0., 0.5)]
        v3.t6_N = lambda r, s: np.array([
            (1-r-s)*(1-2*r-2*s), r*(2*r-1), s*(2*s-1),
            4*(1-r-s)*r, 4*r*s, 4*s*(1-r-s)])
        return name, v3._t6(references[name])
    return name, v3.REGISTRY.get(name)


# ═══════════════ penulisan ═══════════════
def w(fh, label, vals):
    a = np.atleast_1d(np.asarray(vals, dtype=float)).ravel()
    fh.write(f"{label}: " + " ".join(f"{x:.10E}" for x in a) + "\n")

def w_mat(fh, label, M):
    M = np.atleast_2d(np.asarray(M, dtype=float))
    for i, row in enumerate(M, start=1):
        w(fh, f"{label}/r{i:02d}", row)


# ═══════════════ dump satu elemen ═══════════════
def dump_element(fh, v3, e, fem, conn, node_xy, tag, mode, with_attrs=False, eid=None):
    """v3 = modul print_pernode_stress_2001_v3 (dipakai: konstanta, basis, Bm/Bb)."""
    if eid is None:
        eid = e.id if hasattr(e, 'id') else '?'
    is_q8 = (tag == 'Q8')
    node_rs  = v3.Q8_NODE_RS  if is_q8 else v3.T6_NODE_RS
    node_lbl = v3.Q8_NODE_LBL if is_q8 else v3.T6_NODE_LBL
    ctr      = v3.Q8_CENTER   if is_q8 else v3.T6_CENTER
    Nfun     = v3.q8_N        if is_q8 else v3.t6_N

    fh.write(f"\n# ================= ELEM {eid}  TYPE {tag}  MODE {mode} =================\n")
    w(fh, f"E{eid}/CONN", conn)
    for k, nid in enumerate(conn, start=1):
        x, y = node_xy[nid]
        w(fh, f"E{eid}/XY_GLOBAL/n{k}", [x, y, 0.0])

    u_global = fem.u
    idx = np.asarray(e.global_dof_indices())
    T = np.asarray(e.T_matrix())
    u_g = u_global[idx]
    u_e = T @ u_g
    w(fh, f"E{eid}/GLOBAL_DOF_INDEX", idx)
    w_mat(fh, f"E{eid}/T_MATRIX", T)
    w(fh, f"E{eid}/U_GLOBAL", u_g)
    w(fh, f"E{eid}/U_LOCAL", u_e)
    w_mat(fh, f"E{eid}/U_LOCAL_PERNODE", u_e.reshape(-1, 6))

    # material (padanan EM / EB / ET / ZS di RECOV249 MYSTRAN)
    C = v3.C_
    w_mat(fh, f"E{eid}/CMAT", C)
    w(fh, f"E{eid}/EM_row1", C[0, :])
    w(fh, f"E{eid}/ET_G", [C[2, 2], 0.0])
    w(fh, f"E{eid}/ZS", [-v3.H_/2, v3.H_/2])

    # frame di pusat (frame lokal elemen untuk output)
    c1, c2, c3 = v3.get_basis(e, ctr[0], ctr[1])
    w(fh, f"E{eid}/FRAME_CTR", np.concatenate([c1, c2, c3]))

    pts = [('CTR', ctr)] + list(zip(node_lbl, node_rs))
    for lbl, (r, s) in pts:
        base = f"E{eid}/{lbl}"
        w(fh, f"{base}/RS", [r, s])
        w(fh, f"{base}/N", Nfun(r, s))
        e1, e2, e3 = v3.get_basis(e, r, s)
        w(fh, f"{base}/FRAME_PT", np.concatenate([e1, e2, e3]))
        Bm, Bb = v3.get_Bm_Bb(e, r, s)
        Bm = np.asarray(Bm); Bb = np.asarray(Bb)
        w_mat(fh, f"{base}/BM", Bm)
        w_mat(fh, f"{base}/BB", Bb)
        if hasattr(e, '_compute_Bs'):
            try:
                Bs = np.asarray(v3.get_Bs(e, r, s)) if hasattr(v3, 'get_Bs') else np.asarray(e._compute_Bs(r, s))
                w_mat(fh, f"{base}/BS", Bs)
                w(fh, f"{base}/STRAIN_S", Bs @ u_e)
            except Exception as ex:
                fh.write(f"# {base}/BS gagal: {ex}\n")
        eps = Bm @ u_e                       # regangan membran (frame titik)
        kap = Bb @ u_e                       # kurvatur
        sig = C @ eps
        mom = C @ kap * v3.H_**3/12.
        w(fh, f"{base}/STRAIN_M", eps)
        w(fh, f"{base}/KAPPA", kap)
        w(fh, f"{base}/STRESS_M_PT", sig)
        w(fh, f"{base}/MOMENT_PT", mom)
        # padanan RECOV249: STRAIN1 = regangan membran; serat = eps -/+ z*kappa
        w(fh, f"{base}/STRAIN_Z1", eps - (v3.H_/2)*kap)
        w(fh, f"{base}/STRAIN_Z2", eps + (v3.H_/2)*kap)
        res = v3.eval_point(e, u_global, r, s, ctr)
        w(fh, f"{base}/STRESS_LOC", res['sig_loc'])
        w(fh, f"{base}/MOMENT_LOC", res['mom_loc'])
        w(fh, f"{base}/STRESS_GLB", res['sig_glb'])
        w(fh, f"{base}/MOMENT_GLB", res['mom_glb'])
        z1, z2 = v3.fibers(res['sig_loc'], res['mom_loc'])
        w(fh, f"{base}/FIBER_Z1_LOC", z1)
        w(fh, f"{base}/FIBER_Z2_LOC", z2)

    if with_attrs:
        for name, val in sorted(vars(e).items()):
            try:
                a = np.asarray(val, dtype=float)
            except (TypeError, ValueError):
                continue
            if a.size == 0 or a.size > 4000:
                continue
            if a.ndim <= 1:
                w(fh, f"E{eid}/ATTR/{name}", a)
            elif a.ndim == 2:
                w_mat(fh, f"E{eid}/ATTR/{name}", a)


# ═══════════════ diff dua dump ═══════════════
def load_dump(path):
    d = {}
    for line in open(path, errors='ignore'):
        if not line.strip() or line.startswith('#') or ':' not in line:
            continue
        lab, rest = line.split(':', 1)
        try:
            d[lab.strip()] = np.array([float(x) for x in rest.split()])
        except ValueError:
            pass
    return d

def diff_dumps(pa, pb, tol=1e-6, maxshow=40):
    A, B = load_dump(pa), load_dump(pb)
    common = [k for k in A if k in B]
    only_a = [k for k in A if k not in B]; only_b = [k for k in B if k not in A]
    bad = []
    for k in common:
        a, b = A[k], B[k]
        if a.shape != b.shape:
            bad.append((k, np.inf, f"ukuran beda {a.shape} vs {b.shape}")); continue
        den = max(np.abs(a).max(), np.abs(b).max(), 1e-30)
        err = np.abs(a-b).max()
        if err > tol*den and err > 1e-14:
            bad.append((k, err/den, f"maxabs={err:.3e}"))
    print(f"label sama: {len(common)} | hanya di A: {len(only_a)} | hanya di B: {len(only_b)} "
          f"| beda (tol rel {tol:g}): {len(bad)}")
    print("Urutan pertama yang beda (mengikuti urutan file A -> tracing dari hulu ke hilir):")
    for k, e, msg in bad[:maxshow]:
        print(f"  {k:<40s} rel={e:.3e}  {msg}")
    if only_a[:5]: print("  contoh hanya di A:", only_a[:5])
    if only_b[:5]: print("  contoh hanya di B:", only_b[:5])
    return len(bad)


# ═══════════════ main ═══════════════
def main(argv):
    if '--diff' in argv:
        i = argv.index('--diff'); tol = 1e-6
        if '--tol' in argv: tol = float(argv[argv.index('--tol')+1])
        diff_dumps(argv[i+1], argv[i+2], tol); return
    if not argv or argv[0].startswith('--'):
        print(__doc__); return
    pyname = argv[0]
    modes = ['membrane', 'bending']
    if '--mode' in argv: modes = [argv[argv.index('--mode')+1]]
    only_elems = None
    if '--elem' in argv:
        j = argv.index('--elem') + 1; only_elems = []
        while j < len(argv) and not argv[j].startswith('--'):
            only_elems.append(int(argv[j])); j += 1
    with_attrs = '--attrs' in argv

    sys.path.insert(0, '.')
    import print_pernode_stress_2001_v3 as v3
    pyname, entry = resolve_reference(v3, pyname)
    if entry is None:
        print(f"'{pyname}' tidak ada. Pilihan: SIMOT6, MITC6, MH6T, REZAIEE, {list(v3.REGISTRY)}"); return
    cls, nxy, conn, bnd, x0, y0, tag = entry
    if tag == 'T6':
        v3.FLIP_SHEAR_IF_NEG_E3 = False
    for mode in modes:
        fem, elems = v3.build_and_solve(cls, nxy, conn, bnd, x0, y0, mode)
        out = f"debug_{pyname}_{mode}.txt"
        with open(out, 'w') as fh:
            fh.write(f"# dump_python_debug: {pyname} {mode}\n")
            fh.write(f"# reference: {cls.__module__}.{cls.__name__}\n")
            fh.write(f"# E={v3.E_} NU={v3.NU_} H={v3.H_}\n")
            for eid in sorted(elems):
                if only_elems and eid not in only_elems: continue
                dump_element(fh, v3, elems[eid], fem, conn[eid], nxy, tag, mode, with_attrs, eid=eid)
        print(f"ditulis: {out}")

if __name__ == '__main__':
    main(sys.argv[1:])

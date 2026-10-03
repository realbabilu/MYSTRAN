"""Strict A/B/C coverage and stress comparison using compare_f06_vs_python."""
from pathlib import Path
import sys, os, subprocess, json, contextlib
import numpy as np
root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root/'python'))
import compare_f06_vs_python as compare
import print_pernode_stress_2001_v3 as v3
from dump_python_debug import resolve_reference
folder = root/'test/t6_upgrade/stress_recovery'
folder.mkdir(parents=True, exist_ok=True)
env = os.environ.copy()
env['PATH'] = str(root/'test')+';C:\\gcc\\bin;'+env['PATH']
results = []
for family in ('SIMOT6', 'MITC6', 'MH6T', 'REZAIEE'):
    deck = folder/f'{family}.dat'
    source = (root/'test/duel3a_mitc6_1.dat').read_text()
    source = source.replace('PARAM,TRIA6TYP,MITC6', f'PARAM,TRIA6TYP,{family}')
    source = source.replace('DEBUG,249,1', '$DEBUG,249,1')
    source = source.replace('  FORCES(CENTER,CORNER)= ALL',
                            '  FORCES(CENTER,CORNER)= ALL\n  STRESS(CENTER,CORNER)= ALL')
    deck.write_text(source)
    run = subprocess.run([str(root/'MYSTRAN/Binaries/mystran.exe'),deck.name],
                         cwd=folder, env=env, capture_output=True, text=True, timeout=90)
    deck.with_suffix('.stdout.txt').write_text(run.stdout+run.stderr)
    D = compare.parse_mystran(str(deck.with_suffix('.F06')))
    with deck.with_suffix('.compare.txt').open('w',encoding='utf-8') as out, contextlib.redirect_stdout(out):
        compare.main([str(deck.with_suffix('.F06')), '--python', family])
    _, entry = resolve_reference(v3, family)
    cls, nxy, conn, bnd, x0, y0, tag = entry
    v3.FLIP_SHEAR_IF_NEG_E3 = False
    for sc, mode in ((1,'membrane'),(2,'bending')):
        fem, elems = v3.build_and_solve(cls, nxy, conn, bnd, x0, y0, mode)
        errors = {'A':[], 'B':[], 'C':[], 'force':[]}
        missing = []
        acc = {}
        def check(stage, label, actual, expected):
            if actual is None:
                missing.append([stage,label]); return
            expected = np.asarray(expected)
            errors[stage].append(float(np.max(np.abs(actual-expected))/max(np.max(np.abs(expected)),1e-8)))
        for eid, e in elems.items():
            for p in v3.collect_element(e,fem,conn[eid],mode,tag):
                gid = p['nid']
                actual = D['elstr'].get(sc,{}).get(eid) if gid is None else D['elgrd'].get(sc,{}).get(eid,{}).get(gid)
                check('A' if gid is None else 'B', [eid,gid], actual,
                      v3.fibers(p['sig_loc'],p['mom_loc']))
                if gid is not None:
                    acc.setdefault(gid,[]).append(v3.fibers(p['sig_glb'],p['mom_glb']))
                if sc == 2:
                    force = D['elfrc'].get(sc,{}).get(eid,{}).get('CENTER' if gid is None else gid)
                    check('force',[eid,gid], None if force is None else force[3:6],p['mom_loc'])
        for gid, vals in acc.items():
            gp = D['gp'].get(sc,{}).get(gid,{})
            check('C',gid,None if 'Z1' not in gp or 'Z2' not in gp else [gp['Z1'],gp['Z2']],np.mean(vals,axis=0))
        record = dict(family=family,mode=mode,exit_code=run.returncode,
                      counts={k:len(v) for k,v in errors.items()},
                      max_rel={k:max(v,default=0) for k,v in errors.items()},missing=missing)
        record['pass_gate'] = not missing and all(x < 1e-4 for x in record['max_rel'].values())
        print(json.dumps(record),flush=True)
        results.append(record)
        (folder/'results.json').write_text(json.dumps(results,indent=2))
print('PASS',sum(r['pass_gate'] for r in results),'/',len(results))

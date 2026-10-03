"""T6-only displacement gate: identical Python models exported to MYSTRAN."""
from pathlib import Path
import sys, re, json, subprocess, argparse, os
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'python'))
import test_patch2001_q8_t6_v4c as patch
import battle_2_002_quadratic_clamped as clamped
import battle_2_003_quadratic_curved as curved
import battle_2_004_twistedv2 as twisted
from pyNastran.bdf.field_writer_16 import print_card_16
import compare_f06_vs_python as compare

def model_for(problem, family, nx, lc):
    cls = clamped.T6_ELEMENTS[family]
    if problem == '2-001':
        cls, sfn, tag, xy, conn, bnd, x0, y0 = patch.REGISTRY[family]
        *_, model = patch.build_and_solve(cls, xy, conn, bnd, x0, y0,
            'membrane' if lc == 1 else 'bending', sfn, tag=tag)
        return model
    if problem == '2-002':
        model, fb, ft, tb, tt = clamped.build_t6(nx, cls)
        clamped.apply_load_case(model, fb, ft, tb, tt, lc)
    elif problem == '2-003':
        model, fi, fo, ti, to = curved.build_t6(nx, cls, True)
        curved.apply_bc_load(model, fi, fo, ti, to, lc)
    else:
        model, fixed, tips = twisted.build_t6(nx, cls, alternate_diag=True)
        twisted.apply_bc_load(model, fixed, tips, lc)
    model.build()
    model.solve_static()
    return model

def export(model, family, path):
    e = model.elements[0]
    cards = ['SOL 101\nCEND\nECHO = NONE\nSUBCASE 1\n SPC = 1\n LOAD = 2\n DISPLACEMENT = ALL\nBEGIN BULK\n',
             print_card_16(['PARAM', 'TRIA6TYP', family]),
             print_card_16(['PARAM', 'POST', -1]),
             print_card_16(['PARAM', 'AUTOSPC', 'Y' if model.autospc else 'N']),
             print_card_16(['MAT1', 1, float(e.E), None, float(e.nu)]),
             print_card_16(['PSHELL', 1, 1, float(np.mean(e.h)), 1, None, 1])]
    for nid, node in model.nodes.items():
        cards.append(print_card_16(['GRID', nid, None, float(node.x), float(node.y), float(node.z)]))
    for i, element in enumerate(model.elements, 1):
        cards.append(print_card_16(['CTRIA6', i, 1, *[n.nid for n in element.nodes]]))
    for bc in model.bcs:
        cards.append(print_card_16(['SPC', 1, bc.node_id, bc.dof_local + 1, 0.0]))
        if bc.value:
            cards.append(print_card_16(['SPCD', 2, bc.node_id, bc.dof_local + 1, float(bc.value)]))
    for load in model.loads:
        vec = [0., 0., 0.]
        vec[load.dof_local % 3] = 1.
        cards.append(print_card_16(['FORCE' if load.dof_local < 3 else 'MOMENT',
            2, load.node_id, None, float(load.value), *vec]))
    cards.append('ENDDATA\n')
    path.write_text(''.join(cards), encoding='ascii')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--phase', default='baseline')
    ap.add_argument('--mesh', type=int, nargs='+', default=[2, 4])
    ap.add_argument('--families', nargs='+', default=['SIMOT6','MITC6','MH6T','REZAIEE'])
    ap.add_argument('--problems', nargs='+', default=['2-001','2-002','2-003','2-004'])
    args = ap.parse_args()
    folder = ROOT / 'test' / 't6_upgrade' / args.phase
    folder.mkdir(parents=True, exist_ok=True)
    results = []
    for problem in args.problems:
        for nx in ([0] if problem == '2-001' else args.mesh):
            for family in args.families:
                for lc in range(1, 7 if problem == '2-002' else 3):
                    model = model_for(problem, family, nx, lc)
                    deck = folder / f'{problem}_{family}_n{nx}_lc{lc}.dat'
                    export(model, family, deck)
                    env = os.environ.copy()
                    env['PATH'] = str(ROOT/'test') + ';C:\\gcc\\bin;' + env['PATH']
                    run = subprocess.run([str(ROOT/'MYSTRAN/Binaries/mystran.exe'), deck.name],
                        cwd=folder, env=env, capture_output=True, text=True, timeout=90)
                    deck.with_suffix('.stdout.txt').write_text(run.stdout+run.stderr, encoding='utf-8')
                    f06 = deck.with_suffix('.F06')
                    record = dict(problem=problem, nx=nx, family=family, lc=lc, exit_code=run.returncode)
                    try:
                        disp = compare.parse_mystran(str(f06))['disp'][1]
                        expected = np.array([model.get_displacement(nid) for nid in sorted(model.nodes)])
                        actual = np.array([disp[nid] for nid in sorted(model.nodes)])
                        scale = np.max(np.abs(expected), axis=0)
                        err = np.max(np.abs(actual-expected),axis=0)
                        relative = err / np.maximum(scale, 1e-10)
                        record.update(pass_gate=bool(np.all(err <= 1e-10 + 1e-4*scale)),
                            max_rel=float(np.max(relative)), abs_by_dof=err.tolist(),
                            scale_by_dof=scale.tolist(), nodes=len(disp))
                    except Exception as exc:
                        record.update(pass_gate=False, error=str(exc))
                    results.append(record)
                    print(json.dumps(record), flush=True)
                    (folder/'results.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    print('PASS',sum(r['pass_gate'] for r in results),'/',len(results),flush=True)

if __name__ == '__main__':
    main()

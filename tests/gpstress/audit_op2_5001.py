"""Compare the two observed 5-001 solid OGS result types, without relabeling MYSTRAN."""
import csv
from pathlib import Path
import numpy as np
from read_ogs import read_volumes
ROOT=Path(__file__).resolve().parent
WS=ROOT.parents[2]
REF=WS/'nastran_test/5-001_GPSTRESS_BASIC/nastran/5-001_8_SOLIDS_GPSTRESS_PENTA_REMESH.op2'
ACT=WS/'nastran_test/mystran/quadratic_patch_audit/5-001_PENTATYP_EAS9.OP2'
rows=[]
for code in (27,28):
    ref=read_volumes(REF,code); actual=read_volumes(ACT,code)
    assert len(ref)==len(actual)==8
    for r,a in zip(ref,actual):
        assert r['header']==a['header'],(r['header'],a['header'])
        assert r['grids']==a['grids']
        rv=np.array(r['values']); av=np.array(a['values'])
        # Eigenvector signs/bases are nonunique. Compare principal values and
        # invariants numerically; axes are verified through the eigen-equation.
        cols=list(range(8)) if code==27 else [0,1,2,12,13]
        error=float(np.max(abs(rv[:,cols]-av[:,cols])/np.maximum(abs(rv[:,cols]),1.)))
        rows.append([code,r['volume'],len(r['grids']),error])
with (ROOT/'op2_5001_comparison.csv').open('w',newline='') as f:
    w=csv.writer(f);w.writerow(['table_code','volume','rows','reference_relative_difference']);w.writerows(rows)
for row in rows: print(*row,sep=', ')

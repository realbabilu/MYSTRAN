from pathlib import Path
import json,sys,hashlib,contextlib
import numpy as np
P=Path(__file__).resolve().parent;sys.path.insert(0,'H:/Projects/MYSTRAN/shell_linear_benchmark_20261023');from run_suite import run;from results import parse_f06
np.in1d=np.isin
from numpy.char import chararray
np.chararray=chararray
from pyNastran.op2.op2 import read_op2
rows=json.loads((P/'first_run_results.json').read_text());out=[]
for r in rows:
 if r['tag']!='restrained_gradient' or r['version']!='after':continue
 text=Path(r['deck']).read_text().replace('BEGIN BULK','BEGIN BULK\nPARAM,POST,-1\nPARAM,POSTEXT,YES');d=P/'op2_check'/r['formulation']/'model.dat';d.parent.mkdir(parents=True,exist_ok=True);d.write_text(text);s=run(dict(case='op2_gradient',formulation=r['formulation'],deck=str(d),sha256=hashlib.sha256(d.read_bytes()).hexdigest()),P/'bin/mystran.exe',60,True);assert s['status']=='completed',s
 f=parse_f06(d.with_suffix('.F06'),'MYSTRAN');fib=np.array([v[1:5] for v in f['stress'][1] if abs(v[1])>0]);force=np.array([v[1:9] for v in f['forces'][1]])
 with (d.parent/'op2_reader.log').open('w') as log,contextlib.redirect_stdout(log):op=read_op2(str(d.with_suffix('.OP2')),debug=False,log=None)
 family='ctria3' if r['formulation'].startswith('CTRI') else 'cquad4';st=getattr(op,family+'_stress')[1].data[0,:,:4];ef=getattr(op.op2_results.force,family+'_force')[1].data[0,:,:8]
 assert st.shape==fib.shape and ef.shape==force.shape,(r['formulation'],st.shape,fib.shape)
 se=float(np.max(abs(st-fib))/max(1,np.max(abs(fib))));fe=float(np.max(abs(ef-force))/max(1,np.max(abs(force))));out.append(dict(formulation=r['formulation'],stress_relative_error=se,force_relative_error=fe,passed=se<2e-6 and fe<2e-6));assert out[-1]['passed'],out[-1]
(P/'op2_checks.json').write_text(json.dumps(out,indent=2));print(len(out),'F06/OP2 thermal output checks passed')

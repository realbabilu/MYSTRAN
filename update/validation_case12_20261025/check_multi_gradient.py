from pathlib import Path
import json,sys,hashlib
import numpy as np
P=Path(__file__).resolve().parent;sys.path.insert(0,'H:/Projects/MYSTRAN/shell_linear_benchmark_20261023');from run_suite import run;from results import parse_f06
out=[]
for r in json.loads((P/'first_run_results.json').read_text()):
 if r['tag']!='restrained_zero' or r['version']!='after':continue
 form=r['formulation'];tri=form.startswith('CTRI');n=3 if tri else 4;card='CTRIAR' if form.startswith('CTRIAR') else 'CTRIA3' if tri else 'CQUADR' if form.startswith('CQUADR') else 'CQUAD4'
 lines=Path(r['deck']).read_text().splitlines();param=next(l for l in lines if l.startswith('PARAM,'));coords=[[0,0,0],[2,0,0],[0,2,0]] if tri else [[0,0,0],[2,0,0],[2,2,0],[0,2,0]]
 text='SOL 101\nCEND\nECHO=NONE\nSPC=1\nSUBCASE 1\n TEMP(LOAD)=2\n STRESS(CENTER)=ALL\n FORCE(CENTER)=ALL\n DISPLACEMENT=ALL\nSUBCASE 2\n TEMP(LOAD)=3\n STRESS(CENTER)=ALL\n FORCE(CENTER)=ALL\n DISPLACEMENT=ALL\nBEGIN BULK\n'+param+'\nMAT1,1,29000.,,0.3,,6.5-6,50.\nPSHELL,1,1,1.,1,,1\n'
 for e in range(2):
  for j,c in enumerate(coords):text+=f'GRID,{e*n+j+1},,'+','.join(str(float(v+(10*e if k==0 else 0))) for k,v in enumerate(c))+'\n'
  ids=list(range(e*n+1,(e+1)*n+1));text+=f'{card},{e+1},1,'+','.join(map(str,ids))+'\nSPC1,1,12345,'+str(ids[0])+'\nSPC1,1,123456,'+','.join(map(str,ids[1:]))+'\n'
 text+='TEMPP1,2,1,150.,-100.\nTEMPP1,2,2,250.,60.\nTEMPP1,3,1,50.,100.\nTEMPP1,3,2,50.,-60.\nENDDATA\n'
 d=P/'multi_gradient'/form/'model.dat';d.parent.mkdir(parents=True,exist_ok=True);d.write_text(text);s=run(dict(case='mixed_element_temperatures_TREF',formulation=form,deck=str(d),sha256=hashlib.sha256(d.read_bytes()).hexdigest()),P/'bin/mystran.exe',60,True);assert s['status']=='completed',s;f=parse_f06(d.with_suffix('.F06'),'MYSTRAN');err=0
 for sc,vs in f['stress'].items():
  for eid,z,sx,sy,tau in vs:
   if abs(z)==0:continue
   mean=([100,200] if sc==1 else [0,0])[eid-1];grad=([-100,60] if sc==1 else [100,-60])[eid-1];expected=29000*6.5e-6/(1-.3)*(-mean+z*grad);err=max(err,abs(sx-expected),abs(sy-expected),abs(tau))
 out.append(dict(formulation=form,status=s['status'],max_stress_error=err,passed=err<1e-4));assert out[-1]['passed'],out[-1]
(P/'multi_gradient_checks.json').write_text(json.dumps(out,indent=2));print(len(out),'two-element / two-subcase / nonzero TREF checks passed')

from pathlib import Path
import sys,json,csv
import numpy as np
P=Path(__file__).resolve().parent;sys.path.insert(0,'H:/Projects/MYSTRAN/shell_linear_benchmark_20261023');from results import kgg
rows=json.loads((P/'results.json').read_text());assert len(rows)==396 and all(r['status']=='completed' for r in rows)
b={(r['formulation'],r['tag'],r['version']):r for r in rows};forms=sorted({r['formulation'] for r in rows});checks=[];ann=[]
def arrays(r,field):
 return np.array([v for sc,vals in sorted(r[field].items()) for k,v in sorted(vals.items(),key=lambda t:int(t[0]))]) if field=='displacements' else np.array([v for sc,vals in sorted(r[field].items()) for v in vals])
def test(name,form,err,tol):checks.append(dict(check=name,formulation=form,max_error=float(err),tolerance=float(tol),passed=bool(err<=tol)))
for form in forms:
 r=b[form,'annulus','after'];u=arrays(r,'displacements');amp=max(1e-12,np.max(abs(u)))
 for tag,factor in [('reverse',-1),('double',2),('zero',0)]:
  test(tag,form,np.max(abs(arrays(b[form,tag,'after'],'displacements')-factor*u)),3e-6*amp+1e-10)
 m=arrays(b[form,'mean_only','after'],'displacements');c=arrays(b[form,'combined','after'],'displacements');test('combined_superposition',form,np.max(abs(c-u-m)),3e-6*max(amp,np.max(abs(c)))+1e-10)
 for obs in r['observations']:
  val=obs['value'];ref=obs['reference'];ann.append(dict(formulation=form,quantity=obs['quantity'],value=val,reference=ref,error_pct=100*abs(val/ref-1)))
  test('annulus_'+obs['quantity'],form,abs(val/ref-1),.03)
 before=b[form,'restrained_zero','before'];after=b[form,'restrained_zero','after'];K0=kgg(Path(before['deck']).with_suffix('.OU1'));K1=kgg(Path(after['deck']).with_suffix('.OU1'));test('KGG_unchanged',form,np.linalg.norm(K1-K0)/np.linalg.norm(K0),1e-14)
 for tag,T,g in [('mean',100,0),('gradient',0,-100),('combined',100,-100),('reverse_gradient',0,100),('zero',0,0)]:
  q=b[form,'restrained_'+tag,'after'];C=29000*6.5e-6/(1-.3);err=0
  for z,sx,sy,tau in arrays(q,'stress')[:,1:]:
   expected=C*(-T+z*g);err=max(err,abs(sx-expected),abs(sy-expected),abs(tau))
  test('restrained_stress_'+tag,form,err,1e-4)
  st=arrays(q,'stress');lo=st[st[:,1]<0][0];hi=st[st[:,1]>0][0];moment=(hi[2:5]-lo[2:5])/12;forces=arrays(q,'forces')[0];test('fiber_moment_'+tag,form,np.max(abs(moment-forces[4:7])),2e-5)
 for tag in sorted(t for f,t,v in b if f==form and t.startswith('regression_') and v=='after'):
  old=b[form,tag,'before'];new=b[form,tag,'after'];test(tag+'_displacements',form,np.max(abs(arrays(old,'displacements')-arrays(new,'displacements'))),1e-10)
summary=dict(runs=len(rows),completed=sum(r['status']=='completed' for r in rows),checks=len(checks),passed=sum(r['passed'] for r in checks),all_passed=all(r['passed'] for r in checks),failures=[r for r in checks if not r['passed']])
(P/'validation.json').write_text(json.dumps(dict(summary=summary,checks=checks),indent=2));print(json.dumps(summary,indent=2))
with (P/'annulus_comparison.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(ann[0]));w.writeheader();w.writerows(ann)

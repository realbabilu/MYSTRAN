from pathlib import Path
import sys,re,json,csv,contextlib,io
import numpy as np
import run_regression as r
from scipy.optimize import linear_sum_assignment
P=r.P;OUT=r.OUT

def read(path):
 m=r.OP2(debug=False,log=r.SimpleLogger(level='error'),mode='msc');m.set_additional_result_tables_to_read({b'PVT0':False})
 m.set_results(['displacements','stress','eigenvectors','grid_point_stresses_volume_direct'])
 with contextlib.redirect_stdout(io.StringIO()):m.read_op2(str(path),build_dataframe=False)
 return m

def vectors(m):
 return [(k,{int(j):np.array(u[:3],float) for j,u in zip(a.node_gridtype[:,0],a.data[0])}) for k,a in m.displacements.items()]
def stress(m):
 return [a for a in m.op2_results.stress.chexa_stress.values()]
def eigen(path,m):
 txt=path.with_suffix('.F06').read_text(errors='replace').split('R E A L   E I G E N V A L U E S')[-1]
 vals=np.array([float(v) for v in re.findall(r'^\s+\d+\s+\d+\s+([+-]?\d+\.\d+E[+-]\d+)',txt,re.M)])
 a=next(iter(m.eigenvectors.values()))
 if getattr(a,'eigrs',None) is not None and (len(vals)!=len(a.data)):vals=np.asarray(a.eigrs,float)
 return vals,np.asarray(a.data[:,:,:3],float).reshape(len(a.data),-1)

def main():
 manifest=json.loads((OUT/'manifest.json').read_text());runs=list(csv.DictReader((OUT/'runs.csv').open()));status={(t['name'],t['solver']):t['status'] for t in runs};rows=[];parity=[];coverage=[]
 for job in manifest:
  name=job['name'];models={}
  nx=read(OUT/name/'nastran/model.op2') if status.get((name,'nastran'))=='COMPLETED' else None
  for solver in ['baseline','EAS9','TIED','TIED_AUTO']:
   if status.get((name,solver))!='COMPLETED':continue
   try:models[solver]=read(OUT/name/solver/'model.OP2')
   except Exception as exc:coverage.append(dict(name=name,solver=solver,status='OP2_READ_FAILED',detail=str(exc)));continue
   coverage.append(dict(name=name,solver=solver,status='OP2_READABLE',chexa_stress_rows=sum(len(a.element_node) for a in stress(models[solver]))))
  if nx is not None:models['nastran']=nx
  if 'baseline' in models and 'EAS9' in models:
   a,b=models['baseline'],models['EAS9'];err=0
   for (_,aa),(_,bb) in zip(vectors(a),vectors(b)):
    assert aa.keys()==bb.keys();err=max(err,max(float(np.max(abs(aa[j]-bb[j]))) for j in aa))
   serr=max((float(np.max(abs(aa.data-bb.data))) for aa,bb in zip(stress(a),stress(b))),default=0)
   if job.get('buckling'):
    av,am=eigen(OUT/name/'baseline/model.OP2',a);bv,bm=eigen(OUT/name/'EAS9/model.OP2',b);err=max(err,float(np.max(abs(av-bv))))
   parity.append(dict(name=name,displacement_or_eigenvalue_max_absolute_delta=err,stress_max_absolute_delta=serr,passed=err==0 and serr==0))
  for solver,m in models.items():
   if solver=='baseline':continue
   base=dict(case=job['case'],name=name,solver=solver,nastran_status='NEW_SOLVE_COMPLETED' if nx is not None else 'LICENSE_BLOCKED_NO_NEW_REFERENCE')
   vv=vectors(m)
   if nx is not None:
    nv=vectors(nx);deltas=[];norms=[]
    for (_,aa),(_,bb) in zip(vv,nv):
     assert aa.keys()==bb.keys()
     deltas.extend(float(np.linalg.norm(aa[j]-bb[j])) for j in aa);norms.extend(float(np.linalg.norm(bb[j])) for j in bb)
    base['nastran_displacement_relative_difference']=max(deltas,default=0)/max(max(norms,default=0),1e-30)
    coverage.append(dict(name=name,solver=solver,status='NEW_NASTRAN_COMPARISON',detail=str(base['nastran_displacement_relative_difference'])))
   if job['case']=='5-001':
    G=np.array(job['grads'][0]);xyz={int(j):np.array(x) for j,x in job['xyz'].items()};u=vv[0][1]
    uerr=max(float(np.max(abs(u[j]-G@xyz[j]))) for j in xyz)
    e=np.array([G[0,0],G[1,1],G[2,2],G[0,1]+G[1,0],G[1,2]+G[2,1],G[0,2]+G[2,0]]);D=np.zeros((6,6));D[:3,:3]=400000;np.fill_diagonal(D[:3,:3],1200000);D[3:,3:]=np.eye(3)*400000;s=D@e
    # OP2 solid stress ordering xx yy zz xy yz xz.
    a=stress(m)[0];serr=float(np.max(abs(a.data[0,:,:6]-s)))
    txt=(OUT/name/solver/'model.F06').read_text(errors='replace');tail=txt.split(' G R I D   P O I N T   S T R E S S E S')[-1]
    gp=[]
    for l in tail.splitlines():
     v=l.split()
     if len(v)==7 and v[0].isdigit():
      try:gp.append(list(map(float,v[1:])))
      except ValueError:pass
    direct=list(m.grid_point_stresses_volume_direct.values())
    if direct:gp=np.concatenate([a.data[0,:,:6] for a in direct]).tolist()
    gperr=float(np.max(abs(np.array(gp)-s))) if gp else None
    rows.append(dict(**base,response='affine_patch',displacement_max_absolute_error=uerr,corner_center_stress_max_absolute_error=serr,gpstress_max_absolute_error=gperr,gpstress_rows=len(gp),status='PASS_ANALYTICAL' if uerr<1e-9 and serr<1e-3 and len(gp)==16 and gperr<1e-3 else 'FAIL'))
   elif job['case'] in ['5-004','5-007']:
    for t in job['probes']:
     u=vv[t['sid']-1][1];v=float(np.mean([u[j][t['axis']] for j in t['joints']]))
     ref=([.005429,.001749][0 if t['pattern']=='IN' else 1] if job['case']=='5-004' else (.094 if t['axis']==0 else -.094))
     rows.append(dict(**base,response=t['pattern']+'_'+str(t['axis']),mystran=v,published_reference=ref,error_pct=100*abs(v-ref)/abs(ref),status='MEASURED_NO_NEW_NASTRAN'))
   elif job['case']=='5-008':
    u=vv[0][1];inner=job['inner'];xyz={int(j):np.array(x) for j,x in job['xyz'].items()};v=float(np.mean([u[j][:2]@(xyz[j][:2]/np.linalg.norm(xyz[j][:2])) for j in inner]));nu=job['nu'];ref=3*(1+nu)/job['E']*(81+(1-2*nu)*9)/72
    rows.append(dict(**base,response='inner_radial_mean',mystran=v,published_reference=ref,error_pct=100*abs(v-ref)/ref,status='MEASURED_ANALYTICAL_PLANE_STRAIN_NO_NEW_NASTRAN'))
    # Inner-face corner stresses transformed from BASIC into radial/hoop axes.
    samples=[]
    for a in stress(m):
     for (eid,gid),v in zip(a.element_node,a.data[0]):
      if int(gid) not in inner:continue
      n=xyz[int(gid)][:2]/np.linalg.norm(xyz[int(gid)][:2]);t=np.array([-n[1],n[0]])
      S=np.array([[v[0],v[3]],[v[3],v[1]]]);samples.append([n@S@n,t@S@t,v[2]])
    if samples:
     for label,value,ref in zip(['inner_sigma_radial','inner_sigma_hoop','inner_sigma_axial'],np.mean(samples,axis=0),[-1,1.25,.25*nu]):
      rows.append(dict(**base,response=label,mystran=float(value),published_reference=ref,error_pct=100*abs(value-ref)/abs(ref),status='MEASURED_ANALYTICAL_PLANE_STRAIN_NO_NEW_NASTRAN'))
   elif job['case']=='5-011':
    a=stress(m)[0];v=a.data[0,np.where(a.element_node[:,1]==0)[0],:6].mean(axis=0)
    for d in [0,1]:
     ref=-job['E']*job['alpha']*100/(1-job['nu']) if job['scope']=='restrained_uniform' else 0 if job['scope']=='free_uniform' else None
     rows.append(dict(**base,response='S'+'XY'[d]+'XY'[d]+'_center_mean',mystran=float(v[d]),published_reference=ref,status='PASS_ANALYTICAL' if ref is not None and abs(v[d]-ref)<.001 else 'MEASURED_GRADIENT' if ref is None else 'FAIL'))
   elif job['case']=='5-010':
    vals,vec=eigen(OUT/name/solver/'model.OP2',m);seen=set()
    # Modes separated by predominant transverse translation.
    for lam,v in zip(vals,vec):
     u=v.reshape(-1,3);key='weak_about_X' if np.linalg.norm(u[:,1])>np.linalg.norm(u[:,0]) else 'strong_about_Y'
     if key in seen:continue
     seen.add(key);ref=.888 if key.startswith('weak') else 127.91
     row=dict(**base,response=key,mystran=float(lam),published_reference=ref,error_pct=100*abs(lam-ref)/ref,status='MEASURED_PRELOAD_NO_NEW_NASTRAN')
     old=Path(job['source']).with_suffix('.op2')
     if old.exists():
      nm=read(old);na=next(iter(nm.eigenvectors.values()));nv=np.array(na.data[:,:,:3],float).reshape(len(na.data),-1);mac=(nv@v)**2/((nv*nv).sum(axis=1)*(v@v));ix=int(np.argmax(mac));row['archived_nastran_MAC']=float(mac[ix]);row['archived_nastran_eigenvalue']=float(na.eigrs[ix]);row['archived_reference_scope']='existing output only; current solve license-blocked'
     if nx is not None and solver!='nastran':
      na=next(iter(nx.eigenvectors.values()));nv=np.array(na.data[:,:,:3],float).reshape(len(na.data),-1)
      mac=(nv@v)**2/((nv*nv).sum(axis=1)*(v@v));ix=int(np.argmax(mac))
      row['nastran_MAC']=float(mac[ix]);row['nastran_matched_eigenvalue']=float(na.eigrs[ix])
     rows.append(row)
 matched={(x['name'],x['response']):x for x in rows if x['solver']=='nastran'}
 for row in rows:
  nr=matched.get((row['name'],row['response']))
  if nr and row.get('mystran') is not None:
   row['nastran']=nr['mystran'];row['nastran_error_pct']=nr.get('error_pct','')
   row['mystran_minus_nastran']=float(row['mystran'])-float(nr['mystran'])
   row['status']=row['status'].replace('NO_NEW_NASTRAN','MATCHED_NASTRAN')
 r.writecsv(OUT/'benchmark_matrix.csv',rows);r.writecsv(OUT/'eas9_baseline_parity.csv',parity);r.writecsv(OUT/'output_coverage.csv',coverage)
 print('results',len(rows),'baseline cases',len(parity),'baseline failures',[p for p in parity if not p['passed']]);print(json.dumps([t for t in rows if t.get('solver')=='TIED_AUTO' and t['case'] in ['5-004','5-007','5-010','5-008']],indent=2))
if __name__=='__main__':main()

from pathlib import Path
import json,sys,re,hashlib,shutil,concurrent.futures,csv
P=Path(__file__).resolve().parent;W=Path('H:/Projects/MYSTRAN/shell_linear_benchmark_20261023');sys.path[:0]=[str(W),str(W.parent/'codex_mod/mitc3_lee_nguyen_20261024')]
from run_suite import run
from results import parse_f06,measures
from build_suite import fields,cards
from report_tools import loaded_displacements_basic
manifest=json.loads((W/'manifest.json').read_text());jobs=[]
def add(caseid,r,label='',kmit=False,nd=None,sap=False):
 text=(W/r['deck']).read_text();out=[];ids=[int(c[1]) for c in cards(text) if c[0]=='GRID']
 for line in text.splitlines():
  f=fields(line)
  if f and f[0]=='PARAM' and len(f)>1 and f[1]=='DARPACK':continue
  if sap and f and f[0]=='SPC1' and f[2]=='246':continue
  if kmit and f and f[0]=='PARAM' and f[1]=='TRIARTYP':line='PARAM,TRIARTYP,KMITC3+'
  if nd and f and f[0]=='EIGRL':
   while len(f)<5:f.append('')
   f[4]=str(nd);line=','.join(f)
  f=fields(line)
  if f and f[0]=='SPC1' and 'THRU' not in f and len(f)>9:out.extend(','.join(f[:3]+f[i:i+6]) for i in range(3,len(f),6))
  else:out.append(line)
 if sap:
  j=out.index('ENDDATA');out[j:j]=['SPC1,1,246,'+','.join(map(str,ids[i:i+6])) for i in range(0,len(ids),6)]
 d=P/'runs'/f'case_{caseid:02d}'/(label or r['case'])/r['formulation']/'model.dat';d.parent.mkdir(parents=True,exist_ok=True);d.write_text('\n'.join(out)+'\n');jobs.append((caseid,label,dict(r,deck=str(d),sha256=hashlib.sha256(d.read_bytes()).hexdigest(),formulation='CTRIAR_KMITC3plus' if kmit else r['formulation'])))
for r in manifest:
 c=r['case'];form=r['formulation']
 if form=='CQUADR_DKMQ24' and c.startswith('prob_2_002_nx'):add(1,r)
 if form=='CTRIAR_DKMT18' and c.startswith('prob_2_002_nx'):add(2,r)
 if form=='CTRIAR_DKMT18' and c.startswith('prob_2_005_clamped'):add(3,r)
 if form=='CTRIA3_MITC3plus' and c.startswith('prob_2_005_simple'):add(4,r)
 if form=='CTRIAR_DKMT18' and c.startswith('prob_2_005_simple'):add(5,r)
 if c=='prob_2_007_thick_cquadr_simo' and form in ['CQUADR_DKMQ24','CQUAD4_DSQK','CTRIA3_DSG3','CTRIA3_MITC3plus']:
  did={'CQUADR_DKMQ24':6,'CQUAD4_DSQK':7,'CTRIA3_DSG3':8,'CTRIA3_MITC3plus':9}[form];add(did,r,'n8')
  for n in [16,32]:add(did,dict(r,deck=str(W/f'investigations/convergence_2_007/n{n}/{form}/model.dat')),f'n{n}')
 if form=='CTRIAR_DKMT18' and c.startswith('prob_2_010_'):add(10,r)
 if form=='CQUADR_DKMQ24' and c=='prob_2_013a_temp_cquadr_simo':add(11,r)
 if c=='prob_2_014b_tempgrad_cquadr_simo':add(12,r,kmit=form=='CTRIAR_MITC3plusHB')
 if c=='prob_2_016c_thick_beamseg_bend' and form in ['CTRIA3_DSG3','CTRIAR_MITC3plusHB']:
  did=13 if form=='CTRIA3_DSG3' else 14
  for nd in [1,6]:add(did,r,f'ND{nd}_default',kmit=did==14,nd=nd)
 if c.startswith('prob_2_017_nx') and form in ['CQUADR_DKMQ24','CTRIAR_MITC3plusHB']:
  add(15 if form=='CQUADR_DKMQ24' else 16,r,kmit=form=='CTRIAR_MITC3plusHB',nd=1,sap=True)
def task(j):
 did,label,r=j;s=run(r,P/'bin/mystran.exe',300,True);data=parse_f06(Path(r['deck']).with_suffix('.F06'),'MYSTRAN');data['stress']={sc:[v for v in vals if abs(v[1])>0] for sc,vals in data['stress'].items()};obs=measures(r,data)
 if r['group']=='2007' and s['status']=='completed':
  b=loaded_displacements_basic(Path(r['deck']),data)
  for m in obs:
   for f in cards(Path(r['deck']).read_text()):
    if f[0]=='FORCE':
     dof=next(i for i in range(3) if float(f[5+i])!=0)
     if m['quantity']=='loaded_U'+('X','Y','Z')[dof]:m['value']=float(b[int(f[2])][dof])
 return dict(discussion=did,label=label,case=r['case'],formulation=r['formulation'],status=s['status'],errors=s['errors'],deck=r['deck'],sha256=r['sha256'],solver_sha256=s['solver_sha256'],observations=obs,modes=data['modes'])
if __name__=='__main__':
 shutil.copy2('H:/Projects/MYSTRAN/source/Binaries/libopenblas.dll',P/'bin/libopenblas.dll');results=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  for res in pool.map(task,jobs):
   results.append(res);(P/'results.json').write_text(json.dumps(results,indent=2))
   if len(results)%10==0:print(len(results),'/',len(jobs),'completed',sum(r['status']=='completed' for r in results),flush=True)
 with (P/'comparison.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['discussion','label','case','formulation','status','loadcase','quantity','computed','reference','ratio','deck']);w.writeheader()
  for r in results:
   for m in r['observations']:w.writerow(dict(discussion=r['discussion'],label=r['label'],case=r['case'],formulation=r['formulation'],status=r['status'],loadcase=m['loadcase'],quantity=m['quantity'],computed=m['value'],reference=m['reference'],ratio=m['value']/m['reference'] if m['reference'] and m['value'] is not None else '',deck=r['deck']))
 print('FINAL',len(results),{s:sum(r['status']==s for r in results) for s in set(r['status'] for r in results)},flush=True)

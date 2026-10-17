import numpy as np
import check_features as f
import check_operators as o
import collect_regression as c
import run_regression as r
src=(r.OUT/'5007_8x8/TIED_AUTO/model.dat').read_text();rows=[]
base=c.vectors(c.read(r.OUT/'5007_8x8/TIED_AUTO/model.OP2'))[0][1]
for ax in [0,1]:
 lines=src.splitlines();out=[];i=0;perm=o._perm_for(ax)
 while i<len(lines):
  l=lines[i];v=l.split(',')
  if v[0]=='CHEXA':
   fields=v[1:];i+=1
   while i<len(lines) and lines[i].startswith('+,'):fields+=lines[i].split(',')[1:];i+=1
   cc=fields[2:];assert len(cc)==8
   out.append(r.card('CHEXA',fields[:2]+[cc[k] for k in perm]));continue
  out.append(l);i+=1
 ff,d=f.run('dome8_cyclic_'+str(ax),'\n'.join(out)+'\n');my=c.vectors(c.read(d/'model.OP2'))[0][1]
 delta=max(float(np.max(abs(my[j]-base[j]))) for j in base);assert delta<1e-8,delta
 rows.append(dict(test='dome8_AUTO_cyclic_'+str(ax),max_absolute_displacement_delta=delta,status='PASS'))
r.writecsv(f.OUT/'orientation_matrix.csv',rows);print(rows)

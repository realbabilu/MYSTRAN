import numpy as np
from collections import Counter
from pathlib import Path
from itertools import product
P=Path(__file__).resolve().parent
XYZ={1:(.249,.342,.192),2:(.626,.288,.288),3:(.850,.649,.263),4:(.273,.750,.230),5:(.320,.186,.643),6:(.677,.305,.683),7:(.788,.693,.644),8:(.165,.745,.702),9:(0,0,0),10:(1,0,0),11:(1,1,0),12:(0,1,0),13:(0,0,1),14:(1,0,1),15:(1,1,1),16:(0,1,1)}
HEX={1:[9,1,12,4,13,5,16,8],2:[9,10,1,2,13,14,5,6],3:[1,2,4,3,5,6,8,7],4:[4,3,12,11,8,7,16,15],5:[2,10,3,11,6,14,7,15],6:[5,6,8,7,13,14,16,15],7:[9,10,12,11,1,2,4,3]}
BIN2NAS=[0,1,3,2,4,5,7,6]
HED=[(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
TED=[(0,1),(1,2),(0,2),(0,3),(1,3),(2,3)]
PED=[(0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(0,3),(1,4),(2,5)]
YED=[(0,1),(1,2),(2,3),(3,0),(0,4),(1,4),(2,4),(3,4)]
TRIS={0:[(0,1,2),(0,2,3)],1:[(1,2,3),(1,3,0)]}
DIAG={0:(0,2),1:(1,3)}
def bi(c):return c[0]+2*c[1]+4*c[2]
def face(a,v):
 p,q=[k for k in (0,1,2) if k!=a]; out=[]
 for pv,qv in [(0,0),(1,0),(1,1),(0,1)]:
  c=[0,0,0]; c[a],c[p],c[q]=v,pv,qv;out.append(bi(c))
 return out
class Mesh:
 def __init__(self):self.xyz={k:np.array(v,float) for k,v in XYZ.items()};self.n=100;self.keys={};self.els=[]
 def avg(self,ids,tag):
  key=(tag,tuple(sorted(ids)))
  if key not in self.keys:
   self.xyz[self.n]=np.mean([self.xyz[k] for k in ids],axis=0);self.keys[key]=self.n;self.n+=1
  return self.keys[key]
 def mid(self,a,b):return self.avg([a,b],'E')
 def add(self,n,ed=None):
  if ed:n=n+[self.mid(n[i],n[j]) for i,j in ed]
  self.els.append(n)
def split_wedges():
 def opts(h):
  b=HEX[h];out={}
  for a,bit in product(range(3),range(2)):
   out[a,bit]={frozenset(b[i] for i in face(a,v)):frozenset(b[face(a,v)[i]] for i in DIAG[bit]) for v in (0,1)}
  return out
 os={h:opts(h) for h in HEX}; faces={h:{frozenset(HEX[h][i] for i in face(a,v)) for a,v in product(range(3),range(2))} for h in HEX};hs=list(HEX)
 def rec(sel):
  if len(sel)==len(hs):return sel.copy()
  h=hs[len(sel)]
  for o,diags in os[h].items():
   if all(all(diags.get(k)==os[hh][oo].get(k) for k in faces[h]&faces[hh]) for hh,oo in sel.items()):
    sel[h]=o;r=rec(sel)
    if r:return r
    del sel[h]
  return None
 ans=rec({});assert ans;return ans
SPLIT=split_wedges()
def orient_tet(n,xyz):
 v=np.array([xyz[i] for i in n]);d=np.linalg.det(np.column_stack((v[1]-v[0],v[2]-v[0],v[3]-v[0])))
 return n if d>0 else [n[0],n[1],n[3],n[2]]
def orient_penta(n,xyz):
 p=np.array([xyz[i] for i in n]);d=np.linalg.det(np.column_stack((p[1]-p[0],p[2]-p[0],p[3]-p[0])))
 return n if d>0 else [n[0],n[2],n[1],n[3],n[5],n[4]]
def orient_pyramid(n,xyz):
 p=np.array([xyz[i] for i in n]);d=np.linalg.det(np.column_stack((p[1]-p[0],p[3]-p[0],p[4]-p[0])))
 return n if d>0 else [n[0],n[3],n[2],n[1],n[4]]
def build(typ):
 m=Mesh();quad=typ.endswith(('10','15','20','13'))
 for h,b in HEX.items():
  if typ.startswith('CHEXA'):
   m.add([b[i] for i in BIN2NAS],HED if quad else None)
  elif typ.startswith('CPENTA'):
   a,bit=SPLIT[h];p,q=[i for i in (0,1,2) if i!=a];corn=[]
   for pv,qv in [(0,0),(1,0),(1,1),(0,1)]:
    row=[]
    for av in (0,1):
     c=[0,0,0];c[a],c[p],c[q]=av,pv,qv;row.append(b[bi(c)])
    corn.append(row)
   for tr in TRIS[bit]:
    n=orient_penta([corn[t][0] for t in tr]+[corn[t][1] for t in tr],m.xyz)
    m.add(n,PED if quad else None)
  elif typ.startswith('CTETRA'):
   c=m.avg(b,'H')
   for a,v in product(range(3),range(2)):
    f=[b[i] for i in face(a,v)];fc=m.avg(f,'F')
    for i in range(4):m.add(orient_tet([c,fc,f[i],f[(i+1)%4]],m.xyz),TED if quad else None)
  else:
   c=m.avg(b,'H')
   for a,v in product(range(3),range(2)):
    f=[b[i] for i in face(a,v)];n=orient_pyramid(f+[c],m.xyz)
    m.add(n,YED if quad else None)
 return m
TYPES=['CTETRA4','CTETRA10','CPENTA6','CPENTA15','CHEXA8','CHEXA20','CPYRAM5','CPYRAM13']
G=np.array([[2,1,1],[1,2,1],[1,1,2]],float)/2000
# Nastran large-field, 16-wide float fields, continuation by '*' free large-field continuation.
def card(name,fields):
 # Prefer free-field Nastran BDF to avoid truncating decimal SPCD values.
 vals=[str(x) for x in fields]
 parts=[vals[i:i+7] for i in range(0,len(vals),7)]
 return '\n'.join([name+','+','.join(parts[0])]+['+,'+','.join(x) for x in parts[1:]])+'\n'

def main():
 lines=['$ SAP2000 Problem 5-001 / MacNeal-Harder - Eight independent solid patches\n',
 '$ DO NOT CONNECT THE PATCHES. One homogeneous PSOLID / MAT1 per group.\n',
 '$ Enforcement: SPC1 on surface + SPCD in static LOAD case.\n',
 '$ Note: validate SPCD nonzero support on the selected MYSTRAN build.\n',
 'SOL 101\n','CEND\n','TITLE = 5-001 EIGHT INDEPENDENT SOLID PATCHES\n',
 'SUBCASE 1\n','  LABEL = AFFINE_DISPLACEMENT_PATCH\n','  SPC = 1\n','  LOAD = 1\n','  DISPLACEMENT(PLOT) = ALL\n','  STRESS(PLOT) = ALL\n','BEGIN BULK\n','PARAM,POST,-1\n','MAT1,1,1.0E6,,0.25,1.0\n']
 csv=['type,box_index,elements,nodes,boundary_nodes,monitor_node,x_offset,y_offset,monitor_ux,monitor_uy,monitor_uz']
 all_m=[];spcd=[]
 for i,typ in enumerate(TYPES):
  m=build(typ);dx=(i%4)*2.0;dy=(i//4)*2.0;gidbase=(i+1)*100000;eidbase=(i+1)*100000;pid=i+1
  used=sorted({nid for ns in m.els for nid in ns});bound=[]
  lines += [f'$ -------- BOX {i+1}: {typ}, {len(m.els)} elements --------\n',f'PSOLID,{pid},1\n']
  for nid in used:
   pt=m.xyz[nid]; gid=gidbase+nid; xy=pt+np.array([dx,dy,0.]);lines.append(card('GRID',[gid,'',*(f'{v:.12g}' for v in xy)]))
   if np.any(np.isclose(pt,0,atol=1e-11)) or np.any(np.isclose(pt,1,atol=1e-11)):
    bound.append(gid);u=G@pt
    for component,val in enumerate(u,1):
     if abs(val)>1e-14:spcd.append((gid,component,val))
  elemname=typ.rstrip('0123456789')
  if typ.startswith('CPYRAM'):elemname='CPYRAM'
  for j,ns in enumerate(m.els,1):lines.append(card(elemname,[eidbase+j,pid]+[gidbase+n for n in ns]))
  # Some GRID boundary dofs have prescribed exact zero; put entire face into SPC1 (all 3 comps)
  lines.append(f'$ Boundary nodes for {typ}: {len(bound)}\n')
  for k in range(0,len(bound),7):lines.append(card('SPC1',[1,123]+bound[k:k+7]))
  csv.append(f'{typ},{i+1},{len(m.els)},{len(used)},{len(bound)},{gidbase+8},{dx},{dy},0.00089,0.00118,0.00116')
  all_m.append((typ,m,dx,dy,gidbase))
 lines.append('$ SPCD applied via LOAD = 1; exact displacement u=G*(x-offset)\n')
 for gid,c,u in spcd:lines.append(card('SPCD',[1,gid,c,f'{u:.12g}']))
 lines.append('ENDDATA\n')
 # Stress and boundary expected; SPCD load set shares SID with LOAD case.
 (P/'5-001_8_SOLIDS.dat').write_text(''.join(lines))
 (P/'5-001_reference.csv').write_text('\n'.join(csv)+'\n')
 readme='''# 5-001 eight independent solid patches

Eight translates of the MacNeal-Harder 7-irregular-hexa unit-cube geometry. Original XYZ and connectivity copied from the project's `solid3d_patch_test_macneal.py`. Each hexa is subdivided conformingly: CTETRA 24/hex, CPENTA 2/hex, CPYRAM 6/hex. For CPYRAM13 add exactly 8 edge midpoint nodes; no face center (that would make CPYRAM14). Shared edges use a unique midpoint GRID.

## Deck

`5-001_8_SOLIDS.dat` SOL 101, 8 disconnected regions, offsets (0,0), (2,0), (4,0), (6,0), (0,2), (2,2), (4,2), (6,2). Local deformation u=G x_local, G=(1/2000)*[[2,1,1],[1,2,1],[1,1,2]]. E=1e6, nu=0.25. Expected per-region S=[2000,2000,2000,400,400,400], local joint 8 u=[0.00089,0.00118,0.00116].

All true boundary nodes (including midside nodes) are constrained by SPC1 123 and nonzero imposed components added using SPCD SID=1 with case control LOAD=1. **Important**: This tests MYSTRAN's enforced-displacement implementation too; validate whether your MYSTRAN SOL101 accepts nonzero SPCD specified as a static LOAD SID. If unsupported use an equivalent enforced-displacement method or a known-good solver; never silently replace nonzero SPCD by zero SPC. SPCD and SPC1 sharing SID=1 is intentional (different namespaces). Avoid false passes: validate free interior nodes, stress at integration points and element orientations, not merely prescribed boundary nodes.

Nastran free-field continuation is used (`+,`), and every free-field GRID/element/SPCD is parseable by pyNastran. CPYRAM13 nodal order adopted: 1-4 base, 5 apex, 6-9 base edge midsides, 10-13 apex edge midsides. Check the corresponding MYSTRAN CPYRAM13 connectivity convention.

## Acceptance

For each region compare analytical displacement on all **unconstrained** nodes, and all element integration-point stresses. Tolerance initially 1e-7 relative (floating-point) for pure linear consistency patch. Where formulation under evaluation has known deficiencies, report FAIL transparently. Test results are not asserted here: generated decks are unexecuted.

Run `python make_5001.py` to rebuild. Reference index `5-001_reference.csv`.
'''
 (P/'README.md').write_text(readme)
 print('\n'.join(csv));print('deck lines',len(lines),'spcd',len(spcd))
if __name__=='__main__':main()

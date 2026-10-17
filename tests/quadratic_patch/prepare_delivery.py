"""Prepare explicit MYSTRAN-order equivalents of the user's fixed DAT."""
from pathlib import Path
import shutil

ROOT=Path(__file__).resolve().parent
WORK=ROOT.parents[2]
DEST=WORK/'nastran_test/mystran/quadratic_patch_audit'
DEST.mkdir(exist_ok=True)
fixed=WORK/'nastran_test/mystran/5-001_8_SOLIDS_GPSTRESS_PENTA_REMESH_FIXED_QUADRATIC.dat'
lines=fixed.read_text().splitlines();out=[];i=0
while i<len(lines):
    line=lines[i]
    if not line.startswith(('CHEXA,','CPENTA,')):
        out.append(line);i+=1;continue
    typ,*fields=line.split(',');i+=1
    while i<len(lines) and lines[i].startswith('+'):
        fields+=lines[i].split(',')[1:];i+=1
    fields=[f for f in fields if f.strip()]
    nodes=fields[2:]
    if typ=='CHEXA' and len(nodes)==20:
        nodes=nodes[:12]+nodes[16:20]+nodes[12:16]
    elif typ=='CPENTA' and len(nodes)==15:
        nodes=nodes[:9]+nodes[12:15]+nodes[9:12]
    fields=fields[:2]+nodes
    out.append(typ+','+','.join(fields[:8]))
    for j in range(8,len(fields),8):out.append('+,'+','.join(fields[j:j+8]))
(DEST/'5-001_FIXED_QUADRATIC_MYSTRAN_ORDER.dat').write_text('\n'.join(out)+'\n')
for pid,name in [(4,'CPENTA15_PATCH'),(6,'CHEXA20_PATCH')]:
    shutil.copy2(ROOT/f'pid{pid}_manual.dat',DEST/f'{name}.dat')
for name in ['python_patch.csv','python_shape_checks.csv','fortran_matrices_before.csv','fortran_matrices_after.csv',
             'fortran_patch_before.csv','fortran_patch_after.csv','mystran_solver_patch.csv','mystran_solver_gauss.csv',
             'mystran_solver_nodal.csv','node_dof_ordering.csv','boundary_constraints.csv','quadrature_weights.csv']:
    if (ROOT/name).exists():shutil.copy2(ROOT/name,DEST/name)
print(DEST)

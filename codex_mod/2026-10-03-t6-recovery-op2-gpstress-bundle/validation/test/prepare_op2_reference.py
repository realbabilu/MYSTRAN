from pathlib import Path
import re
root=Path(__file__).resolve().parents[1]
folder=root/'test/t6_upgrade/op2_reference'
folder.mkdir(parents=True,exist_ok=True)
s=(root/'test/t6_upgrade/stress_recovery/MITC6.dat').read_text()
s=s.replace('DISPLACEMENT = ALL','DISPLACEMENT(PLOT,PRINT) = ALL')
s=s.replace('STRESS(CENTER,CORNER)','STRESS(PLOT,PRINT,CENTER,CORNER)')
s=s.replace('FORCES(CENTER,CORNER)','FORCES(PLOT,PRINT,CENTER,CORNER)')
s=s.replace('BEGIN BULK','BEGIN BULK\nPARAM,POST,-1')
(folder/'mystran_t6.dat').write_text(s)
n=re.sub(r'^PARAM,TRIA6TYP,.*\n','',s,flags=re.M)
n=n.replace('STRESS(PLOT,PRINT,CENTER,CORNER)','STRESS(PLOT,PRINT,CORNER)')
n=n.replace('FORCES(PLOT,PRINT,CENTER,CORNER)','FORCE(PLOT,PRINT,CORNER)')
n=n.replace('STRFIELD=ALL','$STRFIELD=ALL')
n=re.sub(r'^\$?DEBUG,.*\n','',n,flags=re.M)
(folder/'nastran_t6.dat').write_text(n)
print(folder)

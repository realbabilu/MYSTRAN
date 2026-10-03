"""Check standalone CENTER/CORNER and hidden GPSTRESS recovery requests."""
from pathlib import Path
import sys, os, subprocess, json
root = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'python'))
from compare_f06_vs_python import parse_mystran
folder = root/'test/t6_upgrade/output_selection'
folder.mkdir(parents=True,exist_ok=True)
env = os.environ.copy()
env['PATH'] = str(root/'test')+';C:\\gcc\\bin;'+env['PATH']
results = []
for family in ('SIMOT6','MITC6','MH6T','REZAIEE'):
    base = (root/f'test/t6_upgrade/stress_recovery/{family}.dat').read_text()
    for style,gp in (('CENTER',False),('CORNER',False),('CENTER',True),('GPONLY',True)):
        source = base.replace('STRESS(CENTER,CORNER)',f'STRESS({style})')
        source = source.replace('FORCES(CENTER,CORNER)', 'FORCES(CENTER)')
        if not gp:
            source = source.replace('GPSTRESS=ALL','$GPSTRESS=ALL')
        if style == 'GPONLY':
            source = source.replace('  STRESS(GPONLY)= ALL', '$ STRESS OMITTED')
            source = source.replace('  STRESS(GPONLY) = ALL', '$ STRESS OMITTED')
        deck = folder/f'{family}_{style}_{gp}.dat'
        deck.write_text(source)
        run = subprocess.run([str(root/'MYSTRAN/Binaries/mystran.exe'),deck.name],cwd=folder,env=env,capture_output=True,text=True,timeout=90)
        deck.with_suffix('.stdout.txt').write_text(run.stdout+run.stderr)
        D = parse_mystran(str(deck.with_suffix('.F06')))
        for sc in (1,2):
            counts = [len(D['elstr'].get(sc,{})),sum(len(v) for v in D['elgrd'].get(sc,{}).values()),len(D['gp'].get(sc,{}))]
            expected = [10 if style=='CENTER' else 0,60 if style=='CORNER' else 0,25 if gp else 0]
            record = dict(family=family,style=style,gp=gp,sc=sc,counts=counts,expected=expected,
                          pass_gate=run.returncode==0 and counts==expected)
            results.append(record)
            print(json.dumps(record),flush=True)
(folder/'results.json').write_text(json.dumps(results,indent=2))
print('PASS',sum(r['pass_gate'] for r in results),'/',len(results))

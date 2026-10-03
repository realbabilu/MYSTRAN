"""Check all four T6 OP2 files and the Siemens native table contract."""
from pathlib import Path
import subprocess,sys,json
root=Path(__file__).resolve().parents[1]
folder=root/'test/t6_upgrade'
for family in ('SIMOT6','MITC6','MH6T','REZAIEE'):
    path=folder/f'stress_recovery/{family}.OP2'
    run=subprocess.run([sys.executable,str(root/'test/read_t6_op2.py'),str(path)],capture_output=True,text=True)
    path.with_suffix('.pynastran.txt').write_text(run.stdout+run.stderr,encoding='utf-8')
    assert run.returncode==0,run.stdout+run.stderr
    summary=json.loads(path.with_suffix('.pynastran.json').read_text())
    for key in ('stress_SC1','stress_SC2'):
        assert summary[key]['element_type']==75 and summary[key]['num_wide']==70
        assert summary[key]['shape']==[1,80,8] and summary[key]['f06_samples']==80
        assert summary[key]['elements']==list(range(51,61))
    assert summary['force_SC2']['num_wide']==38 and summary['force_SC2']['shape']==[1,40,8]
    assert summary['force_SC2']['f06_samples']==40
    for sc in (1,2): assert summary[f'GPSTRESS_600_SC{sc}']['f06_samples']==50
    print(family,'native OES/OEF/OGS + F06 values PASS')
ref=folder/'op2_reference'
reference=json.loads((ref/'nastran_t6.pynastran.json').read_text())
actual=json.loads((ref/'mystran_t6.pynastran.json').read_text())
for key in ('stress_SC1','stress_SC2','force_SC2'):
    for field in ('element_type','num_wide','shape','elements'):
        assert reference[key][field]==actual[key][field],(key,field)
print('Siemens Nastran CTRIA6 table layout and element IDs MATCH')

"""Check bundle integrity; optionally compare snapshots to the original workspace."""
from pathlib import Path
import argparse,ast,hashlib,json,subprocess
bundle=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser()
ap.add_argument('--workspace',type=Path)
args=ap.parse_args()
manifest=json.loads((bundle/'manifest.json').read_text(encoding='utf-8'))
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
count=0
for line in (bundle/'CHECKSUMS.sha256').read_text().splitlines():
    expected,relative=line.split('  ',1)
    path=(bundle/relative).resolve()
    assert path.is_relative_to(bundle),relative
    assert sha(path)==expected,relative
    count+=1
for record in manifest['files']:
    snapshot=bundle/record['bundle_path']
    assert snapshot.stat().st_size==record['bytes'] and sha(snapshot)==record['sha256'],record['bundle_path']
    if args.workspace:
        original=args.workspace.resolve()/record['original']
        assert sha(original)==record['sha256'],str(original)
python_files=list((bundle/'validation').rglob('*.py'))+list((bundle/'tools').glob('*.py'))
for path in python_files: ast.parse(path.read_text(encoding='utf-8-sig'),filename=str(path))
if args.workspace:
    repo=args.workspace.resolve()/'MYSTRAN'
    names=subprocess.check_output(['git','diff','--name-only',manifest['baseline'],'--','Source'],cwd=repo).decode().splitlines()
    assert sorted(names)==sorted(manifest['tracked_modified']), 'Source delta changed after packaging'
    subprocess.run(['git','apply','--reverse','--check',str(bundle/'SOURCE_FROM_931dde4.patch')],cwd=repo,check=True)
print('PASS:',count,'file checksums;',len(manifest['files']),'snapshot records;',len(python_files),'Python syntax checks; baseline source patch verified' if args.workspace else 'PASS: bundle integrity and Python syntax')

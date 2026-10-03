"""Restore sibling validation assets without silently overwriting local work."""
from pathlib import Path
import argparse,shutil
bundle=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser()
ap.add_argument('--workspace',type=Path,required=True)
ap.add_argument('--overwrite',action='store_true')
ap.add_argument('--evidence',action='store_true')
args=ap.parse_args()
workspace=args.workspace.resolve()
assert (workspace/'MYSTRAN/Source').is_dir(),'Expected workspace/MYSTRAN/Source; clone the repo as MYSTRAN first.'
pairs=[]
for folder in ('python','test'):
    for src in (bundle/'validation'/folder).rglob('*'):
        if src.is_file(): pairs.append((src,workspace/folder/src.relative_to(bundle/'validation'/folder)))
if args.evidence:
    for src in (bundle/'evidence').rglob('*'):
        if src.is_file(): pairs.append((src,workspace/'test'/src.relative_to(bundle/'evidence')))
conflicts=[str(dst) for src,dst in pairs if dst.exists() and src.read_bytes()!=dst.read_bytes()]
if conflicts and not args.overwrite:
    raise SystemExit('Different existing files; no copies performed. Review before --overwrite:\n'+'\n'.join(conflicts))
for src,dst in pairs:
    dst.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(src,dst)
print('Restored',len(pairs),'validation/evidence files into',workspace)

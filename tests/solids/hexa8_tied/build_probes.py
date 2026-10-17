"""Build the production and frozen-baseline operator probes with the configured Windows toolchain."""
from pathlib import Path
import os,subprocess
P=Path(__file__).resolve().parent;ROOT=P.parents[2]
env=os.environ.copy();env['PATH']='C:/gcc/bin;C:/gcc/openblas32/bin;'+env['PATH']
common=['C:/gcc/bin/gfortran.exe','-static','-fcheck=all','-ffree-line-length-none','-I',str(ROOT/'Binaries/mod'),'-J',str(P)]
libs=[str(ROOT/'build/CMakeFiles/mystran.dir/objects.a'),'@'+str(ROOT/'build/CMakeFiles/mystran.dir/linkLibs.rsp')]
for name,sources in [('probe',['probe.f90']),('baseline_probe',['baseline_module.f90','baseline_probe.f90'])]:
 subprocess.run(common+[str(P/s) for s in sources]+libs+['-o',str(P/(name+'.exe'))],cwd=ROOT,env=env,check=True)

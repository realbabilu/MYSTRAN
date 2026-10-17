"""Build and reproduce the Python -> kernel -> full MYSTRAN acceptance audit."""
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT.parents[1]
env=os.environ.copy();env['PATH']='C:/gcc/bin;'+str(SOURCE/'Binaries')+os.pathsep+env['PATH']
subprocess.run(['cmake','--build',str(SOURCE/'build'),'-j','8'],check=True,env=env)
subprocess.run(['C:/gcc/bin/gfortran.exe','-static','-fcheck=all','-ffree-line-length-none',
                '-I',str(SOURCE/'Binaries/mod'),str(ROOT/'quadratic_probe.f90'),
                str(SOURCE/'build/CMakeFiles/mystran.dir/objects.a'),
                '@'+str(SOURCE/'build/CMakeFiles/mystran.dir/linkLibs.rsp'),'-o',str(ROOT/'quadratic_probe.exe')],
               check=True,cwd=SOURCE,env=env)
shutil.copy2(SOURCE/'Binaries/mystran.exe',ROOT/'mystran_audit.exe')
(ROOT/'MYSTRAN.INI').write_text(f'{"L5A":8s}{"":8s}{"KEEP":8s}\n')
for phase in ([],['after'],['geometry'],['solver']):
    subprocess.run([sys.executable,str(ROOT/'audit_quadratic.py'),*phase],check=True,cwd=ROOT,env=env)
subprocess.run([sys.executable,str(ROOT/'prepare_delivery.py')],check=True,cwd=ROOT,env=env)

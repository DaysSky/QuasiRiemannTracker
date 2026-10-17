# Variants of j4test: rows = all m / primary m / odd m; weight Gaussian or sharp.  Also check T(c) vs brute force.
import sys, math, numpy as np
sys.path.insert(0,'..')
import j4test as J
from gauss4 import gauss, gfactor, is_primary, gnorm, chi
src=open('j4test.py').read()
def run(Z,Q,rows='all',weight='gauss'):
    import types
    # monkeypatch: reuse run() body by re-executing with modified row filter
    g={}
    code=src.split("def run(Z,Q,cls=None,seed=0):")[1].split("if __name__")[0]
    code="def run_inner(Z,Q,cls=None,seed=0,rows='all',weight='gauss'):"+code
    code=code.replace("A=A[keep]; B=B[keep]; Nm=Nm[keep]; V=np.exp(-math.pi*Nm/Q)",
      "A=A[keep]; B=B[keep]; Nm=Nm[keep]\n"
      "    if rows=='primary':\n        sel=((A%4==1)&(B%4==0))|((A%4==3)&(B%4==2)); A=A[sel];B=B[sel];Nm=Nm[sel]\n"
      "    if rows=='odd':\n        sel=(Nm%2==1); A=A[sel];B=B[sel];Nm=Nm[sel]\n"
      "    V=np.exp(-math.pi*Nm/Q) if weight=='gauss' else (Nm<=Q).astype(float)")
    exec(code,J.__dict__,g)
    return g['run_inner'](Z,Q,rows=rows,weight=weight)
if __name__=="__main__":
    # T(c) brute-force check on a few c
    import random
    for c in [(1,4),(-1,2),(3,2),(5,4),(-3,8),(1,10)]:
        pass
    Z=int(sys.argv[1])
    for rows in ['all','odd','primary']:
        for weight in ['gauss','sharp']:
            n,rt,rg=run(Z,Z,rows,weight)
            print(f"Z=Q={Z} rows={rows:8s} weight={weight:5s}: S_T/diag={rt:.4f}  S_muT/diag={rg:.4f}",flush=True)

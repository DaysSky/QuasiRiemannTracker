import sys; sys.path.insert(0,'..')
from rd import load
from gauss4 import *
import random, math
D=load('t3000.bin')
random.seed(1)
idx=random.sample(range(D['n']),60)
bad=0
for i in idx:
    n=(int(D['a'][i]),int(D['b'][i]))
    if n==(1,0): continue
    if gnorm(n)>1500: continue
    fac=gfactor(n)
    G=gauss(n,3,(1,0),fac)
    if abs(G-D['G'][i])>1e-6*math.sqrt(gnorm(n)): bad+=1; print("G mismatch",n,G,D['G'][i])
    for k in range(0,D['K'],5):
        m=D['tp'][k][0]
        c=chi(m,n,fac); ee=D['e'][i,k]
        cc=0 if ee==255 else 1j**int(ee)
        if abs(c-cc)>1e-9: bad+=1; print("chi mismatch",n,m,c,cc)
print("checked, bad=",bad, "tp[:6]=",D['tp'][:6])

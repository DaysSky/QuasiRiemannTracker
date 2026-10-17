import sys, math, numpy as np
import j4test3 as J
src=open('j4test3.py').read()
code=src.split("def run(Z,Q,rows='primary',weight='sharp',inert=True):")[1].split("if __name__")[0]
code=code.replace("RT+=Tc*w*val; RG+=mu*Tc*w*val","al=complex(*c)/math.sqrt(N)\n        for key,beta in (('T',Tc),('muT',mu*Tc),('abarT',Tc/al)):\n            ACC[key+'_cbar']+=beta*w*val; ACC[key+'_c']+=beta*w*np.conj(val)")
code=code.replace("RT=np.zeros(len(A),dtype=complex); RG=np.zeros(len(A),dtype=complex); diag=0.0","ACC={k:np.zeros(len(A),dtype=complex) for k in ['T_cbar','T_c','muT_cbar','muT_c','abarT_cbar','abarT_c']}; diag=0.0")
code=code.replace("return len(cs),np.sum(V*abs(RT)**2)/diag,np.sum(V*abs(RG)**2)/diag","return {k:np.sum(V*abs(v)**2)/diag for k,v in ACC.items()}")
g={}
exec("def run(Z,Q,rows='primary',weight='sharp',inert=True):"+code,J.__dict__,g)
for Z in [int(a) for a in sys.argv[1:]]:
    r=g['run'](Z,Z)
    print(Z,{k:round(v,5) for k,v in r.items()},flush=True)

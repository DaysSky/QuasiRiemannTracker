# Build Dirichlet-series coefficient data for index r = m^j and evaluate smooth sums.
import sys, math, numpy as np
sys.path.insert(0,'..')
from rd import load
from gauss4 import gauss, gfactor, gnorm, is_primary

def W(x):
    x=np.asarray(x,float); out=np.zeros_like(x)
    m=(x>0.5)&(x<1.0); y=4*x[m]-3
    out[m]=np.exp(1-1/(1-y*y)); return out

class Data:
    def __init__(s,fn):
        D=load(fn); s.__dict__.update(D)
        s.Nf=s.N.astype(float); s.alpha=(s.a+1j*s.b)/np.sqrt(s.Nf)
        s.Gt=s.G/np.sqrt(s.Nf)
        s.Zmax=int(s.N.max())
    def chipow(s,k,j):
        e=s.e[:,k].astype(int); z=(e==255)
        v=1j**((j*np.where(z,0,e))%4); v=np.where(z,0,v); return v
    def series(s,k,j,kind='B'):
        """return list of (norms, coefs) pieces: B part and (for kind D) the extra c0*m^{j+1} part."""
        if j==0:
            return [(s.Nf, s.alpha*s.Gt)]
        (ma,mb),Nm,Gm=s.tp[k]
        cp=s.chipow(k,j)
        pieces=[(s.Nf, s.alpha*s.Gt*cp)]
        if kind=='D':
            m=(ma,mb); am=complex(ma,mb)/math.sqrt(Nm)
            fac=[m]
            ok=(s.e[:,k]!=255)
            if j==1:
                gq=gauss(m,2,(1,0),fac); extra=am**2*s.chipow(k,1)*s.Gt*gq
            elif j==2:
                g1=gauss(m,1,(1,0),fac); sgn=np.where((((s.N-1)//4)%2==1)&(((Nm-1)//4)%2==1),-1,1)
                extra=am**3*sgn*s.Gt*math.sqrt(Nm)*g1
            elif j==3:
                extra=am**4*(-Nm)*s.chipow(k,3)*s.Gt
            extra=np.where(ok,s.alpha*extra,0)
            pieces.append((s.Nf*Nm**(j+1), extra))
        return pieces
    def smooth(s,pieces,Zs):
        out=np.zeros(len(Zs),complex)
        for i,Z in enumerate(Zs):
            for (Nn,cf) in pieces:
                sel=Nn<=Z
                if sel.any(): out[i]+=np.sum(cf[sel]*W(Nn[sel]/Z))
        return out
    def sigma(s,Zs):
        return np.array([math.sqrt(np.sum(W(s.Nf[s.Nf<=Z]/Z)**2)) for Z in Zs])

def primary_elems(maxN):
    R=int(math.isqrt(maxN))+1
    return [(a,b) for a in range(-R,R+1) for b in range(-R,R+1) if 0<a*a+b*b<=maxN and is_primary((a,b))]

def zstar(s,pieces,Zs,twist=4):
    """S*(Z)=sum_d alpha(d)^twist N(d) S_D(Z/N(d)^4)"""
    ds=primary_elems(int(max(Zs)**0.25)+1)
    out=np.zeros(len(Zs),complex)
    for d in ds:
        Nd=gnorm(d); ad=complex(*d)/math.sqrt(Nd)
        Zd=np.array(Zs)/Nd**4
        if Zd.max()<0.5: continue
        out+= (ad**twist)*Nd*s.smooth(pieces,Zd)
    return out

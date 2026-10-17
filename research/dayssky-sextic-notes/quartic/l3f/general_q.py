# S52 bookkeeping at general c-length Y with Q = Y^q (q >= 2/3: c-sum at or below self-dual). log_Y units.
exec(open('fixa6.py').read().split("lam=13/8")[0])
def Xexp_q(P,eta,q):
    main=(3/4)*q+0.5+lpmin2(P,P[:,0]+eta*P[:,2],0.25,0.25)
    dual=0.5+q*0.5+lpmin2(P,P[:,0]+q*P[:,1]+eta*P[:,2],0.0,0.2499)
    nos=0.5+q*0.5+lpmin2(P,P[:,0]+q*P[:,1]+eta*P[:,2],0.2501,5.0)
    return min(nos,max(main,dual))
def book_q(P,q,n=21):
    H=2*q-1; etas=np.linspace(0,H,n); ex=[Xexp_q(P,e,q) for e in etas]; best=(9,None)
    for i,k0 in enumerate(etas):
        offT=max((1-q)+etas[j]+(H-etas[j])/4+2*ex[j] for j in range(i+1))
        s=max(2+(H-k0)/4,(1+offT)/2,2)
        if s<best[0]: best=(s,k0)
    return best
P=sample(1.5365,"all")
for q in (2/3,0.68,0.70,0.72,0.75,0.8):
    s,k0=book_q(P,q)
    # (L3)-type exponent for this box in Q-units: S/Y normalised by Y = Q^{1/q}: S_Q-exp = s/q ; compare with main (s0=2.0243 at q=2/3)
    print("q=%.3f (Y=Q^%.3f): S<<Y^%.4f = Q^%.4f   [main box: Q^%.4f]"%(q,1/q,s,s/q,2.0243*1.5))

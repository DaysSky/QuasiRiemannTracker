import numpy as np, time
A={}; exec(open('recur7.py').read().split('if __name__')[0],A)
B={}; exec(open('recur9.py').read().split('if __name__')[0],B)
rng=np.random.default_rng(1)
tab=np.array([A['S'](q,None) for q in A['QS']])-0.05*rng.random(len(A['QS']))
# sigma agreement where q in range
tC=rng.random(2000)*5; tD=rng.random(2000)*3
a=A['sigma_vec'](tC,tD,tab); b=B['sigma_vec'](tC,tD,tab)
inr=(tD/np.maximum(tC,1e-12)<=1.6)
print("sigma max diff (in range):",np.abs(a-b)[inr].max(), " out of range diffs (expected, b>=a):",(b-a)[~inr].min())
# b_twist agreement for a few points with tab=None
for (x,y,tK,tE,th,pol) in [(0.3,0.2,0.25,0.0,0.5,False),(0.45,0.1,0.3,0.02,0.7,True),(0.2,0.5,0.1,0.0,1.0,False)]:
    print(A['b_twist'](x,y,tK,tE,None,{},th=th,polar=pol), B['b_twist'](x,y,tK,tE,None,{},th=th,polar=pol))
for (x,y,pol) in [(0.3,0.2,False),(0.5,0.1,True),(0.4,0.6,False)]:
    print(A['okd'](x,y,None,pol), B['okd'](x,y,None,pol))


d=np.abs(a-b); i=np.argmax(d*inr); q=tD[i]/tC[i]
print("worst:",tC[i],tD[i],"q=",q,"qp=",q/(3*q-1) if q>1/3 else None,a[i],b[i])
qq=tD/np.maximum(tC,1e-12); qp=np.where((qq>1/3)&(qq<2/3),qq/np.maximum(3*qq-1,1e-12),qq)
ok=(qp<=1.6)
print("max diff where reflected q also in table range:",d[ok].max())

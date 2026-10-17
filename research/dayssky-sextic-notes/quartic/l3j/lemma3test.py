# Stress test of THEOREM_ZF_PROOF Lemma 3 (approximate positivity in a narrow annulus).
# For random complex a_c, u_c in [0,eta], r^2 in a range, B=2: compute
#   D = sum_{c1,c2} a1 conj(a2) G(u1+u2),  G(s)=(1+r^2 e^{-s})^{-B},
#   P_n = sum_c a_c (u_c/eta)^n,
# and report the worst ratio  (G(0)|P_0|^2 - 3D) / (B*eta*G(0)*sum_{n>=1} 2^{-n}|P_n|^2)  (should be bounded by a modest constant).
import numpy as np
rng=np.random.default_rng(1)
B=2.0
worst=0; worst_case=None
for trial in range(20000):
    eta=rng.choice([0.05,0.1,0.2])
    n=rng.integers(5,60)
    u=rng.uniform(0,eta,n)
    # adversarial-ish coefficients: random, or designed to make P_0 small
    a=rng.normal(size=n)+1j*rng.normal(size=n)
    if trial%3==0:
        a=a-np.mean(a)  # P_0 ~ 0
    if trial%3==1:
        a=a*np.exp(1j*rng.uniform(0,2*np.pi))*(u/eta-0.5)  # oscillating in u
    r2=10**rng.uniform(-3,3)
    G=lambda s:(1+r2*np.exp(-s))**(-B)
    S=u[:,None]+u[None,:]
    D=np.real(np.sum(a[:,None]*np.conj(a[None,:])*G(S)))
    P=[np.sum(a*(u/eta)**k) for k in range(0,12)]
    lhs=G(0)*abs(P[0])**2-3*D
    rhs=B*eta*G(0)*sum(2.0**(-k)*abs(P[k])**2 for k in range(1,12))
    if rhs>0:
        ratio=lhs/rhs
        if ratio>worst: worst=ratio; worst_case=(eta,n,r2,trial%3)
    elif lhs>1e-12:
        print("violation with rhs=0",lhs)
print("worst ratio (needs to be <= some constant C):",worst,"at (eta,n,r2,type)=",worst_case)
# also check D >= 0 fails sometimes (so positivity is genuinely approximate)
neg=0
for trial in range(5000):
    eta=0.2; n=30; u=rng.uniform(0,eta,n); a=rng.normal(size=n)+1j*rng.normal(size=n); a-=np.mean(a)
    r2=10**rng.uniform(-1,2); G=lambda s:(1+r2*np.exp(-s))**(-B)
    D=np.real(np.sum(a[:,None]*np.conj(a[None,:])*G(u[:,None]+u[None,:])))
    if D<0: neg+=1
print("fraction of random zero-mean cases with D<0:",neg/5000)

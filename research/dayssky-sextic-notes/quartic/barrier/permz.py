# per-test-prime z statistics for untwisted / alpha-twisted series at indices m, m^2, m^3
import sys; sys.argv=[sys.argv[0]]
exec(open('notheta.py').read().split('print("data')[0])
G=s.Gt
print(" m        N(m)   z(r=m)  z(r=m^2) z(r=m^3) | alpha: z(m) z(m^2) z(m^3) | c(r=m) c(r=m^3) [complex]")
for k in range(s.K):
    (ma,mb),Nm,_=s.tp[k]
    out=[];cc=[]
    for j in (0,1):
        for pw in (1,2,3):
            c,se,z,ka,mx=fit((s.alpha**j)*G*s.chipow(k,pw)); out.append(z); cc.append(c)
    print(f"{str((ma,mb)):10s} {Nm:5d}  " + " ".join(f"{v:7.2f}" for v in out[:3]) + " | " + " ".join(f"{v:6.2f}" for v in out[3:]) + f" | {cc[0]:.3f} {cc[2]:.3f}")

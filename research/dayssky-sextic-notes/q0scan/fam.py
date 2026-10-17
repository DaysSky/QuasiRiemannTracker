from fractions import Fraction as F
def rho(d): return F(1,6) + (F(5,6)-d)*5/(5-2*d)
def BL(q0):
    Ka = (2*q0/9)/(1-4*q0/3); r0 = rho(q0); c = Ka*(5-6*r0)
    return (1-r0+c)/(q0+c), (q0+F(9,2)*c)/(q0+c), Ka, r0
if __name__ == "__main__":
  for q in [F(65,100),F(66,100),F(67,100),F(17,25),F(69,100),F(70,100),F(71,100),F(72,100),F(73,100),F(74,100)]:
      B,L,Ka,r0 = BL(q)
      print(f"q0={float(q):.3f} B={float(B):.6f} Lam={float(L):.5f} Ka={float(Ka):.4f} rho0={float(r0):.5f}  cap@m=.417 (B-2m)/L={float((B-F(834,1000))/L):.5f}")

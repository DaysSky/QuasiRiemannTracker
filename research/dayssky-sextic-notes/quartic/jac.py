from gauss4 import *
import math
# J(chi_pi,chi_pi) = sum_{x mod pi} chi(x) chi(1-x)
for p in [5,13,17,29,37,41,53]:
    for pi in primary_prime_over(p):
        J=0
        for x in residues(pi):
            a=chi(x,pi,[pi]); b=chi((1-x[0],-x[1]),pi,[pi])
            J+=a*b
        print(p,pi,"J=",complex(round(J.real,6),round(J.imag,6)),"chi(-1)=",chi((-1,0),pi,[pi]), "pi mod4", (pi[0]%4,pi[1]%4))

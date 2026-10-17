# Halasz-Montgomery / Huxley large values for the sextic family (exponents base U).
# Witness: polynomial of length N=U^nu with |D|>=N^sigma (unnormalized, |b_n|<=1). Gram off-diagonal:
#   pointwise G = N^{1/2} * max|L(1/2, pair)|, pair conductor ~ U^2: Burgess-Wu exponent 3/16 per norm -> U^{3/8};
#   'lindelof' -> U^0.
# HM: count <= N^{2-2sig}  provided N^{2 sig} >= 2 N G.  Huxley with pair 4th moment: count <= N^{2-2sig} + U N^{3-4 sig}.
import numpy as np
def hm_ok(sig, nu, gexp):           # need 2 sig*nu >= nu + nu/2 + gexp
    return 2*sig*nu >= 1.5*nu + gexp
for sig in (0.70, 0.75, 0.80, 0.85):
    for gname, g in (('Burgess', 3/8), ('Lindelof', 0.0)):
        nus = [nu for nu in np.linspace(0.05, 6, 600) if hm_ok(sig, nu, g)]
        nmin = min(nus) if nus else None
        print(f"sigma={sig} {gname:8s}: HM valid for nu >= {nmin if nmin is None else round(nmin,3)}"
              + ("" if nmin is None else f"  -> count U^{round(nmin*(2-2*sig),3)} (DH {round(2-2*sig,3)}; need nu<=1 for DH)"))
    k = np.linspace(0.05, 6, 1200)                  # Huxley with power N^k
    hux = np.maximum(k*(2-2*sig), 1 + k*(3-4*sig))
    print(f"   Huxley(pair 4th moment) best count exponent {hux.min():.3f}  (DH {2-2*sig:.3f}, current rho {1-10*(2*sig-1)/(3*(5-2*(2*sig-1))):.3f})")

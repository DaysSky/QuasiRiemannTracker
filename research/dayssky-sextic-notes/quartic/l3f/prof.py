import sys; sys.argv=['x','all']
exec(open('hull3.py').read().split("which=sys.argv")[0])
P=sample(1.5417,"all")
H=1/3
for eta in np.linspace(0,H,9):
    X=Xexp(P,eta)
    need=(3-1/3-eta-(H-eta)/4)/2
    print("eta=%.3f  X=%.4f (main %.4f dual %.4f nos %.4f)  needed for s=2: %.4f  slope main=%.3f"%(eta,X[0],X[1],X[2],X[3],need,(X[1]-1)/eta if eta>0 else 0))

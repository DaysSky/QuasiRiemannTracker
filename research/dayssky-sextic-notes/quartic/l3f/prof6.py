exec(open('fixa6.py').read().split("lam=13/8")[0])
P=sample(1.5365,"all")
H=1/3
for eta in np.linspace(0,H,9):
    X=Xexp2(P,eta); need=(3-1/3-eta-(H-eta)/4)/2
    print("eta=%.3f X=%.4f main=%.4f dual=%.4f nos=%.4f need=%.4f"%(eta,X[0],X[1],X[2],X[3],need))

import sys, numpy as np, cmath
sys.path.insert(0,'..'); sys.path.insert(0,'.')
import gauss4 as g4
def load(fn):
    with open(fn,'rb') as f:
        NR,J=np.frombuffer(f.read(8),dtype=np.int32)
        rows=np.frombuffer(f.read(8*6*NR),dtype=np.float64).reshape(NR,6)
        win=np.frombuffer(f.read(),dtype=np.float64).reshape(J,NR,3)
    return NR,J,rows,win
NR,J,rows,win=load(sys.argv[1])
# brute-force check of row Gauss sums for small squarefree rows
cnt=0; maxerr=0
for i in range(NR):
    N,sq,a,b,gr,gi=rows[i]
    if sq!=1 or N>400 or N==1: continue
    n=(int(a),int(b))
    gb=g4.gauss(n,j=3)   # chi-bar Gauss sum
    G=complex(gr,gi)
    maxerr=max(maxerr, abs(G-gb))
    cnt+=1
    if cnt<=4: print(n, N, G/np.sqrt(N), gb/np.sqrt(N))
print("checked",cnt,"rows; min-convention max err",maxerr)

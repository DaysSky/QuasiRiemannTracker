import sys, j4test3 as J
for Z in [int(a) for a in sys.argv[1:]]:
    n,rt,rg=J.run(Z,Z,rows='all',weight='gauss',inert=True)
    print(f"Z=Q={Z}: observed S_T/diag-1={rt-1:.4f}  S_muT/diag-1={rg-1:.4f}",flush=True)

import capc4, sys
for cI in [float(a) for a in sys.argv[1:]]:
    print(f"c_I = {cI}: theta = {capc4.theta(4.498/0.7497, cI, 5/6):.5f}", flush=True)

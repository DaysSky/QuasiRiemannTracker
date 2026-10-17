import sys; sys.path.insert(0,'/tmp/claude-0/-home-user-QuasiRiemannTracker/726b9ed3-631a-533e-b4f5-2a9b6d3bc0e0/scratchpad/combo/gen')
from fmodel2 import model
Brad=142100527/143610000; L1=19312/4787
for dB in (0,2e-4,5e-4,1e-3,2e-3,4e-3):
    print("dB=%.4f"%dB, model(L1,Brad+dB))

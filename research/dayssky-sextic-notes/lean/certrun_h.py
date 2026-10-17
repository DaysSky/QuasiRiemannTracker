from model_h import *
import sys, time
FULL = (Q(0), Q(3,4), Q(0), Q(1,2))
def complement_boxes(full, reg):
    a,b,c,d = full; ra,rb,rc,rd = reg
    out = []
    if a < ra: out.append((a, ra, c, d))
    if rb < b: out.append((rb, b, c, d))
    if c < rc: out.append((ra, rb, c, rc))
    if rd < d: out.append((ra, rb, rd, d))
    return out
t = time.time()
res = {}
newrow = row_old(h, Q(1, 10**6), True)
print('deg new row', newrow.deg())
tot = []
for bx in complement_boxes(FULL, REG):
    tot += certify(newrow, bx, 2, 3)
print('newrow outside region boxes', len(tot), time.time()-t); res['newrow_out'] = tot
for name, pol, mdeg, ndeg in [('newrow_in', row_new(h, Q(1,10**6), True), 2, 3), ('V1', V1, 2, 3), ('V2', V2, 2, 3), ('V3', V3, 2, 3), ('V4', V4, 2, 3), ('V5', V5, 2, 3), ('V6', V6, 2, 4)]:
    print(name, 'deg', pol.deg())
    bs = certify(pol, REG, mdeg, ndeg); res[name] = bs
    print(name, 'boxes', len(bs), 'min', float(min(min(min(r) for r in co) for _, co in bs)))
for name, d, v, new in [('oldJunction', dm, Q(1, 10**5), False), ('oldHalf', Q(1,2), Q(1, 1000), False)]:
    bs = certify(row_old(d, v, new), FULL, 2, 3); res[name] = bs
    print(name, 'boxes', len(bs), 'min', float(min(min(min(r) for r in co) for _, co in bs)))
import pickle; pickle.dump(res, open('certs.pkl', 'wb'))
print('total time', time.time()-t)

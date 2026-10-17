"""Split B167Part1 (Cycle61 GeometryMainLowSelected) into three modules at declaration
boundaries. Statements and proofs are unchanged; every piece repeats the section's
namespace/open/variable preamble. The section's local instance is declared once (piece A)
and re-activated with `attribute [local instance]` in pieces B and C. B167Part1 becomes a
thin module importing the last piece, so B167Part2's import stays valid."""
import sys, os

d = sys.argv[1]
src = os.path.join(d, 'B167Part1.lean')
L = open(src, encoding='utf-8').read().split('\n')
idx = {i + 1: l for i, l in enumerate(L)}          # 1-based


def find(prefix):
    hits = [i for i, l in idx.items() if l.startswith(prefix)]
    assert len(hits) == 1, (prefix, hits)
    return hits[0]


begin = find('/- BEGIN ORIGINAL Cycle61.Physical.InitialMass.GeometryMainLowSelected')
sec = begin + 1
assert idx[sec] == 'section'
inst = find('local instance geometryMainLowSelectedDecidable')
cut1 = find('private theorem low_scalar_product_transport')
cut2 = find('theorem actual_geometry_low_paid_fields') - 2   # include its docstring
assert idx[cut2].startswith('/-- The payment dictionary'), idx[cut2]
close_nc = find('end OAI.SevenEighths.Cycle61InitialMass') - 1
assert idx[close_nc] == 'end'
endmark = find('/- END ORIGINAL Cycle61.Physical.InitialMass.GeometryMainLowSelected')
assert idx[endmark - 1] == 'end'

banner = [l for l in L[:begin - 1] if l.startswith('/-')][:1]   # packaging banner only
imports = [l for l in L[:begin - 1] if l.startswith('import ')]
pre = [idx[i] for i in range(sec, inst)]                         # section .. local notation
inst_line = idx[inst]
body = lambda a, b: [idx[i] for i in range(a, b)]
closer = ['end', 'end OAI.SevenEighths.Cycle61InitialMass', '', 'end']
tail = [idx[i] for i in range(close_nc, endmark + 1)]            # end/end ns/#print axioms/end/END

def note(k):
    return (f'/- B167 part 1{k} of 1A-1C: declarations of the original source section '
            f'Cycle61.Physical.InitialMass.GeometryMainLowSelected, unchanged, split for the '
            f'600-second per-module replay limit. -/')

A = banner + [note('A')] + imports + ['', idx[begin]] + pre + [inst_line] + body(inst + 1, cut1) + closer
B = (banner + [note('B')] + ['import Candidate.Bundles.B167Part1A', ''] + pre +
     ['attribute [local instance] geometryMainLowSelectedDecidable', ''] + body(cut1, cut2) + closer)
C = (banner + [note('C')] + ['import Candidate.Bundles.B167Part1B', ''] + pre +
     ['attribute [local instance] geometryMainLowSelectedDecidable', ''] + body(cut2, close_nc) + tail)

# every original body line appears exactly once across the pieces
orig = body(inst + 1, close_nc)
assert body(inst + 1, cut1) + body(cut1, cut2) + body(cut2, close_nc) == orig

for name, lines in (('B167Part1A', A), ('B167Part1B', B), ('B167Part1C', C)):
    open(os.path.join(d, name + '.lean'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
wrapper = banner + ['/- B167 part 1: the original section GeometryMainLowSelected is split into '
                    'B167Part1A, B167Part1B and B167Part1C. -/', 'import Candidate.Bundles.B167Part1C', '']
open(src, 'w', encoding='utf-8').write('\n'.join(wrapper))
print('cuts', inst, cut1, cut2, close_nc, 'pieces written')

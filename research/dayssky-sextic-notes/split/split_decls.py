"""Split one single-section Candidate module (BXXXPartK.lean) at declaration boundaries.

usage: split_decls.py DIR MODULE CUT_DECL [CUT_DECL ...]

Each CUT_DECL is the name of a top-level declaration that starts a new piece; the cut
is placed above its docstring/attribute lines (the block that ends at the declaration
line and is preceded by a blank line). Pieces MODULE{A,B,...} each repeat the section's
preamble (section ... local notation); a `local instance NAME` is declared in piece A
and re-activated by `attribute [local instance] NAME` in later pieces. Piece k imports
piece k-1; MODULE.lean becomes a thin module importing the last piece. Declarations,
statements and proofs are unchanged and appear exactly once, in the original order.
"""
import os
import re
import sys

d, mod, cuts = sys.argv[1], sys.argv[2], sys.argv[3:]
path = os.path.join(d, mod + '.lean')
L = open(path, encoding='utf-8').read().split('\n')


def one(pred, what):
    hits = [i for i, l in enumerate(L) if pred(l)]
    assert len(hits) == 1, (what, hits)
    return hits[0]


begin = one(lambda l: l.startswith('/- BEGIN ORIGINAL'), 'begin')
assert L[begin + 1] == 'section'
ns_end = one(lambda l: re.match(r'^end OAI\.', l) is not None, 'namespace end')
ns_name = L[ns_end][4:]
assert L[ns_end - 1] == 'end'
endmark = one(lambda l: l.startswith('/- END ORIGINAL'), 'end marker')
inst_hits = [i for i, l in enumerate(L) if l.startswith('local instance ')]
notation = max(i for i, l in enumerate(L[:ns_end]) if l.startswith('local notation') or l.startswith('variable'))
pre_end = max([notation] + inst_hits) + 1          # first body line
pre = [l for i, l in enumerate(L[begin + 1:pre_end], begin + 1) if i not in inst_hits]
inst_names = [re.match(r'local instance (\w+)', L[i]).group(1) for i in inst_hits]
for i in inst_hits:
    assert i < pre_end

starts = []
for name in cuts:
    k = one(lambda l: re.match(r'^(private |protected |noncomputable )*(theorem|lemma|def|abbrev) '
                               + re.escape(name) + r'\b', l) is not None, name)
    while L[k - 1].strip() != '':
        k -= 1
    starts.append(k)
assert starts == sorted(starts) and starts[0] > pre_end
bounds = [pre_end] + starts + [ns_end - 1]
pieces_body = [L[a:b] for a, b in zip(bounds, bounds[1:])]
assert sum(pieces_body, []) == L[pre_end:ns_end - 1]

banner = [l for l in L[:begin] if l.startswith('/-')][:1]
imports = [l for l in L[:begin] if l.startswith('import ')]
closer = ['end', 'end ' + ns_name, '', 'end']
tail = L[ns_end - 1:endmark + 1]
letters = 'ABCDEFGH'
n = len(pieces_body)
for k, body in enumerate(pieces_body):
    name = mod + letters[k]
    note = (f'/- {mod}{letters[k]} ({k + 1} of {n}): consecutive declarations of one original '
            f'source section, unchanged, split for the 600-second per-module replay limit. -/')
    head = banner + [note] + (imports if k == 0 else [f'import Candidate.Bundles.{mod}{letters[k - 1]}']) + ['']
    if k == 0:
        sec = [L[begin]] + L[begin + 1:pre_end]               # original preamble incl. instances
    else:
        sec = pre + [f'attribute [local instance] {x}' for x in inst_names] + ['']
    end = tail if k == n - 1 else closer
    open(os.path.join(d, name + '.lean'), 'w', encoding='utf-8').write('\n'.join(head + sec + body + end) + '\n')
wrapper = banner + [f'/- {mod}: its declarations are split into ' + ', '.join(mod + letters[k] for k in range(n)) + '. -/',
                    f'import Candidate.Bundles.{mod}{letters[n - 1]}', '']
open(path, 'w', encoding='utf-8').write('\n'.join(wrapper))
print(mod, 'split into', n, 'pieces at lines', [b + 1 for b in bounds])

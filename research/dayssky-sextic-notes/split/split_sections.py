"""Split a Candidate bundle into one module per original source section.

BXXXPart1 carries the bundle's original imports; BXXXPart(k) imports BXXXPart(k-1).
BXXX.lean becomes a thin module importing the last part, so its public name and
transitive content are unchanged. Usage: split_sections.py BUNDLE_DIR B167 [B176 ...]
"""
import os
import sys

bundle_dir = sys.argv[1]
for b in sys.argv[2:]:
    path = os.path.join(bundle_dir, b + '.lean')
    lines = open(path, encoding='utf-8').read().split('\n')
    i = 0
    head = []
    while not lines[i].startswith('/- BEGIN ORIGINAL'):
        head.append(lines[i])
        i += 1
    sections, cur = [], None
    for ln in lines[i:]:
        if ln.startswith('/- BEGIN ORIGINAL'):
            assert cur is None
            cur = [ln]
            continue
        if cur is None:
            assert ln.strip() == '', repr(ln)
            continue
        cur.append(ln)
        if ln.startswith('/- END ORIGINAL'):
            sections.append(cur)
            cur = None
    assert cur is None
    imports = [l for l in head if l.startswith('import ')]
    banner = [l for l in head if l.startswith('/-')]
    assert len(imports) + len(banner) == len([l for l in head if l.strip()]), head
    n = len(sections)
    rebuilt = '\n\n'.join('\n'.join(s) for s in sections)
    assert rebuilt.strip() == '\n'.join(lines[i:]).strip(), b
    for k, sec in enumerate(sections, 1):
        note = (f'/- {b} part {k} of {n}: one original source section, unchanged, split out '
                f'for the 600-second per-module replay limit. -/')
        imp = imports if k == 1 else [f'import Candidate.Bundles.{b}Part{k - 1}']
        text = '\n'.join(banner + [note] + imp + [''] + sec) + '\n'
        open(os.path.join(bundle_dir, f'{b}Part{k}.lean'), 'w', encoding='utf-8').write(text)
    wrapper = banner + [
        f'/- {b}: its {n} original source sections are {b}Part1 ... {b}Part{n}, '
        f'chained in their original order. -/',
        f'import Candidate.Bundles.{b}Part{n}', '']
    open(path, 'w', encoding='utf-8').write('\n'.join(wrapper))
    print(b, n, 'sections, bodies preserved')

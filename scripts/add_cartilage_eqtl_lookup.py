# -*- coding: utf-8 -*-
"""Fix: insert the cartilage-eQTL Results paragraph (saved first), then renumber."""
import sys, io, re, copy
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
from docx import Document
from docx.text.paragraph import Paragraph

F = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002\submission_PeerJ_v2_20261121\01_Manuscript.docx"

PARA = ('Because the expression instruments derive from blood, we also asked whether VDR is genetically '
        'regulated in cartilage. In the largest available cartilage cis-eQTL resource, generated from 204 '
        'osteoarthritic cartilage samples[38], VDR was not a significant eGene: its lead variant rs34149349 '
        '(chr12:48,200,006, hg19) had a nominal P of 5.7 \u00d7 10\u207b\u2074 but a permutation P of 0.17 '
        'and a permutation-based q value of 0.27, against a threshold of q < 0.05 applied to 18,756 tested '
        'genes, of which 3,109 were significant. The blood-derived instruments therefore have no validated '
        'cartilage counterpart.')

# ---------- insert if absent ----------
d = Document(F)
P = list(d.paragraphs)
if any('rs34149349' in p.text for p in P):
    print('paragraph already present')
else:
    head = next(i for i, p in enumerate(P) if p.text.strip() == 'Tissue-level validation in independent bulk OA cartilage')
    body = next(i for i in range(head + 1, len(P)) if P[i].text.strip())
    el = copy.deepcopy(P[body]._p)
    P[body]._p.addnext(el)
    para = Paragraph(el, P[body]._parent)
    para.runs[0].text = PARA
    for r in para.runs[1:]:
        r.text = ''
    d.save(F)
    print('paragraph inserted and saved after paragraph', body)

# ---------- renumber ----------
d = Document(F)
P = list(d.paragraphs)
ref_head = next(i for i, p in enumerate(P) if p.text.strip() == 'References')
supp_head = next(i for i, p in enumerate(P) if p.text.strip() == 'Supporting information')
REFS = {}
for i in range(ref_head + 1, supp_head):
    t = P[i].text.strip()
    if not t:
        continue
    m = re.match(r'^(\d+)\.\s*(.+)$', t, flags=re.S)
    assert m, 'bad ref line %d' % i
    REFS[int(m.group(1))] = m.group(2).strip()
assert sorted(REFS) == list(range(1, 39)), sorted(REFS)
print('refs:', len(REFS))

s = next(i for i, p in enumerate(P) if p.text.strip().lower() == 'introduction')
NUM = re.compile(r'\[([0-9,\s\u2013-]+)\]')


def expand(chunk):
    out = []
    for part in re.split(r'[,\s]+', chunk.strip()):
        if not part:
            continue
        if re.fullmatch(r'\d+', part):
            out.append(int(part))
        elif re.fullmatch(r'\d+[-\u2013]\d+', part):
            a, b = re.split(r'[-\u2013]', part)
            out.extend(range(int(a), int(b) + 1))
    return out


first, seq = {}, []
for i in range(s, ref_head):
    for m in NUM.finditer(P[i].text):
        for n in expand(m.group(1)):
            if n not in first:
                first[n] = (i, m.start())
                seq.append((first[n], n))
seq.sort(key=lambda x: x[0])
order = [n for _, n in seq]
print('cited numbers:', len(order))
missing = [n for n in REFS if n not in order]
assert not missing, 'uncited references: %s' % missing
mapping = {old: new for new, old in enumerate(order, start=1)}

for i in range(s, ref_head):
    text = P[i].text
    if '[' not in text:
        continue
    nt = text
    for m in reversed(list(NUM.finditer(text))):
        nums = expand(m.group(1))
        if not nums:
            continue
        mapped = sorted({mapping.get(n, n) for n in nums})
        parts, k = [], 0
        while k < len(mapped):
            j = k
            while j + 1 < len(mapped) and mapped[j + 1] == mapped[j] + 1:
                j += 1
            if j - k >= 2:
                parts.append('%d-%d' % (mapped[k], mapped[j]))
            else:
                parts.extend(str(x) for x in mapped[k:j + 1])
            k = j + 1
        rep = '[%s]' % ', '.join(parts)
        if rep != m.group(0):
            nt = nt[:m.start()] + rep + nt[m.end():]
    if nt != text:
        P[i].runs[0].text = nt
        for r in P[i].runs[1:]:
            r.text = ''

new_list = [REFS[old] for old in order]
for i, txt in enumerate(new_list):
    P[ref_head + 1 + i].runs[0].text = '%d.\t%s' % (i + 1, txt)
    for r in P[ref_head + 1 + i].runs[1:]:
        r.text = ''
d.save(F)
print('renumbered; cartilage ref is now [%d]' % mapping[38])

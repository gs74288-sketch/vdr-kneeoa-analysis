# -*- coding: utf-8 -*-
"""Repair the supplementary tables after the bad rename:
   S6 <- pseudotime (from the 11-07 backup)
   S7 <- GSE324993 detection (rebuilt from CSV)
   S8 <- scDRS (currently sitting in S6_Table.docx with a mangled caption)
"""
import sys, io, os, shutil
import pandas as pd
from docx import Document

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

BASE = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002\submission_PeerJ_v2_20261121\Supplementary"
OLD7 = r"C:\Users\liuha\Desktop\孟德尔\revision_20261002\submission_PeerJ_20261107\Supplementary\S7_Table.docx"
TAB = r"C:\Users\liuha\OpenClawVDR\repo\tables"


def set_caption(path, caption):
    d = Document(path)
    first = next(p for p in d.paragraphs if p.text.strip())
    first.runs[0].text = caption
    for r in first.runs[1:]:
        r.text = ''
    d.save(path)
    return caption


# ---- 1. current S6_Table.docx actually holds the scDRS table -> becomes S8 ----
cur = os.path.join(BASE, 'S6_Table.docx')
d = Document(cur)
print('current S6 caption was:', next(p.text for p in d.paragraphs if p.text.strip())[:120])
print('  tables inside:', [(len(t.rows), len(t.columns)) for t in d.tables])

CAP8 = ('Table S8. Single-cell disease relevance scores for knee osteoarthritis in GSE324993. '
        'Upper panel: cells, mean and median normalised score by chondrocyte state. Lower panel: '
        'numbers of regulatory (RegC) and other cells with the within-donor difference in score, '
        'by donor and group.')
set_caption(cur, CAP8)
os.replace(cur, os.path.join(BASE, 'S8_Table.docx'))
print('  -> S8_Table.docx (scDRS)')

# ---- 2. pseudotime table from the 11-07 backup -> S6 ----
shutil.copyfile(OLD7, os.path.join(BASE, 'S6_Table.docx'))
set_caption(os.path.join(BASE, 'S6_Table.docx'),
            'Table S6. Pseudotime root-node sensitivity using donor-state pseudobulk models.')
print('  -> S6_Table.docx (pseudotime, restored)')

# ---- 3. rebuild the GSE324993 detection table -> S7 ----
by_state = pd.read_csv(os.path.join(TAB, 'GSE324993_VDR_by_state.csv'))
per_donor = pd.read_csv(os.path.join(TAB, 'GSE324993_VDR_per_donor.csv'))

doc = Document()
CAP = ('Table S7. VDR detection by chondrocyte state and donor in the GSE324993 atlas. '
       'Upper panel: cells analysed, cells with detectable VDR transcript and the detection '
       'fraction, by annotated state. Lower panel: regulatory (RegC) and other cells with their '
       'detection fractions and the within-donor difference, by donor and group.')
doc.add_paragraph(CAP)

doc.add_paragraph('By chondrocyte state')
t1 = doc.add_table(rows=1, cols=4)
t1.style = 'Table Grid'
hdr = t1.rows[0].cells
for c, v in zip(hdr, ['State', 'Cells', 'VDR detected', 'Detection fraction']):
    c.text = v
for _, r in by_state.iterrows():
    cells = t1.add_row().cells
    cells[0].text = str(r['state'])
    cells[1].text = '%d' % r['n']
    cells[2].text = '%d' % r['detected']
    cells[3].text = '%.3f' % r['frac']

doc.add_paragraph('')
doc.add_paragraph('By donor')
t2 = doc.add_table(rows=1, cols=7)
t2.style = 'Table Grid'
hdr = t2.rows[0].cells
for c, v in zip(hdr, ['Donor', 'Group', 'RegC cells', 'RegC fraction',
                      'Other cells', 'Other fraction', 'Difference (pp)']):
    c.text = v
for _, r in per_donor.iterrows():
    cells = t2.add_row().cells
    cells[0].text = str(r['sample'])
    cells[1].text = str(r['group'])
    cells[2].text = '%d' % r['reg_n']
    cells[3].text = '%.3f' % r['reg_frac']
    cells[4].text = '%d' % r['other_n']
    cells[5].text = '%.3f' % r['other_frac']
    cells[6].text = '%+.1f' % (r['delta'] * 100)

doc.save(os.path.join(BASE, 'S7_Table.docx'))
print('  -> S7_Table.docx (GSE324993, rebuilt)')

print('\nfinal:')
for f in sorted(os.listdir(BASE)):
    print('  ', f)

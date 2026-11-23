#!/usr/bin/env python
"""Convert the manuscript Markdown files to DOCX (no pandoc required)."""
import re, sys, os
from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

INLINE = re.compile(r"(\*\*[^*]+\*\*|\*[^*]+\*|\^[^\^]+\^)")

def add_runs(par, text, size=None, bold=False, italic=False):
    for tok in INLINE.split(text):
        if not tok:
            continue
        b, i, sup = bold, italic, False
        if tok.startswith("**") and tok.endswith("**") and len(tok) > 4:
            tok = tok[2:-2]; b = True
        elif tok.startswith("*") and tok.endswith("*") and len(tok) > 2:
            tok = tok[1:-1]; i = True
        elif tok.startswith("^") and tok.endswith("^") and len(tok) > 2:
            tok = tok[1:-1]; sup = True
        r = par.add_run(tok)
        r.bold = b; r.italic = i
        if sup:
            r.font.superscript = True
        if size:
            r.font.size = Pt(size)

def set_cell_bg(cell, color):
    tcPr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd"); shd.set(qn("w:val"), "clear")
    shd.set(qn("w:color"), "auto"); shd.set(qn("w:fill"), color)
    tcPr.append(shd)

def add_table(doc, rows, font=7):
    if not rows:
        return
    ncol = max(len(r) for r in rows)
    t = doc.add_table(rows=len(rows), cols=ncol)
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for ri, row in enumerate(rows):
        for ci in range(ncol):
            cell = t.cell(ri, ci)
            cell.text = ""
            p = cell.paragraphs[0]
            txt = row[ci] if ci < len(row) else ""
            add_runs(p, txt, size=font)
            if ri == 0:
                for r in p.runs:
                    r.bold = True
                set_cell_bg(cell, "D9E2F3")
            for r in p.runs:
                r.font.size = Pt(font)
    doc.add_paragraph()


def __style_headings(doc):
    """Force black bold Times New Roman on heading styles (Word defaults are blue Calibri Light)."""
    from docx.shared import RGBColor
    for lvl, size in ((1, 14), (2, 13), (3, 11.5)):
        st = doc.styles["Heading %d" % lvl]
        st.font.name = "Times New Roman"
        st.font.size = Pt(size)
        st.font.bold = True
        st.font.color.rgb = RGBColor(0, 0, 0)
        rpr = st.element.get_or_add_rPr()
        rf = rpr.find(qn("w:rFonts"))
        if rf is None:
            rf = OxmlElement("w:rFonts"); rpr.append(rf)
        for a in ("w:ascii", "w:hAnsi", "w:eastAsia", "w:cs"):
            rf.set(qn(a), "Times New Roman")

def convert(md_path, docx_path, base_size=11):
    lines = open(md_path, encoding="utf-8").read().split("\n")
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(base_size)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "SimSun")
    __style_headings(doc)
    for s in doc.sections:
        s.top_margin = Cm(2.0); s.bottom_margin = Cm(2.0)
        s.left_margin = Cm(2.2); s.right_margin = Cm(2.2)
    i = 0
    while i < len(lines):
        line = lines[i].rstrip()
        s = line.strip()
        if not s:
            i += 1
            continue
        if s.startswith("|"):
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-{2,}:?", c) for c in cells if c):
                    rows.append(cells)
                i += 1
            add_table(doc, rows)
            continue
        if s.startswith("# "):
            h = doc.add_heading(level=1)
            add_runs(h, s[2:], size=15)
            for r in h.runs:
                r.bold = True
            i += 1
            continue
        if s.startswith("## "):
            h = doc.add_heading(level=2)
            add_runs(h, s[3:], size=13)
            for r in h.runs:
                r.bold = True
            i += 1
            continue
        if s.startswith("### "):
            h = doc.add_heading(level=3)
            add_runs(h, s[4:], size=11.5)
            for r in h.runs:
                r.bold = True
            i += 1
            continue
        if s.startswith("> "):
            p = doc.add_paragraph()
            add_runs(p, s[2:])
            for r in p.runs:
                r.italic = True
            i += 1
            continue
        m = re.match(r"^(\d+)\.\s+(.*)$", s)
        if m:
            p = doc.add_paragraph(style="List Number")
            add_runs(p, m.group(2), size=base_size - 1)
            i += 1
            continue
        if s.startswith("- "):
            p = doc.add_paragraph(style="List Bullet")
            add_runs(p, s[2:])
            i += 1
            continue
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        add_runs(p, s)
        i += 1
    doc.save(docx_path)
    return docx_path

if __name__ == "__main__":
    md, out = sys.argv[1], sys.argv[2]
    size = int(sys.argv[3]) if len(sys.argv) > 3 else 11
    convert(md, out, size)
    print("wrote", out)

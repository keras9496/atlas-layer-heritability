"""원고 v2 Word 파일 만들기. 본문은 manuscript_main.md, 보충은 manuscript_supp.md에 있다(이 스크립트에는 본문이 없다).
 표기: # / ## 제목, %키 값 표지 정보, ::figure 파일|라벨|제목 … ::end, ::table 라벨|제목[|landscape] (다음 줄 설명, | 표 행) … ::end
       (landscape: 보충 자료에서 그 표를 가로 페이지에 둔다. 넓은 표용, 2026-10-17)
       **굵게**, *기울임*, ^위첨자^, ~아래첨자~, \\* 글자 그대로의 별표, [@key; @key] 인용 (본문 → 보충 순 첫 등장 번호)
 서식: 1e v4와 같은 Candara 12 pt, 줄 번호, 쪽 번호.
출력: 01_Main_Text_v2.docx (제출 형식: 그림·표 자리 표시, 끝에 참고문헌 → Figure legends → Tables)
      02_Supplementary_v2.docx, 03_Figures_v2.pptx, review_copy_v2.docx (그림·표를 본문 안에 넣은 검토용)
"""
import json, re, html
from pathlib import Path
from docx import Document
from docx.enum.section import WD_ORIENT, WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor
P = Path(__file__).resolve().parent; FIG = P/"figures"; FONT = "Candara"
REFS = json.loads((P/"refs.json").read_text(encoding="utf-8"))

# ------------------------------------------------------------ 원고 파싱
def parse(path):
    meta, blocks, lines = {}, [], (P/path).read_text(encoding="utf-8").splitlines(); i = 0
    while i < len(lines):
        s = lines[i].rstrip()
        if s.startswith("%"): k, _, v = s[1:].partition(" "); meta.setdefault(k, []).append(v); i += 1
        elif s.startswith("## "): blocks.append(("h2", s[3:])); i += 1
        elif s.startswith("# "): blocks.append(("h1", s[2:])); i += 1
        elif s.startswith("::figure "):
            f, label, title = s[9:].split("|"); body = []; i += 1
            while lines[i].strip() != "::end": body.append(lines[i].strip()); i += 1
            blocks.append(("fig", f, label, title, " ".join(body))); i += 1
        elif s.startswith("::table "):
            parts = s[8:].split("|"); label, title = parts[:2]; opt = parts[2] if len(parts) > 2 else ""; legend, rows = [], []; i += 1
            while lines[i].strip() != "::end":
                t = lines[i].strip()
                if t.startswith("|"): rows.append([c.strip() for c in t.strip("|").split("|")])
                elif t: legend.append(t)
                i += 1
            blocks.append(("table", label, title, " ".join(legend), rows, opt)); i += 1
        elif not s.strip(): i += 1
        else:
            buf = []
            while i < len(lines) and lines[i].strip() and not lines[i].startswith(("#", "::", "%")): buf.append(lines[i].strip()); i += 1
            blocks.append(("p", " ".join(buf)))
    return meta, blocks
MAIN_META, MAIN = parse("manuscript_main.md"); SUPP_META, SUPP = parse("manuscript_supp.md")

# ------------------------------------------------------------ 인용 번호 (본문 등장 순 → 보충)
ORDER = []
CIT = re.compile(r"\[(@[^\]]+)\]")
def texts(blocks):
    for b in blocks:
        for x in b[1:]:
            if isinstance(x, str): yield x
            elif isinstance(x, list): yield from (c for r in x for c in r)
for t in list(texts(MAIN)) + list(texts(SUPP)):
    for m in CIT.finditer(t):
        for k in [k.strip().lstrip("@") for k in m.group(1).split(";")]:
            assert k in REFS, f"참고문헌 없음: {k}"
            if k not in ORDER: ORDER.append(k)
def cite(m):
    n = sorted(ORDER.index(k.strip().lstrip("@")) + 1 for k in m.group(1).split(";")); out, i = [], 0
    while i < len(n):
        j = i
        while j + 1 < len(n) and n[j + 1] == n[j] + 1: j += 1
        out.append(f"{n[i]}–{n[j]}" if j - i >= 2 else ",".join(map(str, n[i:j + 1]))); i = j + 1
    return "[" + ",".join(out) + "]"
GROUP_FIRST = {"roadmap2015"}                       # 논문 표기가 "단체, 개인 저자…"인 경우
GROUP_ONLY = {"kg2015", "gtex2020", "encode2020", "pgc2011"}    # 논문 표기가 단체 이름인 경우
def fmt_ref(r, key):
    g, pp = r.get("groups", []), r["authors"]
    au = g[:1] if key in GROUP_ONLY or not pp else (g[:1] + pp if key in GROUP_FIRST else pp)
    a = ", ".join(au[:6]) + (" et al." if len(au) > 6 else "")
    title = re.sub(r"<[^>]+>", "", html.unescape(r["title"])).strip().rstrip(".")
    title = re.sub(r"(\w)\s+-\s*(\w)", r"\1-\2", re.sub(r"\s+", " ", title))   # "disease <b>-</b> risk" → "disease-risk"
    s = f"{a}{'' if a.endswith('.') else '.'} {title}{'' if title.endswith(('?', '!')) else '.'} *{r['journal']}*. {r['year']}"
    if r.get("volume"): s += f";{r['volume']}"
    if r.get("pages"): s += ":" + r["pages"].replace("-", "–")
    s += "." + (f" doi:{r['doi']}" if r.get("doi") else (f" ({r['note']})" if r.get("note") else ""))
    return s

# ------------------------------------------------------------ Word 서식
def new_doc(line_numbers=True):
    doc = Document(); st = doc.styles["Normal"]; st.font.name = FONT; st.font.size = Pt(12)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), FONT); st.paragraph_format.space_after = Pt(6); st.paragraph_format.line_spacing = 1.5
    for lv, sz in [(1, 14), (2, 12)]:
        h = doc.styles[f"Heading {lv}"]; h.font.name = FONT; h.font.size = Pt(sz); h.font.bold = True; h.font.italic = False; h.font.color.rgb = RGBColor(0, 0, 0)
        h.element.rPr.rFonts.set(qn("w:eastAsia"), FONT); h.element.rPr.rFonts.set(qn("w:ascii"), FONT); h.element.rPr.rFonts.set(qn("w:hAnsi"), FONT)
        h.paragraph_format.space_before = Pt(14 if lv == 1 else 10); h.paragraph_format.space_after = Pt(4)
    sec = doc.sections[0]
    for side in ["left_margin", "right_margin", "top_margin", "bottom_margin"]: setattr(sec, side, Inches(1))
    if line_numbers:
        ln = OxmlElement("w:lnNumType"); ln.set(qn("w:countBy"), "1"); ln.set(qn("w:restart"), "continuous"); ln.set(qn("w:distance"), "283"); sec._sectPr.append(ln)
    fp = sec.footer.paragraphs[0]; fp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for t, v in [("begin", None), (None, "PAGE"), ("end", None)]:
        r = fp.add_run()
        if t: e = OxmlElement("w:fldChar"); e.set(qn("w:fldCharType"), t)
        else: e = OxmlElement("w:instrText"); e.set(qn("xml:space"), "preserve"); e.text = v
        r._r.append(e)
    return doc
TOK = re.compile(r"(\*\*.+?\*\*|\*.+?\*|\^.+?\^|~.+?~)")
STAR = ""
def runs(par, text, size=None, bold=False):
    text = CIT.sub(cite, text).replace("\\*", STAR)
    for part in TOK.split(text):
        if not part: continue
        if part.startswith("**"): r = par.add_run(part[2:-2]); r.bold = True
        elif part.startswith("*") and len(part) > 1: r = par.add_run(part[1:-1]); r.italic = True
        elif part.startswith("^"): r = par.add_run(part[1:-1]); r.font.superscript = True
        elif part.startswith("~"): r = par.add_run(part[1:-1]); r.font.subscript = True
        else: r = par.add_run(part)
        r.text = r.text.replace(STAR, "*"); r.font.name = FONT
        if bold: r.bold = True
        if size: r.font.size = Pt(size)
    return par
def para(doc, text, align=None, size=None, bold=False, after=None, spacing=None):
    p = doc.add_paragraph(); runs(p, text, size, bold)
    if align == "c": p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    if align == "j": p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    if after is not None: p.paragraph_format.space_after = Pt(after)
    if spacing: p.paragraph_format.line_spacing = spacing
    return p
def heading(doc, text, lv):
    h = doc.add_heading(level=lv); runs(h, text); return h
def pagebreak(doc): doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)
def border(cell, **kw):
    tcPr = cell._tc.get_or_add_tcPr(); b = tcPr.find(qn("w:tcBorders"))
    if b is None: b = OxmlElement("w:tcBorders"); tcPr.append(b)
    for edge, val in kw.items():
        e = OxmlElement(f"w:{edge}"); e.set(qn("w:val"), val); e.set(qn("w:sz"), "8"); e.set(qn("w:color"), "000000"); b.append(e)
def render_table(doc, label, title, legend, rows, size=8.5, width=6.5):
    p = para(doc, f"**{label} | {title}** {legend}", size=10, align="j", after=4, spacing=1.15); p.paragraph_format.keep_with_next = True
    ncol = len(rows[0])
    def vis(s): return len(re.sub(r"[\*\^~\\]", "", CIT.sub("[00]", s)))
    L = [max(vis(r[j]) if j < len(r) else 0 for r in rows) for j in range(ncol)]
    cw = .0068 * size   # 글자 하나의 대략 너비(인치)
    nat = [cw * (l + 2) + .16 for l in L]; mn = [cw * (min(l, (20 if width <= 6.5 else 14) if j == 0 else 10) + 2) + .16 for j, l in enumerate(L)]
    if sum(nat) <= width: widths = [x * width / sum(nat) for x in nat]
    else: f = (width - sum(mn)) / (sum(nat) - sum(mn)); widths = [a + (b - a) * max(f, 0) for a, b in zip(mn, nat)]
    t = doc.add_table(rows=len(rows), cols=ncol); t.alignment = WD_TABLE_ALIGNMENT.CENTER; t.autofit = False
    for i, row in enumerate(rows):
        for j in range(ncol):
            c = t.cell(i, j); c.width = Inches(widths[j]); pp = c.paragraphs[0]; pp.paragraph_format.space_after = Pt(1); pp.paragraph_format.line_spacing = 1.0
            runs(pp, row[j] if j < len(row) else "", size, bold=(i == 0))
            if i == 0: border(c, top="single", bottom="single")
            if i == len(rows) - 1: border(c, bottom="single")
    doc.add_paragraph().paragraph_format.space_after = Pt(4)
def legend_par(doc, label, title, legend, size=11):
    return para(doc, f"**{label} | {title}** {legend}", align="j", after=10, size=size, spacing=1.15)
def picture(doc, f, width=6.5):
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER; p.paragraph_format.keep_with_next = True; p.add_run().add_picture(str(FIG/f), width=Inches(width))

def title_block(doc, meta, supp=False):
    if supp: para(doc, "**Supplementary Information**", align="c", size=16, after=6); para(doc, "for", align="c", size=11, after=6)
    para(doc, f"**{meta['title'][0]}**", align="c", size=12 if supp else 16, after=14, spacing=1.15)
    para(doc, meta["authors"][0], align="c", after=8)
    for a in meta.get("affil", []): para(doc, a, align="c", size=10, after=2, spacing=1.15)
    para(doc, meta["corr"][0], align="c", size=10, after=14)
    if "keywords" in meta: para(doc, "**Keywords:** " + meta["keywords"][0], size=10)
    if "contents" in meta: para(doc, "**Contents:** " + meta["contents"][0], size=10)

def build_main(inline):
    doc = new_doc(); title_block(doc, MAIN_META); pagebreak(doc); legends, tables = [], []
    for b in MAIN:
        if b[0] == "h1":
            if b[1] in ("Data availability",) and not inline: pass
            heading(doc, b[1], 1)
        elif b[0] == "h2": heading(doc, b[1], 2)
        elif b[0] == "p": para(doc, b[1], align="j")
        elif b[0] == "fig":
            _, f, label, title, legend = b
            if inline: picture(doc, f); legend_par(doc, label, title, legend, size=10)
            else: para(doc, f"**[Insert {label} here]**", align="c", after=12); legends.append((label, title, legend))
        elif b[0] == "table":
            _, label, title, legend, rows, _opt = b
            if inline: render_table(doc, label, title, legend, rows)
            else: para(doc, f"**[Insert {label} here]**", align="c", after=12); tables.append((label, title, legend, rows))
    pagebreak(doc); heading(doc, "References", 1)
    for i, k in enumerate(ORDER, 1):
        p = para(doc, f"{i}. " + fmt_ref(REFS[k], k), size=10, after=2, spacing=1.15); p.paragraph_format.left_indent = Inches(.3); p.paragraph_format.first_line_indent = Inches(-.3)
    if not inline:
        pagebreak(doc); heading(doc, "Figure legends", 1)
        for l in legends: legend_par(doc, *l)
        pagebreak(doc); heading(doc, "Tables", 1)
        for t in tables: render_table(doc, *t)
    return doc
def orient(doc, landscape):
    """새 쪽에서 시작하는 구역. 쪽 번호·줄 번호 설정은 앞 구역에서 이어진다."""
    s = doc.add_section(WD_SECTION.NEW_PAGE); w, h = s.page_width, s.page_height
    if landscape != (w > h): s.page_width, s.page_height = h, w
    s.orientation = WD_ORIENT.LANDSCAPE if landscape else WD_ORIENT.PORTRAIT
def build_supp():
    doc = new_doc(); title_block(doc, SUPP_META, supp=True)
    for k, b in enumerate(SUPP):
        if b[0] == "h1":
            nxt = SUPP[k + 1] if k + 1 < len(SUPP) else None
            if nxt and nxt[0] == "table" and nxt[5] == "landscape": orient(doc, True)
            else: pagebreak(doc)
            heading(doc, b[1], 1)
        elif b[0] == "h2": heading(doc, b[1], 2)
        elif b[0] == "p": para(doc, b[1], align="j")
        elif b[0] == "table" and b[5] == "landscape":
            if not (k > 0 and SUPP[k - 1][0] == "h1"): orient(doc, True)
            render_table(doc, *b[1:5], size=7.5, width=9.0); orient(doc, False)
        elif b[0] == "table": render_table(doc, *b[1:5], size=8)
        elif b[0] == "fig": picture(doc, b[1]); legend_par(doc, b[2], b[3], b[4], size=10)
    return doc

build_main(False).save(P/"01_Main_Text_v2.docx")
build_main(True).save(P/"review_copy_v2.docx")
build_supp().save(P/"02_Supplementary_v2.docx")
# 그림 PPTX (한 슬라이드에 한 그림)
from pptx import Presentation
from pptx.util import Inches as PI, Pt as PPt
prs = Presentation(); prs.slide_width, prs.slide_height = PI(13.333), PI(7.5)
for b in MAIN + SUPP:
    if b[0] != "fig": continue
    s = prs.slides.add_slide(prs.slide_layouts[6]); tb = s.shapes.add_textbox(PI(.4), PI(.15), PI(12.5), PI(.5)); tb.text_frame.text = f"{b[2]} | {b[3]}"
    tb.text_frame.paragraphs[0].runs[0].font.size = PPt(16); tb.text_frame.paragraphs[0].runs[0].font.bold = True
    from PIL import Image
    w, h = Image.open(FIG/b[1]).size; maxw, maxh = 12.5, 6.6; sc = min(maxw / (w / 300), maxh / (h / 300)); ww, hh = w / 300 * sc, h / 300 * sc
    s.shapes.add_picture(str(FIG/b[1]), PI((13.333 - ww) / 2), PI(.75 + (6.6 - hh) / 2), PI(ww), PI(hh))
prs.save(P/"03_Figures_v2.pptx")
words = sum(len(re.sub(r"\[@[^\]]+\]", "", b[1]).split()) for b in MAIN if b[0] == "p")
abstract = next(len(MAIN[i + 1][1].split()) for i, b in enumerate(MAIN) if b == ("h1", "Abstract"))
unused = [k for k in REFS if k not in ORDER]
print(f"참고문헌 {len(ORDER)}개 인용" + (f", 미인용 {unused}" if unused else "") + f" | 본문 문단 단어 수 약 {words} (초록 {abstract})")

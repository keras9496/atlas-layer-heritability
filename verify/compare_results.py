"""Compare the result files of two runs of the pipeline: those in this copy of the repository (the reference) and those
in another copy (the reproduction; for example one made with tools/make_run_copy.py and run with its own data root).

    python verify/compare_results.py <other copy or run folder> [--rtol 1e-6] [--out verify/report.md] [--data-root D:/data]

For every file under experiments/*/results/ in this copy (the reference, called the manuscript's copy below), the trait table that e20
writes next to its results folder, and the figure files under figures/out/ (PNG files, the values drawn and the table
of Fig. 5b; when the other folder is a run copy), the same file is looked up in the other folder. CSV files are compared column by column: numeric columns within a relative tolerance
(and an absolute tolerance of 1e-9 near zero), other columns exactly, after sorting rows by the common columns.
JSON files are compared entry by entry: the same structure, keys and strings, numbers within the same tolerances.
Text files (.txt, .md) are compared line by line after masking what legitimately differs between two runs
(timestamps, dates, elapsed times and the data-root path, which is recognised when it is named AD_GWAS_Atlas_*, is given
with --data-root or is the value of AGA_DATA); a text file whose only differences are of that kind is
reported as "match (timing/path lines only)". The dated notes RESULTS_<tag>_<date>.md that the summary scripts write
are paired by tag and order; the manuscript's copies carry an interpretation section added by hand after the run, so
a reproduced note whose lines are all contained in the manuscript's note is reported as "match (manuscript copy has
added notes)". Notes with no writer in the pipeline and outputs of scripts that are not part of this repository are
listed under their own status. Files that exist only on one side are listed. The report is written as Markdown.
The exit status is 0 when every file is a match of some kind.
"""
import os, sys, re, json, argparse, difflib
import numpy as np, pandas as pd
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ap = argparse.ArgumentParser(); ap.add_argument("other"); ap.add_argument("--rtol", type=float, default=1e-6); ap.add_argument("--atol", type=float, default=1e-9)
ap.add_argument("--out", default=str(REPO/"verify"/"report.md")); ap.add_argument("--data-root", action="append", default=[]); a = ap.parse_args()
ROOTS = sorted({r.replace("\\", "/").rstrip("/") for r in a.data_root + [os.environ.get("AGA_DATA", ""), os.environ.get("AGA_DATA_POSIX", "")] if r}, key=len, reverse=True)   # data roots to mask in logs
OTHER = Path(a.other).resolve(); OTHER_ROOT = OTHER if (OTHER/"experiments").exists() else None
if (OTHER/"experiments").exists(): OTHER = OTHER/"experiments"
REF = REPO/"experiments"
EXTRA = {"e20_gwas_panel/traits_exp3u.csv"}                      # written by a script outside results/ (e20 u00_inputs.py)
OK = {"same bytes", "match", "match (timing/path lines only)", "match (manuscript copy has added notes)",
      "note written by hand (not a pipeline output)", "output of a script not included (AlphaGenome API)", "differ (expected, no value differs)"}
# result files the pipeline cannot produce (none: the output of the excluded z01_lineages.py is kept in e02_lineage_map/inputs)
NOT_INCLUDED = set()
# text logs whose manuscript copy is known to differ from a clean run for a reason that does not touch any number
EXPECTED = {"e01_targets_sldsc/results/x06_prepare_3b.txt": "the manuscript's copy is the 3-line log of a partial re-run of x06 (only its last block); a full run writes all blocks, with the same values",
            "e13_larger_gwas/results/n01_munge.txt": "the manuscript's copy repeats the LDSC munge summary on the same line (the summary was appended twice by the original run); the values are identical"}
# RESULTS_<tag>_<date>.md tags that a script of this repository writes (the others are notes written by hand)
WRITTEN_TAGS = set()
for s in REF.glob("*/scripts/*.py"):
    WRITTEN_TAGS |= set(re.findall(r"RESULTS_(EXP[0-9a-z]+)_\{", s.read_text(encoding="utf-8", errors="replace")))
DATED = re.compile(r"^RESULTS_(EXP[0-9A-Za-z_]+?)_(\d{4}-\d{2}-\d{2})\.md$")

def norm(s):
    """Mask what differs between two runs of the same code: data root, timestamps, dates, elapsed times."""
    for r in ROOTS: s = s.replace(r, "<DATA>").replace(r.replace("/", "\\"), "<DATA>")
    s = re.sub(r"[A-Za-z]:[\\/]AD_GWAS_Atlas_\w+", "<DATA>", s); s = re.sub(r"/[a-z]/AD_GWAS_Atlas_\w+", "<DATA>", s)
    s = re.sub(r"\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}(:\d{2})?", "<TS>", s); s = re.sub(r"\d{4}-\d{2}-\d{2}", "<DATE>", s)
    s = re.sub(r"\(\s*[\d.]+\s*(분|초|s|sec|min|h)\s*\)", "(<DUR>)", s); s = re.sub(r"[\d.]+\s*(분|초)\b", "<DUR>", s)
    s = re.sub(r"\b\d+:\d{2}:\d{2}\b", "<DUR>", s); s = re.sub(r"(elapsed|소요|경과|took)[^\n]*", r"\1 <DUR>", s)
    return s

def compare_csv(f1, f2):
    if f1.read_bytes() == f2.read_bytes(): return "same bytes", ""         # also covers empty CSVs (e.g. a diffs file with no rows)
    try: A = pd.read_csv(f1); B = pd.read_csv(f2)
    except Exception as e: return "error", f"could not read: {e}"
    if list(A.columns) != list(B.columns): return "columns differ", f"{len(A.columns)} vs {len(B.columns)} columns"
    if len(A) != len(B): return "rows differ", f"{len(A)} vs {len(B)} rows"
    key = [c for c in A.columns if A[c].dtype == object]
    if key: A = A.sort_values(key, kind="stable").reset_index(drop=True); B = B.sort_values(key, kind="stable").reset_index(drop=True)
    bad = []
    for c in A.columns:
        x, y = A[c], B[c]
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            x = x.astype(float).to_numpy(); y = y.astype(float).to_numpy(); nan = np.isnan(x) & np.isnan(y)
            ok = nan | np.isclose(x, y, rtol=a.rtol, atol=a.atol, equal_nan=True)
            if not ok.all():
                with np.errstate(divide="ignore", invalid="ignore"): rel = np.nanmax(np.abs(x - y)[~ok] / np.maximum(np.abs(y[~ok]), a.atol))
                bad.append(f"{c}: {(~ok).sum()} of {len(x)} values differ (max rel {rel:.2e})")
        else:
            ne = ~((x.astype(str) == y.astype(str)) | (x.isna() & y.isna()))
            if ne.any(): bad.append(f"{c}: {ne.sum()} values differ")
    return ("match", f"{len(A)} rows x {len(A.columns)} columns") if not bad else ("differ", "; ".join(bad))

def compare_json(f1, f2):
    """Same structure, keys and strings; numbers within the tolerance (a text comparison would report a difference in the last digit of one value)."""
    if f1.read_bytes() == f2.read_bytes(): return "same bytes", ""
    try: x = json.loads(f1.read_text(encoding="utf-8")); y = json.loads(f2.read_text(encoding="utf-8"))
    except Exception as e: return "error", f"could not read: {e}"
    n = [0]; bad = []; worst = [0.0]
    def walk(u, v, path):
        if isinstance(u, dict) and isinstance(v, dict):
            if list(u) != list(v): bad.append(f"{path or '/'}: keys differ"); return
            for k in u: walk(u[k], v[k], f"{path}/{k}")
        elif isinstance(u, list) and isinstance(v, list):
            if len(u) != len(v): bad.append(f"{path or '/'}: {len(u)} vs {len(v)} items"); return
            for k, (p, q) in enumerate(zip(u, v)): walk(p, q, f"{path}[{k}]")
        elif isinstance(u, (int, float)) and isinstance(v, (int, float)) and not isinstance(u, bool) and not isinstance(v, bool):
            n[0] += 1
            if not np.isclose(u, v, rtol=a.rtol, atol=a.atol, equal_nan=True):
                bad.append(f"{path}: {u!r} vs {v!r}"); worst[0] = max(worst[0], abs(u - v) / max(abs(v), a.atol))
        elif u != v: bad.append(f"{path or '/'}: {u!r} vs {v!r}")
    walk(x, y, "")
    return ("match", f"{n[0]} numbers") if not bad else ("differ", f"{len(bad)} of {n[0]} numbers or entries differ (max rel {worst[0]:.2e}); first: {bad[0][:120]}")

def compare_text(f1, f2, dated=False):
    t1 = f1.read_text(encoding="utf-8", errors="replace").splitlines(); t2 = f2.read_text(encoding="utf-8", errors="replace").splitlines()
    if t1 == t2: return "match", f"{len(t1)} lines"
    n1 = [norm(l) for l in t1]; n2 = [norm(l) for l in t2]
    if n1 == n2: return "match (timing/path lines only)", f"{len(t1)} lines"
    d = [l for l in difflib.unified_diff(n1, n2, lineterm="", n=0) if l[:1] in "+-" and l[:3] not in ("+++", "---")]
    placeholder = {"## 해석", "(작성 예정)", "(결과 확인 후 작성)", ""}                 # the summary scripts end the note with an empty interpretation section
    if dated and all(l[0] == "-" or l[1:].strip() in placeholder for l in d):   # that the manuscript's copy replaces with the written interpretation
        return "match (manuscript copy has added notes)", f"{len(t2)} reproduced lines all present in the manuscript's {len(t1)}"
    return "differ", f"{len(d)} changed lines of {len(t1)} (after masking timestamps/durations/paths)"

rows = []; used = set()
for f in sorted(REF.rglob("*")):
    if not f.is_file() or f.name == ".gitkeep" or ("results" not in f.relative_to(REF).parts and f.relative_to(REF).as_posix() not in EXTRA): continue
    rel = f.relative_to(REF); key = rel.as_posix(); g = OTHER/rel; m = DATED.match(f.name)
    if key in NOT_INCLUDED: rows.append((key, "output of a script not included (AlphaGenome API)", "")); continue
    if m:                                                        # dated note: pair by tag and order with the other side's notes of the same tag
        tag = m.group(1)
        if tag not in WRITTEN_TAGS: rows.append((key, "note written by hand (not a pipeline output)", "")); continue
        mine = sorted(x for x in f.parent.glob(f"RESULTS_{tag}_*.md") if DATED.match(x.name) and DATED.match(x.name).group(1) == tag)
        theirs = sorted(x for x in (OTHER/rel.parent).glob(f"RESULTS_{tag}_*.md") if DATED.match(x.name) and DATED.match(x.name).group(1) == tag)
        i = mine.index(f)
        if i >= len(theirs): rows.append((key, "missing in reproduction", f"no {i + 1}-th note with tag {tag}")); continue
        g = theirs[i]; used.add(g.relative_to(OTHER).as_posix()); s, d = compare_text(f, g, dated=True); rows.append((key, s, f"paired with {g.name}; {d}")); continue
    if not g.exists(): rows.append((key, "missing in reproduction", "")); continue
    used.add(key)
    if f.suffix == ".csv": s, d = compare_csv(f, g)
    elif f.suffix == ".json": s, d = compare_json(f, g)
    elif f.suffix in (".txt", ".md", ".tsv"): s, d = compare_text(f, g)
    else: s, d = ("same bytes" if f.read_bytes() == g.read_bytes() else "bytes differ"), ""
    if s == "differ" and key in EXPECTED: s, d = "differ (expected, no value differs)", EXPECTED[key] + f"; {d}"
    rows.append((key, s, d))
for g in sorted(OTHER.rglob("*")):
    k = g.relative_to(OTHER).as_posix()
    if g.is_file() and g.name != ".gitkeep" and "results" in g.relative_to(OTHER).parts and k not in used and not (REF/k).exists(): rows.append((k, "only in reproduction", ""))
if OTHER_ROOT is not None and (REPO/"figures"/"out").exists():       # figures: PNG files byte by byte, the values drawn in them and the table of Fig. 5b as text or CSV
    for f in sorted((REPO/"figures"/"out").iterdir()):
        if not f.is_file() or f.suffix == ".pdf" or f.name == ".gitkeep": continue
        key = "figures/out/" + f.name; g = OTHER_ROOT/"figures"/"out"/f.name
        if not g.exists(): rows.append((key, "missing in reproduction", "")); continue
        if f.suffix == ".csv": s, d = compare_csv(f, g)
        elif f.suffix == ".json": s, d = compare_json(f, g)
        else: s, d = ("same bytes" if f.read_bytes() == g.read_bytes() else "bytes differ"), ""
        rows.append((key, s, d))

R = pd.DataFrame(rows, columns=["file", "status", "detail"])
counts = R.status.value_counts().to_dict()
def short(p):
    try: return p.relative_to(REPO.parent).as_posix()
    except ValueError: return p.as_posix()
lines = [f"# Reproduction check", "", f"Reference: `{short(REF)}`  ", f"Reproduction: `{short(OTHER)}`  ", f"Tolerance: rtol {a.rtol:g}, atol {a.atol:g}", "",
         "| Status | Files |", "|---|---|"] + [f"| {k} | {v} |" for k, v in counts.items()] + ["", "| File | Status | Detail |", "|---|---|---|"]
lines += [f"| {r.file} | {r.status} | {r.detail} |" for r in R.itertuples()]
Path(a.out).parent.mkdir(parents=True, exist_ok=True); Path(a.out).write_text("\n".join(lines) + "\n", encoding="utf-8")
print("\n".join(f"{k:>48s}: {v}" for k, v in counts.items())); print("report:", a.out)
sys.exit(0 if set(counts) <= OK else 1)

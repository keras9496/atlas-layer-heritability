"""Build the public repository from the working folders (2026-10-30).

Copy-only: reads the curated 'version 2' copy and the manuscript folder, writes only inside this repository.
Original folders are never modified. Re-running overwrites the copies (results/ and scripts/ are rebuilt from the sources).

What it does
  1. Copies, for each of the 18 experiments, the dated protocol, the scripts used in the paper and the result files,
     into experiments/<english name>/ (mapping in RENAME below).
  2. Rewrites paths so that the scripts run from any clone:
       "D:/AD_GWAS_Atlas_data"  -> f"{DATA}" (config.DATA, env AGA_DATA)       in .py
       /d/AD_GWAS_Atlas_data     -> ${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}      in .sh
       sibling experiment folders in Korean -> their English names
       the absolute path of the working folder -> the repository's experiments/ folder
  3. Copies the manuscript build files into paper/ and rewrites make_figures.py in the same way.
  4. Writes tools/rewrite_log.tsv (every substitution, per file) so the changes can be audited.
Usage: python tools/build_repo.py
"""
import re, shutil, csv
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
V2 = Path("C:/Users/User/alphagenome/Atlas_AD_Genomewide 실험 (1-3)/version 2/experiments")
PAPER_SRC = Path("C:/Users/User/alphagenome/Atlas_AD_Genomewide 실험 (1-3)/paper/manuscript_v1")
LOCI_SRC = Path("C:/Users/User/alphagenome/Atlas_AD_GWAS (1e)/results/upstream")

RENAME = {  # working folder -> repository folder (order = order of the analysis stages)
    "실험3_SLDSC": "e01_targets_sldsc", "실험3c_형질세포지도": "e02_lineage_map", "실험3d_모달리티지도": "e03_single_layer_map",
    "실험3e_실측주석조건부": "e04_measured_annotations", "실험3f_층위이질성": "e05_layer_map", "실험3g_역치민감도": "e06_thresholds",
    "실험3h_계통고정비교": "e07_shared_lineages", "실험3i_세포형일치": "e08_cell_type_matched", "실험3j_형질간공분산": "e09_joint_jackknife",
    "실험3k_근위원위": "e10_proximal_distal", "실험3l_활성도조건부": "e11_activity", "실험3m_기존산출물점검": "e12_diagnostics",
    "실험3n_최신GWAS": "e13_larger_gwas", "실험3o_연속구간": "e14_continuous_bins", "실험3p_활성도구성": "e15_activity_constructions",
    "실험3q_활성요소내효과": "e16_within_active", "실험3r_더큰AD_GWAS": "e17_larger_ad_gwas", "실험3s_세포형eQTL조건화": "e18_cell_type_eqtl"}
EXCLUDE_SCRIPTS = {"v06_figure_validation.py", "p05_fig_draft.py"}   # working figures not used in the paper
PAPER_FILES = ["make_figures.py", "build_manuscript.py", "fetch_refs.py", "export_pdf.ps1", "refs_keys.tsv", "refs.json", "numbers.json",
               "manuscript_main.md", "manuscript_supp.md"]

CONFIG_IMPORT = ("import sys as _sys; from pathlib import Path as _Path\n"
                 "_sys.path.insert(0, str(_Path(__file__).resolve().parents[{n}])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)\n")
LOG = []

def rewrite_py(text, depth, fname):
    """depth: number of parents from the file to the repository root (3 for experiments/*/scripts/*.py, 1 for paper/*.py)."""
    out = text; subs = []
    def rep(pat, new, label, regex=False):
        nonlocal out
        n = len(re.findall(pat, out)) if regex else out.count(pat)
        if n: out = re.sub(pat, new, out) if regex else out.replace(pat, new); subs.append((label, n))
    rep('Path("D:/AD_GWAS_Atlas_data")', 'Path(DATA)', "data root")
    rep('"D:/AD_GWAS_Atlas_data/', 'f"{DATA}/', "data path (double quotes)")
    rep("'D:/AD_GWAS_Atlas_data/", "f'{DATA}/", "data path (single quotes)")
    exp_dir = 'Path(__file__).resolve().parents[2]' if depth == 3 else 'Path(__file__).resolve().parents[1]/"experiments"'   # experiments/*/scripts/x.py -> parents[2] is experiments/; paper/x.py -> repo/experiments
    rep('Path("C:/Users/User/alphagenome/AD_GWAS_Atlas_연구")', exp_dir, "working folder -> experiments/")
    rep('Path("C:/Users/User/alphagenome/Atlas_AD_GWAS (1e)/results/upstream")', 'Path(__file__).resolve().parents[1]/"inputs"', "locus files -> inputs/")
    rep('R0/"실험1e_AVI기전"/"results"/"e13_loci_s3"/"loci75_table_S3.csv"', 'Path(__file__).resolve().parents[1]/"inputs"/"loci75_table_S3.csv"', "1e locus table -> inputs/")
    for ko, en in RENAME.items(): rep(ko, en, f"folder {ko}")
    rep("D:/AD_GWAS_Atlas_data", "<DATA>", "data root mentioned in comments or docstrings")   # only prose can remain after the rules above
    if "DATA" in out and "from config import DATA" not in out:
        # insert after the module docstring if there is one, else at the top
        m = re.match(r'\s*("""[\s\S]*?"""|\'\'\'[\s\S]*?\'\'\')\s*\n', out)
        imp = CONFIG_IMPORT.format(n=depth)
        out = out[:m.end()] + imp + out[m.end():] if m else imp + out
        subs.append(("config import", 1))
    for label, n in subs: LOG.append((fname, label, n))
    assert "D:/AD_GWAS_Atlas_data" not in out and "C:/Users" not in out, fname
    return out

def rewrite_sh(text, fname):
    out = text; subs = []
    for pat, new, label in [("/d/AD_GWAS_Atlas_data", "${AGA_DATA_POSIX:-/d/AD_GWAS_Atlas_data}", "data root (posix)"),
                            ("D:/AD_GWAS_Atlas_data", "${AGA_DATA:-D:/AD_GWAS_Atlas_data}", "data root")]:
        n = out.count(pat)
        if n: out = out.replace(pat, new); subs.append((label, n))
    for ko, en in RENAME.items():
        n = out.count(ko)
        if n: out = out.replace(ko, en); subs.append((f"folder {ko}", n))
    for label, n in subs: LOG.append((fname, label, n))
    return out

EXCLUDE_RESULTS = {"p05_fig_draft.png"}   # output of the excluded draft-figure script
def copy_tree_filtered(src, dst):
    """results/: everything except stdout redirects (*.log; the *_main.txt etc. files hold the same text) and EXCLUDE_RESULTS."""
    if dst.exists(): shutil.rmtree(dst)
    dst.mkdir(parents=True)
    for f in sorted(src.rglob("*")):
        if f.is_file() and "__pycache__" not in f.parts and f.suffix != ".log" and f.name not in EXCLUDE_RESULTS:
            d = dst/f.relative_to(src); d.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(f, d)

# 1. experiments
EXP = REPO/"experiments"
for ko, en in RENAME.items():
    src, dst = V2/ko, EXP/en
    for sub in ("scripts", "results"):
        if (dst/sub).exists(): shutil.rmtree(dst/sub)
    dst.mkdir(parents=True, exist_ok=True)
    for f in src.glob("PROTOCOL*"): shutil.copy2(f, dst/f.name)
    copy_tree_filtered(src/"results", dst/"results")
    (dst/"scripts").mkdir()
    for f in sorted((src/"scripts").iterdir()):
        if f.name in EXCLUDE_SCRIPTS or f.name.endswith(".pyc"): continue
        rel = f"experiments/{en}/scripts/{f.name}"
        if f.suffix == ".py": (dst/"scripts"/f.name).write_text(rewrite_py(f.read_text(encoding="utf-8"), 3, rel), encoding="utf-8")
        elif f.suffix == ".sh": (dst/"scripts"/f.name).write_text(rewrite_sh(f.read_text(encoding="utf-8"), rel), encoding="utf-8", newline="\n")
        else: shutil.copy2(f, dst/"scripts"/f.name)
# the pre-stage AD summary-statistics step (run interactively on 2026-09-25, later written down as code/00_prepare_sumstats.py)
pre = Path("C:/Users/User/alphagenome/Atlas_AD_Genomewide 실험 (1-3)/code/00_prepare_sumstats.py")
(EXP/"e01_targets_sldsc"/"scripts"/"x00_prepare_sumstats.py").write_text(
    rewrite_py(pre.read_text(encoding="utf-8"), 3, "experiments/e01_targets_sldsc/scripts/x00_prepare_sumstats.py"), encoding="utf-8")
# track metadata and lineage assignment (output of the excluded z01_lineages.py; needed by e02 z03, e04 w01, e11 l02)
inp2 = EXP/"e02_lineage_map"/"inputs"; inp2.mkdir(exist_ok=True); shutil.copy2("D:/AD_GWAS_Atlas_data/exp3c/track_lineage.csv", inp2/"track_lineage.csv")
# locus lists used by 3b (lead variants of public AD GWAS; no individual data)
inp = EXP/"e01_targets_sldsc"/"inputs"; inp.mkdir(exist_ok=True)
for n in ("loci.csv", "rep25_loci.csv"): shutil.copy2(LOCI_SRC/n, inp/n)
# 75-locus table of the 1e study (Supplementary Table S3 there; lead variants, GWAS Catalog match, nearest genes; used by e12 m04a for block names)
inp3 = EXP/"e12_diagnostics"/"inputs"; inp3.mkdir(exist_ok=True)
shutil.copy2("C:/Users/User/alphagenome/AD_GWAS_Atlas_연구/실험1e_AVI기전/results/e13_loci_s3/loci75_table_S3.csv", inp3/"loci75_table_S3.csv")

# 2. paper
P = REPO/"paper"; P.mkdir(exist_ok=True)
for n in PAPER_FILES:
    s = PAPER_SRC/n
    if n == "make_figures.py": (P/n).write_text(rewrite_py(s.read_text(encoding="utf-8"), 1, f"paper/{n}"), encoding="utf-8")
    else: shutil.copy2(s, P/n)
(P/"figures").mkdir(exist_ok=True)
for f in sorted((PAPER_SRC/"figures").glob("*.png")): shutil.copy2(f, P/"figures"/f.name)

# 3. log
with open(REPO/"tools"/"rewrite_log.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["file", "substitution", "count"]); w.writerows(LOG)
print(f"experiments: {len(RENAME)}; substitutions logged: {len(LOG)} in {len(set(l[0] for l in LOG))} files")

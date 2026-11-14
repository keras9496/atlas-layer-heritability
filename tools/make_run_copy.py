"""Make a second copy of the code in which to run the pipeline, for example to compare two runs
(verify/compare_results.py).

    python tools/make_run_copy.py [target]        default target: <repository>/repro_run

Copies scripts, inputs, config.py, run_all.sh and figures/, without anything the pipeline writes (result files, figures/out
and the trait table of e20); creates empty results/ folders. The default target is ignored by git (.gitignore).
Afterwards:  AGA_DATA=... bash <target>/run_all.sh     and     python verify/compare_results.py <target>
"""
import shutil, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DST = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else REPO/"repro_run"
if DST.exists(): raise SystemExit(f"{DST} exists; remove it first")
SKIP_TOP = {".git", "repro_run", "verify", "tools", "__pycache__"}
WRITTEN = {"experiments/e20_gwas_panel/traits_exp3u.csv"}                                # written by e20 u00_inputs.py outside results/
for item in sorted(REPO.iterdir()):
    if item.name in SKIP_TOP: continue
    if item.is_file(): shutil.copy2(item, DST/item.name) if DST.mkdir(parents=True, exist_ok=True) is None else None; continue
    for f in sorted(item.rglob("*")):
        if not f.is_file() or "__pycache__" in f.parts: continue
        rel = f.relative_to(REPO)
        if "results" in rel.parts: continue                                              # start with empty results folders
        if rel.parts[0] == "figures" and rel.parts[1] == "out": continue                 # written by the figure scripts
        if rel.as_posix() in WRITTEN: continue
        d = DST/rel; d.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(f, d)
for e in (REPO/"experiments").iterdir():
    if e.is_dir(): (DST/"experiments"/e.name/"results").mkdir(parents=True, exist_ok=True)
(DST/"figures"/"out").mkdir(parents=True, exist_ok=True)
n = sum(1 for f in DST.rglob("*") if f.is_file()); print(f"run copy at {DST}: {n} files")

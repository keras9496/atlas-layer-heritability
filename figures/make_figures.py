# -*- coding: utf-8 -*-
"""Draw all figures (fig1-fig5, figS1-figS5) from the result files under experiments/*/results.

    python figures/make_figures.py

Each figure script can also be run on its own (python figures/fig3.py). Outputs go to figures/out/: <name>.png, <name>.pdf and
numbers_<name>.json (the values drawn in the figure); fig5.py also writes eqtl14.csv (the values of its panel b).
"""
import subprocess, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
for name in ["fig1", "fig2", "fig3", "fig4", "fig5", "figS1", "figS2", "figS3", "figS4", "figS5"]:
    print(f"--- {name}", flush=True)
    subprocess.run([sys.executable, "-X", "utf8", str(HERE / f"{name}.py")], check=True)

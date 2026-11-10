"""실험 3p p04: 원고 그림용 값 내보내기 (p03과 같은 τ* 계산, 새 S-LDSC 없음).
 모형 E, EA, EAc, S, SA, EAM(원래 GWAS 14형질)과 E, EA, EAc, S, SA(새 GWAS 5형질)의 종류 × 층위 τ*, 단일 형질 jackknife SE, 대비 c1·c2와 SE.
출력: results/p04_values.csv (trait, gwas, model, kind, layer, tau, se, p1), results/p04_contrasts.csv (trait, gwas, model, c1, c1_se, c2, c2_se)
"""
import sys
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
sys.path.insert(0, str(Path(__file__).resolve().parent))
from p_common import *
EXP = Path(__file__).resolve().parents[1]
rows, crow = [], []
for gwas, traits, models in [("original", ALL, MODELS), ("new", NEW, ["E", "EA", "EAc", "S", "SA"])]:
    for t in traits:
        for model in models:
            R = load(t, model, gwas)
            for kind, d in R.items():
                for m, (tau, dv) in d.items(): se = jse(dv); rows.append(dict(trait=t, gwas=gwas, model=model, kind=kind, layer=m, tau=float(tau), se=se, p1=float(stats.norm.sf(tau / se))))
            c, dv = prof(R["E"]); crow.append(dict(trait=t, gwas=gwas, model=model, c1=float(c[0]), c1_se=jse(dv[:, 0]), c2=float(c[1]), c2_se=jse(dv[:, 1])))
V = pd.DataFrame(rows); C = pd.DataFrame(crow)
assert abs(V.query("trait == 'AD' and gwas == 'original' and model == 'E' and kind == 'E' and layer == 'CAGE'").tau.iloc[0] - 2.91) < .006
V.to_csv(EXP/"results"/"p04_values.csv", index=False); C.to_csv(EXP/"results"/"p04_contrasts.csv", index=False)
print(f"p04: 값 {len(V)}행, 대비 {len(C)}행 저장")

"""실험 3m m06: 원고 보충 표 S14용 값 내보내기 (2026-10-20, 기존 실행만 읽음, 새 S-LDSC 없음).
 AD 대 키 제외 집합(ADH_k5, ADH_k10)과 제외 없음(3j T3all)에서 AD·키의 층위별 τ*와 c1 = CAGE − DNase (단일 형질 jackknife SE).
출력: results/m06_decomp.csv
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
from m_common import load, jse, RUNS, LAY
R3M = Path(f"{DATA}/exp3m/runs"); OUT = Path(__file__).resolve().parents[1]/"results"/"m06_decomp.csv"
rows = []
for k, (model, runs) in {"none": ("T3all", RUNS), "k5": ("T3all__ADH_k5", R3M), "k10": ("T3all__ADH_k10", R3M)}.items():
    for t in ["AD", "HEIGHT"]:
        L = load(t, model, runs); nb = len(L["CAGE"][1])
        for m in LAY: rows.append(dict(k=k, trait=t, term=m, est=L[m][0], se=jse(L[m][1]), nblocks=nb))
        c = (L["CAGE"][0] - L["DNASE"][0], L["CAGE"][1] - L["DNASE"][1]); rows.append(dict(k=k, trait=t, term="c1", est=c[0], se=jse(c[1]), nblocks=nb))
D = pd.DataFrame(rows); D["p1"] = stats.norm.sf(D.est / D.se); D["p2"] = 2 * stats.norm.sf(abs(D.est / D.se))
D.to_csv(OUT, index=False); print(D.round(3).to_string(index=False))

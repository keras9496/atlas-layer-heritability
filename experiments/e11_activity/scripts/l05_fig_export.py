"""실험 3l l05: 원고 그림 5b·c와 보충 표 S15용 값 내보내기 (2026-10-20, 기존 실행만 읽음, 새 S-LDSC 없음).
 l04와 같은 τ* 계산. 모형 E(= 3j T3all), EA, A, EAM에서 종류(E 효과 / A 활성도 / M 실측) × 층위의 τ*와 단일 형질 jackknife SE.
출력: results/l05_values.csv (긴 형식: trait, model, kind, layer, tau, se, p1)
"""
import sys as _sys; from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parents[3])); from config import DATA   # repository root; DATA = data root (env AGA_DATA)
import re
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
EXP = Path(__file__).resolve().parents[1]
Y = Path(f"{DATA}/exp3d"); W = Path(f"{DATA}/exp3e"); J = Path(f"{DATA}/exp3j"); L = Path(f"{DATA}/exp3l"); X = Path(f"{DATA}/exp3")
LAY = ["CAGE", "DNASE", "H3K27ac"]
SD = pd.concat([pd.read_csv(p, index_col=0).sd_common for p in [Y/"annot_sd_common.csv", W/"measured_sd_common.csv", L/"annot_sd_common.csv"]]); SD = SD[~SD.index.duplicated()]
B0 = X/"GRCh38"/"baselineLD_v2.2"; M_TOT = sum(float(open(B0/f"baselineLD.{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
KIND = {"L5_": "E", "A5_": "A", "M_": "M"}
ALL = ["AD", "MONO", "EA", "SCZ", "BIP", "NEURO", "BMI", "IBD", "RA", "LDL", "HEIGHT", "BMD", "NEUT", "TC"]
strip = lambda n: re.sub(r"L2_\d+$", "", n)
def jse(v): v = np.asarray(v); B = len(v); return float(np.sqrt((B - 1) / B * ((v - v.mean()) ** 2).sum()))
rows = []
for t in ALL:
    for model, runs in [("E", (J/"runs", "T3all")), ("EA", (L/"runs", "EA")), ("A", (L/"runs", "A")), ("EAM", (L/"runs", "EAM"))]:
        z = np.load(runs[0]/f"{t}__{runs[1]}.npz"); h2 = float(z["h2_log"])
        for k, n in enumerate(strip(x) for x in z["names"]):
            if n not in SD.index: continue
            kind = [v for p, v in KIND.items() if n.startswith(p)][0]; m = [x for x in LAY if n.endswith("_" + x)][0]; f = SD[n] * M_TOT / h2
            tau, se = float(z["coef"][k] * f), jse(z["part_delete"][:, k] * f)
            rows.append(dict(trait=t, model=model, kind=kind, layer=m, tau=tau, se=se, p1=float(stats.norm.sf(tau / se))))
D = pd.DataFrame(rows)
assert abs(D.query("trait == 'AD' and model == 'E' and layer == 'CAGE'").tau.iloc[0] - 2.91) < .006   # l04 기준값
assert abs(D.query("trait == 'AD' and model == 'EA' and kind == 'A' and layer == 'CAGE'").tau.iloc[0] - 1.91) < .006
D.to_csv(EXP/"results"/"l05_values.csv", index=False)
print(D.pivot_table(index=["trait"], columns=["model", "kind", "layer"], values="tau").round(2).to_string())

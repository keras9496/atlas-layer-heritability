"""참고문헌 메타데이터를 Crossref에서 받아 refs.json으로 쓴다 (DOI 확인 겸).
 입력: refs_keys.tsv (key, doi). 수동 항목(DOI 없음)은 MANUAL에 둔다.
 출력: refs.json, 확인 기록 refs_check.txt
"""
import json, time, urllib.request, urllib.parse
from pathlib import Path
P = Path(__file__).resolve().parent
# choi2026: (2026-10-20 사용자 확인: medRxiv 저자 ByungHyug Choi, Jaeseung Jeong — 아래 값 그대로) Crossref 메타데이터는 제1저자의 성과 이름이 뒤바뀌어 있어(given "Choi", family "byung hyug") 수동으로 둔다. 제목·DOI는 medRxiv에서 확인(2026-10-12).
MANUAL = {"choi2026": dict(authors=["Choi B", "Jeong J"], title="Expression predictions inform molecular interpretation of Alzheimer's disease risk loci",
                           journal="medRxiv", year="2026", volume="", pages="", doi="10.64898/2026.10.07.26364928", note="preprint")}
rows = [l.split("\t") for l in (P/"refs_keys.tsv").read_text(encoding="utf-8").splitlines()[1:] if l.strip()]
out, log = {}, []
for key, doi in rows:
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi.strip())
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "manuscript-refcheck/1.0"})
        m = json.load(urllib.request.urlopen(req, timeout=30))["message"]
    except Exception as e:
        log.append(f"{key}\t{doi}\t실패 {e}"); continue
    au, groups = [], []   # 개인 저자와 단체(컨소시엄) 저자를 따로 둔다. 표기 순서는 build_manuscript.py가 정한다.
    for a in m.get("author", []):
        if "family" in a: au.append(a["family"] + " " + "".join(w[0].upper() for w in a.get("given", "").replace("-", " ").split() if w))
        elif "name" in a and not a["name"].rstrip(":").lower().endswith(("authors", "author")): groups.append(a["name"])
    j = (m.get("container-title") or [""])[0] or (m.get("institution", [{}])[0].get("name", "") if m.get("institution") else "")
    if not j and m.get("type") == "posted-content": j = (m.get("publisher") or "preprint")
    yr = (m.get("published-print") or m.get("published-online") or m.get("issued") or m.get("created"))["date-parts"][0][0]
    out[key] = dict(authors=au, groups=groups, title=" ".join((m.get("title") or [""])[0].split()), journal=j, year=str(yr), volume=m.get("volume", ""), pages=m.get("page", "") or m.get("article-number", ""),
                    doi=doi.strip(), type=m.get("type", ""))
    log.append(f"{key}\t{doi}\t{yr}\t{j}\t{out[key]['title'][:90]}")
    time.sleep(0.3)
out.update(MANUAL)
(P/"refs.json").write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
(P/"refs_check.txt").write_text("\n".join(log), encoding="utf-8")
print(f"{len(out)} refs, 실패 {sum('실패' in l for l in log)}")

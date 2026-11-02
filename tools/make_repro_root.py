"""Create a fresh data root for the reproduction run, reusing only the downloaded inputs.

The reproduction must not overwrite the original computed files, so it writes into a new root (REPRO) and reads the
downloaded inputs through directory junctions into the original root (ORIG). Junctions cost no disk space.
Nothing in ORIG is modified. Run once:  python tools/make_repro_root.py
Then run the pipeline with  AGA_DATA=D:/AD_GWAS_Atlas_repro  AGA_DATA_POSIX=/d/AD_GWAS_Atlas_repro.

Inputs linked (everything else is recomputed by the pipeline):
  Atlas query results            exp3/atlas, exp3c/atlas, exp3l/atlas_active   (AlphaGenome Atlas API; not redistributed)
  Reference panels and LDSC data ref_1000G, exp3/GRCh38 (baselineLD v2.2, plink files, weights), exp3/1000G_Phase3_*
  Gene annotation                ref (GENCODE v46)
  GWAS summary statistics        sumstats (Bellenguez 2022), exp3b/sumstats (Price-lab PASS files, monocyte count),
                                 exp3f/sumstats/raw (validation traits), exp3n/raw, exp3r/raw (larger GWAS)
  Measured annotations           exp3e/raw (FANTOM5, ENCODE)
  Cell-type eQTLs                exp3s/raw (FunGen-xQTL)
  Software                       tools (ldsc, ldsc_bv, plink2)
"""
import subprocess, sys
from pathlib import Path

ORIG = Path(sys.argv[1] if len(sys.argv) > 1 else "D:/AD_GWAS_Atlas_data")
REPRO = Path(sys.argv[2] if len(sys.argv) > 2 else "D:/AD_GWAS_Atlas_repro")
LINKS = ["exp3/atlas", "exp3/GRCh38", "exp3/1000G_Phase3_frq", "exp3/1000G_Phase3_weights_hm3_no_MHC",
         "exp3c/atlas", "exp3l/atlas_active", "ref_1000G", "ref", "sumstats", "exp3b/sumstats", "exp3f/sumstats/raw",
         "exp3n/raw", "exp3r/raw", "exp3e/raw", "exp3s/raw", "tools"]
COPY_FILES = ["exp3/readme_baseline_versions.txt", "exp3/hm3_no_MHC.list.txt",   # the HapMap3 no-MHC SNP list is a download (Zenodo 10515792)
              "exp3c/track_lineage.csv",                                          # track metadata and lineage assignment, written by z01_lineages.py (API script, not in the repository)
              "exp3c/test/part_0000.parquet"]                                     # small Atlas query used only by the chr22 test modes (V01_TEST etc.; form checks, no manuscript number)

def junction(dst, src):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists(): print("exists ", dst); return
    r = subprocess.run(["cmd", "/c", "mklink", "/J", str(dst), str(src)], capture_output=True, text=True)
    print(("linked " if r.returncode == 0 else "FAILED ") + str(dst), r.stderr.strip())

for rel in LINKS:
    assert (ORIG/rel).is_dir(), ORIG/rel
    junction(REPRO/rel, ORIG/rel)
for rel in COPY_FILES:
    (REPRO/rel).parent.mkdir(parents=True, exist_ok=True)
    if not (REPRO/rel).exists(): (REPRO/rel).write_bytes((ORIG/rel).read_bytes()); print("copied ", REPRO/rel)
for d in ["exp3", "exp3b", "exp3c", "exp3d", "exp3e", "exp3f", "exp3g", "exp3h", "exp3i", "exp3j", "exp3k", "exp3l", "exp3m", "exp3n", "exp3o", "exp3p", "exp3q", "exp3r", "exp3s"]:
    (REPRO/d).mkdir(exist_ok=True)
print("done:", REPRO)

# Data layout

All scripts read and write under one data root, given by the environment variable `AGA_DATA` (default
`D:/AD_GWAS_Atlas_data`; the shell runner also takes `AGA_DATA_POSIX`, the same folder in POSIX form, for example
`/d/AD_GWAS_Atlas_data` in Git Bash). Nothing in this repository is read from or written to any other place.

The subfolders are named after the analysis stages of the working folders (`exp3` … `exp3w`). The table lists what
must be present before the pipeline starts (inputs) and what the pipeline writes.

## Inputs (downloaded or queried once)

| Path under `AGA_DATA` | Contents | Source |
|---|---|---|
| `ref_1000G/all_hg38.pgen`, `all_hg38_noannot.pvar.zst`, `hg38_corrected.psam` | 1000 Genomes high-coverage (hg38) plink2 files | plink2 resources, 1000 Genomes 30x on GRCh38 |
| `exp3/GRCh38/baselineLD_v2.2/`, `plink_files/`, `weights/` | baselineLD v2.2 annotations and LD scores (hg38), 1000G EUR plink files, HapMap3 regression weights | [Zenodo 10515792](https://doi.org/10.5281/zenodo.10515792) |
| `exp3/1000G_Phase3_frq/`, `exp3/1000G_Phase3_weights_hm3_no_MHC/` | allele frequencies and weights (hg19 release, used for the M_5_50 counts and the SNP list) | Price lab LDSC resources |
| `exp3/hm3_no_MHC.list.txt` | HapMap3 SNPs without the MHC (1,217,311 rsIDs) | Zenodo 10515792 |
| `ref/gencode.v46.annotation.gtf.gz` | GENCODE v46 (TSS positions, e10) | GENCODE |
| `sumstats/GCST90027158_buildGRCh38.tsv.gz` | Bellenguez et al. 2022 AD stage I | GWAS Catalog |
| `exp3b/sumstats/` | Price lab processed GWAS `PASS_Height1.sumstats`, `PASS_BMI1.sumstats`, `PASS_Bipolar_Disorder.sumstats`, `PASS_IBD.sumstats`, `PASS_LDL.sumstats`, `PASS_Neuroticism.sumstats`, `PASS_Rheumatoid_Arthritis.sumstats`, `PASS_Schizophrenia.sumstats`, `PASS_Years_of_Education2.sumstats`, and `GCST90002393.h.tsv.gz` (monocyte count) | Price lab; GWAS Catalog |
| `exp3f/sumstats/raw/` | validation traits: BMD (GCST006979), neutrophil count (GCST90002398), total cholesterol (GCST90500025) | GWAS Catalog |
| `exp3n/raw/` | larger GWAS: height (Yengo 2022), BMI (Yengo 2018), schizophrenia (PGC3), bipolar disorder (PGC 2024), LDL (GLGC 2021) | consortia (see `experiments/e13_larger_gwas/scripts/n00_download.sh`) |
| `exp3r/raw/GCST90704647.tsv.gz` | 2026 AD consensus meta-analysis without proxy cases | GWAS Catalog (`e17/scripts/r00_download.sh`) |
| `exp3e/raw/fantom5/`, `exp3e/raw/encode/` | FANTOM5 CAGE peaks and enhancers, ENCODE DNase and H3K27ac peaks | downloaded by `e04/scripts/w01_measured_sources.py` (skips files already present) |
| `exp3s/raw/` | FunGen-xQTL fine-mapped eQTLs (NIAGADS NG00184.v1) | NIAGADS (`e18/scripts/s00_download.sh`) |
| `exp3u/raw/sumstats_107/` | panel of independent GWAS in LDSC format: `<name>.sumstats.gz` (106 files; the multiple sclerosis file of the panel needs separate approval and is not used) and `traits_indep107.xlsx` (the table of the panel) | [Zenodo 10515792](https://doi.org/10.5281/zenodo.10515792), `sumstats_indep107.tgz` (md5 `8406e1370e5807a11f9f04174ce9f76e`), extracted under `exp3u/raw/` |
| `exp3/atlas/part_*.parquet` | AlphaGenome Atlas allelic-effect scores, myeloid cell-type tracks, 5.83 M SNPs in 583 parts | Atlas API (`x02_atlas_query.py`, not in this repository) |
| `exp3c/atlas/part_*.parquet` | Atlas allelic-effect scores for all 1,748 tracks (1,312 of them assigned to the 12 lineages) | Atlas API (`z02_atlas_query.py`, not in this repository) |
| `exp3c/track_lineage.csv` | Track metadata with the lineage assigned to each track (Supplementary Table S1) | written by `z01_lineages.py` (not in this repository); the copy in `experiments/e02_lineage_map/inputs/track_lineage.csv` belongs at this path (`run_all.sh` puts it there when it is absent) |
| `exp3c/test/part_0000.parquet` | one small Atlas query (6.6 MB) read by the chr22 test modes of the annotation scripts (`*_TEST=1`; form checks only) | Atlas API (`z02_atlas_query.py`, not in this repository) |
| `exp3l/atlas_active/part_*.parquet` | Atlas activity scores, 784 tracks | Atlas API (`l01_atlas_active.py`, not in this repository) |
| `tools/ldsc/` | abdenlab/ldsc-python3 at commit b6176a4 with `tools/ldsc-python3_local_patch.diff` applied (LD scores, `munge_sumstats.py`) | GitHub |
| `tools/ldsc_bv/` | brainvisa/ldscore fork at commit a09a869 (regression) | GitHub |
| `tools/plink2.exe` | plink2 v2.0.0-a.7.8 | plink |

The Atlas parquet files hold the queried scores per SNP and track (`variant`, track ids, score columns). Their
redistribution is subject to the AlphaGenome terms of use; with one's own Atlas access they can be regenerated with the
scorers, windows and client version given in the Methods of the manuscript.

Three small input tables travel with the repository instead of the data root: `experiments/e01_targets_sldsc/inputs/loci.csv`
and `rep25_loci.csv` (lead variants of the 75 AD loci and the 25 replicated loci, from the public Bellenguez et al. 2022
summary statistics; used by x06 to define the excluded regions) and `experiments/e12_diagnostics/inputs/loci75_table_S3.csv`
(the same loci with their nearest genes; used by m04a to name the blocks). None contains individual-level data.

## Written by the pipeline

| Path | Written by | Contents |
|---|---|---|
| `exp3/targets.parquet`, `hc.*`, `target_pos.txt` | e01 x01 | 5.83 M target SNPs (MAF ≥ 5%, EUR) |
| `exp3/bellenguez_noAPOE.sumstats.gz` | e01 x00 | AD summary statistics, APOE removed, HapMap3 SNPs |
| `exp3*/annot/`, `exp3*/annot_sd_common.csv` | `*_annot*` scripts | annotation files (`.annot.gz`) and LD scores (`.l2.ldscore.gz`) per chromosome; standard deviation of each annotation over common SNPs (for τ*) |
| `exp3*/sldsc/`, `exp3*/runs/` | `*_sldsc*` scripts | S-LDSC outputs: `.results`/`.log` (e01–e08) or `.npz` with coefficients and jackknife delete values (e09 onwards) |
| `exp3t/sldsc/`, `exp3t/runs/` | e19 t01, t08 | LDSC outputs with `--overlap-annot` for the six later GWAS; all-track three-layer models of 14 GWAS in the tissue groups not fitted in e09 |
| `exp3u/sumstats/`, `exp3u/runs/` | e20 u01 | cleaned summary statistics of the 106 GWAS (`T001` … `T107`, without the file that is not used); `.npz` of the tissue-group, three-layer and resampled models |
| `exp3v/annld.npy`, `annld_snps.txt`, `annld_M.csv`, `exp3v/runs/` | e21 v01 | LD scores of the 72 annotations in one array; one `.npz` per GWAS with the 72 single-annotation fits |
| `exp3w/runs/`, `exp3w/runs72/`, `exp3w/annot/`, `exp3w/annot_table.parquet`, `exp3w/annot_sd_common.csv` | e22 w01, w02, w03 | `.npz` of the models of the 14 traits and of the enriched pairs of the panel; single-annotation fits of nine GWAS; activity-bin annotations and their LD scores |
| `exp3b/`, `exp3f/sumstats/`, `exp3n/sumstats/`, `exp3r/sumstats/` | x06, v04, n01, r01 | munged summary statistics |
| `exp3*/test/`, `exp3*/PIPELINE_STATUS*.txt`, `exp3*/STATUS_*.txt` | check stages and fitting scripts | the fits of the check stages; progress notes of long steps |
| `_status/`, `_logs/` | `run_all.sh` | done-markers, lock files and the log of every step |

Sizes after a full run: annotations and LD scores about 12 GB, S-LDSC outputs about 2 GB. The Atlas inputs are 48 GB, the
reference data 14 GB and the summary statistics of the panel of 106 GWAS 1 GB (e20 writes a cleaned copy of the same size).

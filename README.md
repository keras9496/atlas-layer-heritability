# Regulatory-layer profiles of complex-trait heritability from AlphaGenome Atlas predictions

Analysis code and result files for the preprint

> Byung Hyug Choi and Jaeseung Jeong. **Sequence-based variant-effect predictions reveal trait-specific
> regulatory-layer profiles of complex-trait heritability.** medRxiv (2026).
> doi: [10.64898/2026.10.16.26365786](https://doi.org/10.64898/2026.10.16.26365786)

If you use this code, please cite the preprint.

We grouped AlphaGenome Atlas allelic-effect predictions for 5.83 million common SNPs into 12 lineages and six regulatory
layers (CAGE, DNase, H3K27ac, H3K4me1, H3K4me3, RNA-seq) and partitioned the heritability of 11 complex traits, and
three validation traits, with stratified LD score regression (S-LDSC) conditional on baselineLD v2.2. The repository
holds every script that produced a number, figure or table in the preprint, in the order in which it was run, and
the result files themselves.

## Contents

| Path | What it is |
|---|---|
| `experiments/e01_…` to `e18_…` | One folder per analysis stage: `scripts/`, `results/` (the result files used in the preprint) and, in three stages, `inputs/` (small tables that are inputs rather than downloads: locus lists, track–lineage table) |
| `run_all.sh` | Runs every step in order with dependencies, done-markers and logs; resumable; two copies can share the work |
| `config.py` | The single place where the data root is set (`AGA_DATA`) |
| `data_layout.md` | What must be downloaded, where it goes, and what the pipeline writes |
| `CODE_INDEX.md` | What each script does and which figure or table it feeds |
| `paper/` | `make_figures.py` (draws all figures from the result files and writes `numbers.json`, the numbers quoted in the text) and the figures. The text of the preprint is available only from medRxiv |
| `verify/compare_results.py` | Compares a reproduction run with the result files in this repository; `verify/report.md` is the report of the run described under *Reproduction* |
| `tools/` | `ldsc-python3_local_patch.diff` (local changes to LDSC), `make_run_copy.py` (makes the copy to run in) |
| `requirements.txt` | Package versions |

## Running

1. Install Python 3.13 with the packages in `requirements.txt`, LDSC (two forks, see `data_layout.md`) and plink2.
2. Lay out the data as described in `data_layout.md` under one folder and point `AGA_DATA` at it. The AlphaGenome Atlas
   scores are the only input that cannot be downloaded from a public archive; they are queried from the Atlas API
   (the query scripts are not part of this repository, see below).
3. Make a copy of the repository to run in, so that the result files of the preprint stay untouched:

       python tools/make_run_copy.py            # -> repro_run/ with empty results/ folders

4. Run the pipeline from that copy. In Git Bash on Windows:

       cd repro_run
       AGA_DATA=D:/data AGA_DATA_POSIX=/d/data bash run_all.sh        # everything, in order
       bash run_all.sh list                                           # the steps and their dependencies
       bash run_all.sh e05_                                           # one stage

   The pipeline has 111 steps (analysis steps, the self-checks and dry runs each stage was preceded by, and the figure
   script). Each step writes a done-marker under `$AGA_DATA/_status/` and a log under `$AGA_DATA/_logs/`, so an
   interrupted run continues where it stopped, and the S-LDSC scripts skip fits whose `.results` file exists. The
   single-layer-map fits of `e03` take about two days; a second copy of `run_all.sh` started in another terminal takes
   the steps the first has not claimed, so that the rest of the pipeline runs beside them. The whole pipeline took about
   118 hours of computation on one workstation (24 cores, 96 GB; most steps use a single core). On that Windows
   workstation LDSC very occasionally (3 of about 2,000 fits) failed while reading intact LD score files; `y04_sldsc.py`
   retries such a fit, and for any other script simply starting `run_all.sh` again repeats the step, which skips the
   fits already done.

5. Check the result: `python verify/compare_results.py repro_run` compares every file under `experiments/*/results/`
   of the run with the files of this repository and writes `verify/report.md`. The run also redraws the figures and
   rewrites `numbers.json` in its own `paper/` folder.

All random choices (track resampling, bins, block alignment) use fixed seeds, LD scores and S-LDSC are deterministic,
so a run with the same inputs and package versions reproduces the result files exactly.

## Reproduction

The pipeline in this repository was run once more, from a copy made as above and against the downloaded inputs only,
on 2026-10-30 to 2026-11-02 (one power cut in between; the run resumed from its markers). `verify/report.md` is the
comparison with the result files of the preprint:

| Files | Outcome |
|---|---|
| 133 result tables and text outputs | identical (67 byte-identical; 66 identical in content, numbers within rtol 1e-6) |
| 33 run logs | identical except for timestamps, elapsed times and the data-root path |
| 2 run logs | differ for a documented reason without any value differing (see the report) |
| 1 output of the Atlas query script | not a pipeline output |
| `paper/`: 10 figures and `numbers.json` | byte-identical |

## Data

No data are stored in this repository. Everything the pipeline reads is listed in `data_layout.md` with its source;
in short:

| Input | Where to get it |
|---|---|
| AlphaGenome Atlas allelic-effect and activity scores | AlphaGenome Atlas API (Google DeepMind); the query scripts are not included, see below |
| baselineLD v2.2 annotations, LD scores, 1000G EUR plink files and HapMap3 weights (hg38) | [Zenodo 10515792](https://doi.org/10.5281/zenodo.10515792) |
| 1000 Genomes high-coverage GRCh38 genotypes | plink2 resources |
| GWAS summary statistics (AD, 10 further traits, 3 validation traits, larger GWAS) | GWAS Catalog, Price lab processed GWAS, consortium downloads (`n00_download.sh`, `r00_download.sh`); Supplementary Tables S1 and S16 of the preprint |
| FANTOM5 CAGE, ENCODE DNase and H3K27ac peaks | downloaded by `e04/scripts/w01_measured_sources.py` |
| FunGen-xQTL fine-mapped eQTLs | NIAGADS NG00184 (`s00_download.sh`) |
| GENCODE v46 | GENCODE |

The small tables that are inputs rather than downloads (locus lists, the track–lineage table) are in the `inputs/`
folders of e01, e02 and e12. The result files under `experiments/*/results/` are the tables and text summaries the
preprint's figures and tables were made from (3 MB in all); the S-LDSC run files, LD scores and annotations they were
computed from are not included and are rebuilt by the pipeline.

## Not included

- **AlphaGenome Atlas API access code.** Four scripts that query the Atlas are left out: `x02_atlas_query.py`,
  `z01_lineages.py` (track metadata; its lineage rules are in `e02_lineage_map/scripts/lineage_rules.py` verbatim),
  `z02_atlas_query.py` and `l01_atlas_active.py`. The scorers, windows and client version are described in the Methods
  of the preprint.
- **Data.** Atlas predictions and annotations derived from them, LD scores, S-LDSC run files and GWAS summary
  statistics are not included (`.gitignore`); see *Data* above. Redistribution of Atlas-derived annotations is
  subject to the AlphaGenome terms of use.
- **LDSC itself.** LD scores and munging used [abdenlab/ldsc-python3](https://github.com/abdenlab/ldsc-python3) at
  commit b6176a4 with the local changes in `tools/ldsc-python3_local_patch.diff`; regression used the
  [brainvisa/ldscore](https://github.com/brainvisa/ldscore) fork at commit a09a869, unmodified. Both are GPL-3.0.

The docstrings of some scripts name the dated analysis protocol of their stage (`PROTOCOL_*`); the protocol documents
are not part of this repository, and the design of each stage is described in the Methods of the preprint.

Code comments and the text outputs of the scripts are in Korean; the preprint and this documentation are in English.

## License

MIT (see `LICENSE`). LDSC, which this code calls but does not include, is GPL-3.0.

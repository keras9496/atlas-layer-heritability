# Tissue-group and regulatory-layer partitioning of complex-trait heritability with AlphaGenome Atlas predictions

Analysis code.

AlphaGenome Atlas allelic-effect predictions for 5.83 million common SNPs are grouped into 12 tissue groups and six
regulatory layers (CAGE, DNase, H3K27ac, H3K4me1, H3K4me3, RNA-seq), turned into annotations (the top 5% of SNPs by
predicted effect) and used to partition the heritability of complex traits with stratified LD score regression (S-LDSC)
conditional on baselineLD v2.2. The analyses cover 14 traits, six of them with an earlier and a later GWAS, and a public
panel of 106 GWAS. The repository holds the analysis scripts, in the order in which they are run, and the scripts that
draw the figures. Result files and figures are written by the pipeline and are not tracked in this repository.

## Contents

| Path | What it is |
|---|---|
| `experiments/e01_…` to `e22_…` | One folder per analysis stage: `scripts/`, an empty `results/` folder into which the scripts write, and, in three stages, `inputs/` (small tables that are inputs rather than downloads: locus lists, track–lineage tables) |
| `run_all.sh` | Runs every step in order with dependencies, done-markers and logs; resumable; several copies can share the work |
| `config.py` | The single place where the data root is set (`AGA_DATA`) |
| `data_layout.md` | What must be downloaded, where it goes, and what the pipeline writes |
| `CODE_INDEX.md` | What each script does, which scripts stand behind each analysis, and which result files each figure and table is made from |
| `figures/` | The figure scripts (`fig1.py` … `fig5.py`, `figS1.py` … `figS5.py`, run together by `make_figures.py`), which read only the result files and write to `figures/out/` |
| `tools/` | `ldsc-python3_local_patch.diff` (local changes to LDSC), `make_run_copy.py` (makes a second copy of the code to run in) |
| `verify/compare_results.py` | Compares the result files of two runs |
| `requirements.txt` | Package versions |

Stages e01–e18 build the annotations and LD scores and analyse 11 traits and three further traits with the GWAS first
used. Stages e19–e22 add the later GWAS of six traits to a benchmark of 14 traits (the later GWAS where there are two)
and analyse the panel of 106 GWAS; they reuse the annotations, LD scores and run files of the earlier stages. The scripts
and file names call the tissue groups lineages.

## System requirements

- **Operating system.** The pipeline was run on Windows 11 only, from Git Bash. `run_all.sh` is a bash script; one
  script calls the plink2 binary as `tools/plink2.exe` under the data root.
- **Python.** Python 3.13.5 with the packages in `requirements.txt` (numpy 2.4.2, pandas 2.3.3, scipy 1.17.1,
  pyarrow 23.0.1, matplotlib 3.10.8, requests 2.33.0, openpyxl 3.1.5). Other versions were not tested.
- **External tools.** LDSC in two forks: [abdenlab/ldsc-python3](https://github.com/abdenlab/ldsc-python3) at commit
  b6176a4 with the local changes in `tools/ldsc-python3_local_patch.diff` (LD scores, munging), and the
  [brainvisa/ldscore](https://github.com/brainvisa/ldscore) fork at commit a09a869, unmodified (regression).
  plink2 v2.0.0-a.7.8. All three are kept under the data root (`tools/ldsc/`, `tools/ldsc_bv/`, `tools/plink2.exe`).
- **Hardware.** No non-standard hardware. The pipeline was run on one workstation with 24 cores and 96 GB of memory;
  most steps use a single core. Disk: about 63 GB of inputs and about 14 GB written by the pipeline.

## Installation

1. Install Python and the packages: `pip install -r requirements.txt`.
2. Choose a folder as the data root and point `AGA_DATA` at it. Put the two LDSC forks at the commits above (the patch
   applied to the first) and plink2 under its `tools/` folder.
3. Lay out the data under the data root as described in `data_layout.md`.

Installing the software takes a few minutes. Obtaining the inputs takes much longer: the downloads are about 15 GB,
and the AlphaGenome Atlas scores (48 GB) are queried from the Atlas API (see *Data*).

## Data

No data are stored in this repository. Everything the pipeline reads is listed in `data_layout.md` with its source;
in short:

| Input | Where to get it |
|---|---|
| AlphaGenome Atlas allelic-effect and activity scores | AlphaGenome Atlas API (Google DeepMind); the query scripts are not included, see *Not included* |
| baselineLD v2.2 annotations, LD scores, 1000G EUR plink files and HapMap3 weights (hg38) | [Zenodo 10515792](https://doi.org/10.5281/zenodo.10515792) |
| 1000 Genomes high-coverage GRCh38 genotypes | plink2 resources |
| GWAS summary statistics of the 14 traits (earlier and later GWAS) | GWAS Catalog, Price lab processed GWAS, consortium downloads (`n00_download.sh`, `r00_download.sh`); listed in `data_layout.md` |
| Panel of 106 GWAS in LDSC format | [Zenodo 10515792](https://doi.org/10.5281/zenodo.10515792), `sumstats_indep107.tgz` |
| FANTOM5 CAGE, ENCODE DNase and H3K27ac peaks | downloaded by `e04/scripts/w01_measured_sources.py` |
| FunGen-xQTL fine-mapped eQTLs | NIAGADS NG00184 (`s00_download.sh`) |
| GENCODE v46 | GENCODE |

The small tables that are inputs rather than downloads (locus lists, the track–lineage table and the table of tracks per
tissue group and layer) are in the `inputs/` folders of e01, e02 and e12.

## Running

In Git Bash on Windows, from the repository:

    AGA_DATA=D:/data AGA_DATA_POSIX=/d/data bash run_all.sh        # everything, in order
    bash run_all.sh list                                           # the steps and their dependencies
    bash run_all.sh e05_                                           # one stage

The pipeline has 220 steps (analysis steps, the checks and dry runs that precede each stage, the pieces into which the
fits of e20–e22 are split, and the figure step). Each step writes a done-marker under `$AGA_DATA/_status/` and a log
under `$AGA_DATA/_logs/`, so an interrupted run continues where it stopped, and the S-LDSC scripts skip fits whose
output exists. A runner that is killed leaves the lock file of its step behind; delete `$AGA_DATA/_status/*.lock` before
starting again. A second copy of `run_all.sh` started in another terminal takes the steps the first has not claimed: the
single-layer-map fits of `e03` take about two days and can run beside the rest, and the fits of e20–e22 are split into 5
to 12 pieces each so that as many copies can run them side by side.

**Run time.** On the workstation described above the steps of e01–e18 took about 118 hours of computation and those of
e19–e22 about 88 hours summed over the pieces (about ten hours with ten copies of the runner). On that workstation LDSC
very occasionally (3 of about 2,000 fits in e01–e18) failed while reading intact LD score files; `y04_sldsc.py` retries
such a fit, and for any other script simply starting `run_all.sh` again repeats the step, which skips the fits already
done.

**Checks.** There is no demo on example data: every step needs the inputs. The checks of the code are steps of the
pipeline itself. Stages have check and dry-run steps that run before their fits, for example fits to a trait with
permuted *Z* scores and steps that derive again the values of an earlier stage, and they stop the run when a check
fails.

**Outputs.** Result tables and text summaries are written to `experiments/*/results/`, and the figures (PNG, PDF) and
the values drawn in them (`numbers_<figure>.json`) to `figures/out/`; `CODE_INDEX.md` lists the result files behind every
figure and table. The figure scripts set Arial (regular, bold and italic) and were run with matplotlib 3.10.8. They
check the positions of the text they draw and stop when labels overlap or leave the canvas, so with another font or
another matplotlib version they may stop at such a check, or draw with a substitute font where a script does not check;
matplotlib then prints a warning that the font was not found.

## Determinism

All random choices (track resampling, bins, block alignment) use fixed seeds, and LD scores and S-LDSC use no random
numbers. Two runs with the same inputs and package versions therefore write the same result files, apart from
timestamps, elapsed times and paths in the logs and from the last digits of S-LDSC estimates, which can differ between
fits. On the workstation described above, relative differences between repeated runs stayed below 1e-9 in the tables
computed from the fits and below 1e-6 in the tables of the mixed model, which is fitted by numerical optimisation, and
with the same operating system, font and matplotlib version the figures were reproduced byte for byte.

To compare two runs, make a second copy of the code with `python tools/make_run_copy.py [target]`, run the pipeline in
both, each with its own data root (the done-markers and the fitted models are kept there), and run
`python verify/compare_results.py <target>`: it compares every file under `experiments/*/results/` and `figures/out/`
(numbers within a relative tolerance of 1e-6, logs after masking timestamps, elapsed times and the data-root path) and
writes a report.

## Not included

- **AlphaGenome Atlas API access code.** Four scripts that query the Atlas are left out: `x02_atlas_query.py`,
  `z01_lineages.py` (track metadata; its lineage rules are in `e02_lineage_map/scripts/lineage_rules.py` verbatim and
  its output is `e02_lineage_map/inputs/z01_lineages.txt`), `z02_atlas_query.py` and `l01_atlas_active.py`. The scorers,
  windows and client version are described in the Methods of the manuscript, and `data_layout.md` describes the files
  these scripts write, which are inputs of the pipeline.
- **Data.** Atlas predictions and annotations derived from them, LD scores, S-LDSC run files and GWAS summary
  statistics are not included (`.gitignore`); see *Data* above. Redistribution of Atlas-derived annotations is
  subject to the AlphaGenome terms of use.
- **Result files and figures.** They are written by the pipeline and are not tracked (`.gitignore`).
- **LDSC itself.** Both forks are GPL-3.0; this code calls them but does not include them.
- **The manuscript.** Its text and the scripts that lay out its tables are not part of this repository.

The docstrings of some scripts name the dated analysis protocol of their stage (`PROTOCOL_*`); the protocol documents
are not part of this repository, and the design of each stage is described in the Methods of the manuscript.

Code comments and the text outputs of the scripts are in Korean; this documentation is in English. The test and
summary scripts of e01–e18 label their comparisons with the numbering in use when each stage was run (H1, H2 and so on);
`CODE_INDEX.md` says which scripts stand behind each analysis and each figure and table.

## License

MIT (see `LICENSE`). LDSC, which this code calls but does not include, is GPL-3.0.

# Multi-Objective Gene Clustering with NSGA-II

Multi-objective clustering of cancer gene expression data using the NSGA-II algorithm combined with biological knowledge from Gene Ontology, STRING, KEGG, and DisGeNET databases.

Thesis: *"Clustering multiobjetivo de datos de expresión génica de cáncer incorporando búsqueda local y conocimiento biológico a priori"* — Nicolás Mariángel Toledo, USACH, 2020.

## Requirements

- **R 4.5+** (tested with R 4.5.3)
- **Rtools45** (Windows only, for C++ compilation via Rcpp — required, not optional)
- **8GB RAM** minimum (16GB recommended for full datasets)

## Setup

### 1. Configure Credentials

Copy `.env.example` to `.env` and fill in your credentials:

```bash
cp .env.example .env
```

Edit `.env`:
```
DAVID_EMAIL=your.email@university.edu
DISGENET_API_KEY=your-api-key-here
```

**DAVID** — Register at https://davidbioinformatics.nih.gov/webservice/register.htm (requires organizational email, no gmail/yahoo/hotmail). Used for biological enrichment validation.

**DisGeNET** — Register at https://www.disgenet.com/ for an API key. See the [DisGeNET section](#disgenet-limitations) below for details on access restrictions.

### 2. Install R Packages

**Option A: Using renv (recommended)** — Installs exact package versions from lockfile:
```r
install.packages("renv")
renv::restore()
```

**Option B: From scratch** — Installs latest compatible versions:
```r
source("package_installer.r")
```

Either way takes 15-30 minutes on first run. renv is preferred because it guarantees the same versions that were tested.

### 3. Generate Sample Data

Sample datasets are not included in the repo and must be generated from the full datasets:

```r
source("data_sampling.r")
generate.samples("data/evaluation/Leukemia_GSE28497.csv.gz", number_of_samples = 1, sampling_percent = 5, minimum_row_size = 50, minimum_col_size = 500)
```

### 4. Verify Installation

```
Rscript RUN_TEST.r
```

Expected:
```
1. Loading nsga2.r... OK
2. Loading matrices_evaluation.r... OK
3. Loading evaluator.r... OK
4. Testing dataset load... OK
5. Translating gene IDs... OK
PROJECT IS FUNCTIONAL!
```

### 5. Rtools Path (Windows)

If Rtools45 is not at the default `C:\rtools45`, create `<R_HOME>/etc/Renviron.site`:

```
RTOOLS45_HOME=C:/your/rtools45/path
PATH="${RTOOLS45_HOME}/x86_64-w64-mingw32.static.posix/bin;${RTOOLS45_HOME}/usr/bin;${PATH}"
```

C++ compilation is **required** — `performance.cpp` provides `Rmcalc()` and `clusteringCalc()` which the fitness function depends on with no R fallback.

## Quick Start

```r
source("main.r")
run.default()
```

Runs NSGA-II with PLS local search on a small Leukemia sample (5% of GSE28497) using Gene Ontology. 3 runs, 2800 evaluations each. Expect ~90 seconds.

## Biological Databases

| ID | Database | Data Source | First Run |
|----|----------|-------------|-----------|
| `go` | Gene Ontology | `org.Hs.eg.db` R package (local) | ~10 min per dataset |
| `kegg` | KEGG Pathway | KEGG REST API (live) | ~30 sec per dataset |
| `string` | STRING v12.0 | Downloads from stringdb-downloads.org (~80MB) + Ensembl BioMart API | ~5 min first time, then cached |
| `disgenet_dis` | DisGeNET | Local bulk file (see below) | ~2 min per dataset |

All distance matrices are cached in `cache/` after first computation. Subsequent runs load from cache instantly.

### DisGeNET Limitations

DisGeNET moved to a freemium model after the thesis was written (2020). The bulk download files that the code originally used are now behind a paywall.

**Fallback data included**: The repository includes the original DisGeNET v7 data files from 2020 in `disgenet/`:
- `all_gene_disease_associations.tsv.gz` (21MB) — gene-disease association scores
- `disgenet.all.v7.entrez.gmt` (4MB) — pathway GMT file

These files allow the DisGeNET distance matrix to be computed without a license. The data is from 2020 (v7) and will not include newer gene-disease associations.

**To get current data**: Register at https://www.disgenet.com/ for an API key. The free academic plan has severe rate limits (1 query/day). A paid license is required for bulk access. The `disgenet2r` R package (`remotes::install_gitlab("medbio/disgenet2r")`) can query the API but is impractical for building full distance matrices under the free plan.

## DAVID Biological Validation

DAVID (Database for Annotation, Visualization and Integrated Discovery) is used for post-hoc biological validation of clustering results. It was NOT used during parameter tuning (only hypervolume was used for irace).

The project uses a custom HTTP client (`david_client.r`) instead of the deprecated `RDAVIDWebService` Java package (removed from Bioconductor in 2021). The client calls DAVID's API at `davidbioinformatics.nih.gov` via `httr`.

**Registration**: https://davidbioinformatics.nih.gov/webservice/register.htm — requires organizational email (no gmail/yahoo/hotmail).

**Known behavior**: DAVID's `getTermClusterReport` can timeout on large gene lists. The code handles this with retries and splits clusters >3000 genes via k-means.

## Local Search Algorithms

| ID | Algorithm | Best for |
|----|-----------|----------|
| `lmols` | Large-MOLS | Best overall (thesis finding) |
| `pls` | Pareto Local Search | Good with STRING/DisGeNET |
| `nmols` | Narrow-MOLS | — |
| `mosa` | Simulated Annealing | — |
| `ensemble` | Clustering Ensemble | — |

## Datasets

All 15 datasets from [CuMiDa](https://sbcb.inf.ufrgs.br/cumida) (Curated Microarray Database), originally sourced from GEO.

**Evaluation (5)**: Breast GSE89116, Renal GSE53757, Bladder GSE31189, Brain GSE50161, Prostate GSE6919_U95Av2

**Training (10)**: Leukemia GSE28497/GSE9476, Liver GSE22405/GSE60502, Prostate GSE6919_U95B/U95C, Breast GSE10797, Colorectal GSE44861/GSE77953, Ovary GSE6008

Format: Compressed CSV (`.csv.gz`), rows = samples, columns = genes (probe IDs), values = expression levels. Genes are translated to ENTREZ IDs via GPL platform files in `data/chip_to_entrez_id/`.

## Project Structure

### Core Algorithm
- `nsga2.r` — NSGA-II multi-objective optimization
- `local_search.r` — Local search strategies
- `evaluator.r` — Silhouette, hypervolume, DAVID biological validation
- `david_client.r` — Custom DAVID HTTP client (replaces deprecated RDAVIDWebService)
- `matrices_evaluation.r` — Expression and biological distance matrix computation
- `performance.cpp` — C++ fitness function optimizations (Rcpp)

### Utilities
- `main.r` — Entry point, convenience functions, batch processing, dataset definitions
- `file_utils.r` — File I/O
- `credentials.r` — Loads DAVID email and DisGeNET API key from `.env`
- `data_sampling.r` — Dataset sampling for experiments
- `gpl_chip_to_entrez_id.r` — Gene ID translation (probe ID → ENTREZ)
- `globals.r` — Constants and database identifiers

### Configuration
- `package_installer.r` — Dependency installer
- `.env` / `.env.example` — Credentials (gitignored)
- `run_instance.r` — CLI for batch runs
- `parameter_tuner.r` — irace hyperparameter optimization
- `RUN_TEST.r` — Installation verification

### Data Directories
- `data/training/` — Training datasets (`.csv.gz`)
- `data/evaluation/` — Evaluation datasets
- `data/chip_to_entrez_id/` — GPL platform gene mappings
- `data/training/samples/` — Generated 5% samples (not in repo, generate with `data_sampling.r`)
- `cache/` — Auto-generated distance matrices (`.rda`)
- `stringdb/` — STRING database files (auto-downloaded by STRINGdb package)
- `disgenet/` — DisGeNET v7 bulk data files (included as fallback)

## Experimental Pipeline

The full thesis pipeline (as described in Chapter 4):

1. **Data preparation**: Load CSV → translate to ENTREZ IDs → build expression + biological distance matrices
2. **Parameter tuning** (optional): `source("parameter_tuner.r")` — uses irace with hypervolume metric on training datasets
3. **Run experiments**: `calculate.results(list('GSE89116', 'GSE53757', ...), debug=TRUE)` — 13 runs per config per dataset
4. **Evaluate**: `evaluate.results(...)` — silhouette + hypervolume + DAVID enrichment
5. **Display**: `show.results(...)` — prints all metrics

Best configuration from the thesis: **KEGG + L-MOLS**, position 2 (inside evolutionary cycle). See `best_params` in `main.r`.

## Modernization Notes (March 2026)

This code was originally written for R 4.0.4 (2020-2021). The following changes were made to run on R 4.5.3:

- **RDAVIDWebService** replaced with custom `david_client.r` (HTTP client using `httr`/`xml2`). The Bioconductor package was removed in 2021 due to unmaintained Java/rJava dependency.
- **DAVID URL** updated from `david.ncifcrf.gov` (decommissioned March 2024) to `davidbioinformatics.nih.gov`.
- **STRING** updated from v11.0 to v12.0.
- **data.table** compatibility fixes for R 4.5: `$rn` column handling, join return type changes, `[-1]` indexing behavior.
- **R 4.5 strict mode** fixes: `is.na()` on lists, `if()` with length > 1 vectors.
- **GOSemSim** deprecated parameter fix: `OrgDb` → `annoDb`.
- **Credentials** moved from `mail_list.txt` to `.env` file pattern.
- **Email rotation** removed (DAVID no longer enforces 200 job/day limit).

## License

Academic/Research use. See individual R package licenses.

---

**Original code**: 2020-2021 | **Modernized**: March 2026 (R 4.5.3, custom DAVID client)

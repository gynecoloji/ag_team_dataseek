# Dataseek

A multi-agent system for searching public omics repositories, presenting ranked results, and downloading datasets with metadata and pipeline-ready starter configs.

## What it does

Dataseek orchestrates specialized agents to query multiple omics databases in parallel, merges and ranks results by quality, and delivers datasets with all the context you need to run them through a downstream pipeline.

Three slash commands cover the full workflow:

| Command | What it does |
|---------|-------------|
| `/dataseek-search` | Search GEO, CCLE, Xena, DepMap, and Single Cell Portal |
| `/dataseek-download` | Download datasets with metadata, paper info, and pipeline configs |
| `/dataseek-sample-mine` | Extract per-sample metadata from the linked publication |

---

## Quick start

```bash
# Search for single-cell RNA-seq datasets on glioblastoma
/dataseek-search --omic scRNAseq --disease glioblastoma --organism human

# Download a result
/dataseek-download GSE123456

# Mine the paper for detailed sample metadata
/dataseek-sample-mine GSE123456
```

---

## Supported omic types

`scRNAseq` · `bulkRNAseq` · `scMultiomicseq` · `scGenomicseq` · `bulkGenomicseq` · `ATACseq` · `ChIPseq` · `CRISPR`

---

## Data sources

| Source | What's there | Auth |
|--------|-------------|------|
| GEO (NCBI) | General omics, all types | Optional: `NCBI_API_KEY` |
| CCLE | Cancer cell line expression | None |
| Xena UCSC | TCGA, GTEx, and curated datasets | None |
| DepMap | Cancer dependency and genomics data | None |
| Single Cell Portal | Single-cell datasets (Broad Institute) | Optional: `SCP_TOKEN` |

Source selection is automatic based on omic type. Override with `--source GEO,SCP`.

---

## Agent architecture

```
User
 │
 ├─ /dataseek-search ──► dataseek-coordinator
 │                          ├─ dataseek-geo    → search_geo.py
 │                          ├─ dataseek-ccle   → search_ccle.py
 │                          ├─ dataseek-xena   → search_xena.py
 │                          ├─ dataseek-depmap → search_depmap.py
 │                          └─ dataseek-scp    → search_scp.py
 │
 ├─ /dataseek-download ──► dataseek-download
 │                          └─ download_{source}.py
 │
 └─ /dataseek-sample-mine ──► dataseek-sample-mine
                               ├─ supplement_fetch.py
                               └─ MCP tools (PubMed, bioRxiv)
```

Source agents run in parallel. Results are merged, deduplicated, and ranked by a composite quality score (sample count, metadata completeness, linked paper, recency, disease/tissue match).

---

## Search parameters

```
/dataseek-search --omic       scRNAseq         (required)
                 --disease    glioblastoma
                 --organism   human            (default)
                 --tissue     brain
                 --source     GEO,SCP          (override auto-routing)
                 --min-samples 10
                 --date-from  2022-01-01
                 --date-to    2025-12-31
```

---

## Output structure

```
dataseek/
├── results/
│   └── search_cache/        # Cached search results (JSON)
└── downloads/
    └── {accession}/
        ├── data/            # Raw data files
        ├── metadata/        # study_metadata.json, paper.json, sample_metadata.csv
        └── pipeline_config/ # Ready-to-use pipeline starter files
```

Downloaded datasets include a `sample_metadata.csv` pre-filled from source metadata, with `TODO` fields for anything that needs manual review before running a pipeline.

---

## Sample mining

`/dataseek-sample-mine` reads the linked publication to extract per-sample clinical and experimental metadata that isn't captured in the repository record itself.

It accepts an accession ID, DOI, or PubMed ID:

```bash
/dataseek-sample-mine GSE224681
/dataseek-sample-mine 10.1038/s41586-024-xxxx
/dataseek-sample-mine PMID:38654321
```

Outputs to `results/sample_mining/`:
- `{id}_samples.csv` — all extracted sample fields (flexible schema, all columns kept)
- `{id}_extraction_report.md` — what was found, where, and any gaps

Full-text access is attempted via PMC first, then bioRxiv/medRxiv, then abstract-only + supplementary tables.

---

## Setup

```bash
conda env create -f environment.yml
conda activate dataseek
```

Optional environment variables:
```bash
export NCBI_API_KEY=your_key   # Higher GEO rate limits
export SCP_TOKEN=your_token    # Single Cell Portal access
```

---

## Testing

```bash
conda run -n dataseek pytest tests/ -v
```

---

## Validation & benchmarking

### Test suite overview

The test suite covers all five search sources, all five downloaders, shared utilities, and end-to-end integration. HTTP calls are mocked with the `responses` library so no live API keys are needed for unit and integration tests.

| Module | Test file | What's validated |
|--------|-----------|-----------------|
| `utils.py` | `test_utils.py` | HTTP retry logic, exponential backoff, DOI resolution via Crossref, cache save/load, file download, logger setup |
| `search_geo.py` | `test_search_geo.py` | Entrez query construction, eSummary XML parsing, SOFT sample parsing, characteristic key mapping, 500-sample cap |
| `search_ccle.py` | `test_search_ccle.py` | CCLE DepMap API response parsing, omic filtering |
| `search_xena.py` | `test_search_xena.py` | Xena dataset listing, organism/omic filters |
| `search_depmap.py` | `test_search_depmap.py` | DepMap API response parsing |
| `search_scp.py` | `test_search_scp.py` | SCP REST API v1 response parsing |
| `sample_utils.py` | `test_sample_utils.py` | Sample table schema, CSV writing, summary generation |
| `supplement_fetch.py` | `test_supplement_fetch.py` | PMC supplementary table extraction and schema normalization |
| Download scripts | `test_download_*.py` | Per-source download logic, file writing, metadata output |
| GEO end-to-end | `test_integration.py` | Full CLI invocation of each search script via subprocess, JSON output validity |
| Sample table e2e | `test_sample_integration.py` | `search_geo` → sample table generation → CSV write → summary, all with mocked HTTP |

### Running specific subsets

```bash
# Unit tests only (fast, no network)
conda run -n dataseek pytest tests/ -v -m "not integration"

# Integration tests (invokes scripts as subprocesses, requires conda env)
conda run -n dataseek pytest tests/ -v -m integration

# Single module
conda run -n dataseek pytest tests/test_search_geo.py -v

# With coverage report
conda run -n dataseek pytest tests/ --cov=scripts --cov-report=term-missing
```

### Key validation checkpoints

**Search result schema** — every result returned by any source script is validated against the required field set:

```
accession, source, title, organism, omic_type, platform, disease,
tissue, sample_count, condition_groups, data_files, paper,
metadata_quality, date_submitted
```

**Sample table schema** — `sample_table` objects embedded in results are validated for:
```
total_samples, shown_samples, capped, columns, rows, summary
```

**GEO SOFT characteristic mapping** — verified that raw SOFT keys like `"tumor type"`, `"disease state"`, `"Sex"` correctly map to schema columns (`sample_type`, `disease`, `sex`). Two-pass matching (exact → substring) prevents short keys such as `"age"` from shadowing `"stage"`.

**Sample cap** — `parse_soft_samples` is tested to hard-cap at 500 samples regardless of SOFT file size.

**HTTP reliability** — `fetch_with_retry` is tested across success, two-failures-then-success, and exhausted-retries scenarios. Cache determinism is verified: same query params always produce the same cache filename; different params produce different filenames.

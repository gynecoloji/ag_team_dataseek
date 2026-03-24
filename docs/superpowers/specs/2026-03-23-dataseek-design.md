# Dataseek: Multi-Agent Public Omics Dataset Discovery & Download

## Overview

Dataseek is a multi-agent system that searches public omics repositories, presents summarized results, and downloads selected datasets with metadata and pipeline-ready starter configs. It uses a hybrid architecture: skills orchestrate the workflow, Python scripts handle API queries and file downloads.

## Data Sources

| Source | URL | Coverage |
|--------|-----|----------|
| GEO (NCBI) | ncbi.nlm.nih.gov/geo | Broadest: all omic types, all organisms |
| CCLE | depmap.org (Cancer Cell Line Encyclopedia) | Bulk expression/mutation for cancer cell lines |
| Xena UCSC | xena.ucsc.edu | TCGA, GTEx, PCAWG large cohort data |
| DepMap | depmap.org | CRISPR screens, drug sensitivity, multi-omic cancer cell lines |
| Single Cell Portal (SCP) | singlecell.broadinstitute.org | Single-cell omic datasets (scRNAseq, scATAC, scMultiome) |

## Agent Architecture

7 agents: 1 coordinator, 5 source searchers, 1 downloader.

### Search Phase Agents

**dataseek-coordinator**
- Role: Parse user query, apply source routing, dispatch source agents in parallel, merge/deduplicate results, rank by relevance, present tiered output
- Input: User command parameters
- Output: Tiered search results (top 10 detailed + remaining brief list), cached to `results/search_cache/`

**dataseek-geo**
- Role: Search GEO via NCBI Entrez E-utilities API
- Script: `scripts/search_geo.py`
- Output: `SearchResult[]` JSON

**dataseek-ccle**
- Role: Search CCLE via DepMap Portal REST API
- Script: `scripts/search_ccle.py`
- Output: `SearchResult[]` JSON

**dataseek-xena**
- Role: Search UCSC Xena via REST API + xenaPython
- Script: `scripts/search_xena.py`
- Output: `SearchResult[]` JSON

**dataseek-depmap**
- Role: Search DepMap portal for functional genomics data
- Script: `scripts/search_depmap.py`
- Output: `SearchResult[]` JSON

**dataseek-scp**
- Role: Search Single Cell Portal via REST API
- Script: `scripts/search_scp.py`
- Output: `SearchResult[]` JSON

### Download Phase Agent

**dataseek-download**
- Role: Download data files + metadata for selected datasets, generate pipeline-ready configs and summary reports
- Input: One or more accession IDs
- Scripts: `scripts/download_{source}.py`
- Output: Downloaded files in `downloads/{accession}/`, summary report in `results/reports/`

## Source Routing

The coordinator auto-selects relevant sources based on omic type. User can override with `--source`.

| Omic Type | GEO | CCLE | Xena | DepMap | SCP |
|-----------|-----|------|------|--------|-----|
| scRNAseq | yes | no | no | no | yes |
| scMultiomicseq | yes | no | no | no | yes |
| scGenomicseq | yes | no | no | no | yes |
| bulkRNAseq | yes | yes | yes | no | no |
| bulkGenomicseq | yes | no | yes | yes | no |
| ATACseq | yes | no | no | no | yes* |
| ChIPseq | yes | no | no | no | no |
| CRISPR | no | no | no | yes | no |

*SCP only for single-cell ATAC.

## Slash Commands

### `/dataseek-search`

```
/dataseek-search --omic scRNAseq --disease glioblastoma --organism human
```

**Parameters:**
- `--omic` (required): scRNAseq, bulkRNAseq, scMultiomicseq, scGenomicseq, bulkGenomicseq, ATACseq, ChIPseq, CRISPR
- `--disease`: Free text disease/condition filter
- `--organism`: human (default), mouse, or other
- `--tissue`: Tissue/cell type filter
- `--source`: Override auto source routing (comma-separated: GEO,SCP)
- `--min-samples`: Minimum sample count
- `--date-from` / `--date-to`: Publication date range (YYYY-MM-DD)

### `/dataseek-download`

```
/dataseek-download GSE123456 GSE789012
```

Takes one or more accession IDs from a previous search. Looks up source from search cache.

## Search Workflow

1. Coordinator parses parameters, applies source routing (or user `--source` override)
2. Dispatches relevant source agents in parallel as subagents
3. Each source agent calls its Python search script, which queries the source API and returns structured JSON
4. Coordinator merges all results, deduplicates by title/accession similarity
5. Ranks by relevance (sample count, metadata completeness, recency, disease match)
6. Presents tiered output:
   - **Top 10**: Full summary cards with all 5 info categories
   - **Remaining matches**: One-line list (accession, title, source, sample count)
7. Caches results to `results/search_cache/` for the download phase

## Download Workflow

1. Download agent looks up each accession's source from search cache
2. For each dataset, calls `scripts/download_{source}.py`
3. Script downloads:
   - Raw/processed data files (original format preserved)
   - Sample metadata
   - Paper metadata (via CrossRef DOI lookup)
4. Organizes into directory structure
5. Generates pipeline-ready starter config (sample_metadata.csv, config.yaml if applicable)
6. Produces per-dataset summary report

## Data Contracts

### SearchResult Schema

```json
{
  "accession": "GSE123456",
  "source": "GEO",
  "title": "Single-cell RNA-seq of glioblastoma tumors",
  "organism": "Homo sapiens",
  "omic_type": "scRNAseq",
  "platform": "10x Chromium 3' v3",
  "disease": "glioblastoma",
  "tissue": "brain",
  "sample_count": 12,
  "condition_groups": ["tumor", "normal"],
  "data_files": [
    {"type": "count_matrix", "format": "h5", "size_mb": 450}
  ],
  "paper": {
    "title": "...",
    "authors": "...",
    "journal": "...",
    "doi": "10.1234/example",
    "date": "2025-06-15",
    "abstract": "..."
  },
  "metadata_quality": "good",
  "date_submitted": "2025-05-01"
}
```

`metadata_quality` values: `good` (all key fields present), `partial` (some fields missing), `minimal` (only basic info available).

### Summary Report Contents

Each downloaded dataset gets a `{accession}_summary.md` containing:

1. **Dataset metadata** — accession ID, title, organism, tissue/cell type, disease/condition, sample count
2. **Experimental design** — omic type, platform/technology, single-end/paired-end, number of conditions/groups
3. **Data availability** — downloaded files, formats, sizes
4. **Original paper** — title, authors, journal, DOI, publication date, abstract summary
5. **Quality indicators** — sample sizes per group, metadata completeness, reprocessing notes

## Directory Structure

```
dataseek/
├── results/
│   ├── search_cache/                  # Cached search results (JSON)
│   │   └── 2026-03-23_scRNAseq_glioblastoma.json
│   └── reports/                       # Per-dataset summary reports
│       └── GSE123456_summary.md
├── downloads/
│   └── GSE123456/
│       ├── data/                      # Downloaded data files (original format)
│       │   ├── matrix.h5
│       │   └── barcodes.tsv.gz
│       ├── metadata/
│       │   ├── sample_metadata.csv    # Standardized sample metadata
│       │   ├── study_metadata.json    # Study-level info
│       │   └── paper.json             # Paper info
│       └── pipeline_config/
│           ├── sample_metadata.csv    # Formatted for target pipeline
│           └── config.yaml            # Starter config (if applicable)
├── scripts/
│   ├── search_geo.py
│   ├── search_ccle.py
│   ├── search_xena.py
│   ├── search_depmap.py
│   ├── search_scp.py
│   ├── download_geo.py
│   ├── download_ccle.py
│   ├── download_xena.py
│   ├── download_depmap.py
│   ├── download_scp.py
│   └── utils.py                       # Shared: HTTP, retry, logging, DOI resolver
├── CLAUDE.md                          # Agent orchestration spec
└── environment.yml                    # Conda environment definition
```

## Technical Stack

### Per-Source API Strategy

| Source | API | Auth | Python Libraries |
|--------|-----|------|-----------------|
| GEO | NCBI Entrez E-utilities (esearch, efetch) | API key optional (higher rate limits) | `biopython`, `requests` |
| CCLE | DepMap Portal REST API | No | `requests`, `pandas` |
| Xena UCSC | Xena REST API | No | `xenaPython`, `requests` |
| DepMap | DepMap Portal REST API + bulk downloads | No | `requests`, `pandas` |
| SCP | Single Cell Portal REST API v1 | Bearer token (free registration) | `requests` |

### Conda Environment (`dataseek`)

- `python>=3.10`
- `biopython`
- `requests`
- `pandas`
- `xenaPython`

### Script Interface Pattern

```bash
# Search: returns JSON to stdout
python scripts/search_geo.py --omic scRNAseq --disease glioblastoma --organism human --max-results 50

# Download: downloads files to specified directory
python scripts/download_geo.py --accession GSE123456 --output-dir downloads/GSE123456
```

Scripts are standalone, invoked by agents via subprocess.

### Shared Utilities (`scripts/utils.py`)

- HTTP client with retry logic (exponential backoff, 3 retries)
- Download with progress tracking and resume support (range headers)
- Logging with timestamps
- JSON/CSV standardization helpers
- DOI-to-paper metadata resolver (CrossRef API)

## Pipeline Config Generation

The download agent infers the target pipeline from omic type and generates starter configs:

| Omic Type | Target Pipeline Dir | Config Generated |
|-----------|-------------------|-----------------|
| scRNAseq | `../scRNAseq/` | `sample_metadata.csv` (sample_id, condition, batch) |
| bulkRNAseq | `../bulkRNAseq/` | `sample_metadata.csv` + `config.yaml` |
| ATACseq | `../ATACseq/` | `sample_metadata.csv` |
| ChIPseq | `../ChIPseq/` | `sample_metadata.csv` |
| scMultiomicseq | `../scMultiomicseq/` | `sample_metadata.csv` |
| scGenomicseq | `../scGenomicseq/` | `sample_metadata.csv` |
| bulkGenomicseq | `../bulkGenomicseq/` | `sample_metadata.csv` |

Configs are populated from source metadata where possible; unknown fields are marked for user completion.

## Error Handling

- **Source timeout during search**: Skip source, note in results that source was unreachable
- **Source returns zero results**: Report "no matches" for that source, continue with others
- **Download failure**: Retry 3x with exponential backoff, then report failure with error details
- **Malformed API response**: Log raw response, skip dataset, report parsing error
- **Search cache miss during download**: Re-search for the accession to determine source

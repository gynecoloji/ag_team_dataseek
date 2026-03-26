# Dataseek

Multi-agent system for searching public omics repositories, presenting summarized results, and downloading selected datasets with metadata and pipeline-ready starter configs.

## Architecture

Hybrid skill + Python script approach. Two slash commands (`/dataseek-search`, `/dataseek-download`) orchestrate 7 agents. Skills handle orchestration and user interaction; Python scripts handle API queries, response parsing, and file downloads.

## Data Sources

| Source | API | Auth |
|--------|-----|------|
| GEO (NCBI) | Entrez E-utilities (esearch, efetch, esummary) | Optional API key via `NCBI_API_KEY` env var |
| CCLE | DepMap Portal REST API (`depmap.org/portal/api/download/all`) | None |
| Xena UCSC | xenaPython + REST fallback | None |
| DepMap | DepMap Portal REST API | None |
| Single Cell Portal | REST API v1 (`singlecell.broadinstitute.org/single_cell/api/v1/`) | Bearer token via `SCP_TOKEN` env var |

## Conda Environment

Name: `dataseek`

```bash
conda run -n dataseek python scripts/{script}.py --args
```

Key packages: biopython, requests, pandas, xenaPython, pytest, responses

## Directory Structure

```
dataseek/
├── CLAUDE.md                         # This file (authoritative)
├── environment.yml                   # Conda env definition
├── .claude/skills/
│   ├── dataseek-search/SKILL.md      # Search orchestration skill
│   └── dataseek-download/SKILL.md    # Download orchestration skill
├── scripts/
│   ├── utils.py                      # Shared: HTTP retry, DOI resolver, caching, downloads
│   ├── search_geo.py                 # GEO Entrez search
│   ├── search_ccle.py                # CCLE DepMap Portal search
│   ├── search_xena.py                # Xena UCSC search
│   ├── search_depmap.py              # DepMap search
│   ├── search_scp.py                 # Single Cell Portal search
│   ├── sample_utils.py               # Sample metadata: schema, normalization, CSV, summaries
│   ├── supplement_fetch.py           # PMC supplementary table extraction
│   ├── download_geo.py               # GEO downloader
│   ├── download_ccle.py              # CCLE downloader
│   ├── download_xena.py              # Xena downloader
│   ├── download_depmap.py            # DepMap downloader
│   └── download_scp.py              # SCP downloader
├── tests/                            # pytest test suite
├── results/
│   ├── search_cache/                 # Cached search results (JSON) + per-dataset sample CSVs
│   └── reports/                      # Per-dataset summary reports (.md)
└── downloads/                        # Downloaded datasets
    └── {accession}/
        ├── data/                     # Data files (original format)
        ├── metadata/                 # study_metadata.json, paper.json, sample_metadata.csv
        └── pipeline_config/          # Pipeline-ready starter configs
```

## Agent Definitions

| Agent | Role | Trigger | Input | Output | Script |
|-------|------|---------|-------|--------|--------|
| dataseek-coordinator | Parse query, route sources, dispatch agents, merge/rank results | `/dataseek-search` invocation | User parameters | Tiered search results + cache | N/A (skill logic) |
| dataseek-geo | Search GEO via Entrez | Dispatched by coordinator | omic, organism, disease, tissue, dates | SearchResult[] JSON | `search_geo.py` |
| dataseek-ccle | Search CCLE via DepMap | Dispatched by coordinator | omic, disease | SearchResult[] JSON | `search_ccle.py` |
| dataseek-xena | Search Xena UCSC | Dispatched by coordinator | omic, organism, disease, tissue | SearchResult[] JSON | `search_xena.py` |
| dataseek-depmap | Search DepMap | Dispatched by coordinator | omic, disease | SearchResult[] JSON | `search_depmap.py` |
| dataseek-scp | Search Single Cell Portal | Dispatched by coordinator | omic, organism, disease, tissue | SearchResult[] JSON | `search_scp.py` |
| dataseek-download | Download datasets + generate configs | `/dataseek-download` invocation | Accession ID(s) | Downloaded files + reports | `download_{source}.py` |

## Source Routing

Auto-selected by coordinator based on omic type. User can override with `--source`.

| Omic Type | GEO | CCLE | Xena | DepMap | SCP |
|-----------|-----|------|------|--------|-----|
| scRNAseq | yes | no | no | no | yes |
| scMultiomicseq | yes | no | no | no | yes |
| scGenomicseq | yes | no | no | no | yes |
| bulkRNAseq | yes | yes | yes | no | no |
| bulkGenomicseq | yes | no | yes | yes | no |
| ATACseq | yes | no | no | no | yes* |
| ChIPseq | yes | no | no | no | no |
| CRISPR | yes | no | no | yes | no |

*SCP only for single-cell ATAC.

## SearchResult Schema

All search scripts output JSON arrays of this schema to stdout:

```json
{
  "accession": "GSE123456",
  "source": "GEO",
  "title": "...",
  "organism": "Homo sapiens",
  "omic_type": "scRNAseq",
  "platform": "10x Chromium 3' v3",
  "disease": "glioblastoma",
  "tissue": "brain",
  "sample_count": 12,
  "condition_groups": ["tumor", "normal"],
  "data_files": [{"type": "count_matrix", "format": "h5", "size_mb": 450}],
  "paper": {"title": "...", "authors": "...", "journal": "...", "doi": "...", "date": "...", "abstract": "..."},
  "metadata_quality": "good|partial|minimal",
  "date_submitted": "2025-05-01",
  "sample_table": {
    "total_samples": 12,
    "shown_samples": 12,
    "capped": false,
    "columns": ["sample_id", "tissue_site", "sample_type", "treatment_status", "disease", "cell_type", "age", "sex", "stage"],
    "rows": [["GSM123456", "brain", "primary", "untreated", "glioblastoma", "N/A", "45", "M", "IV"]],
    "summary": {"sample_type": {"primary": 8, "normal": 4}, "treatment_status": {"treated": 3, "untreated": 9}}
  }
}
```

## Script Invocation

All scripts are standalone CLI tools invoked via subprocess:

```bash
# Search: returns JSON to stdout
conda run -n dataseek python scripts/search_geo.py --omic scRNAseq --disease glioblastoma --organism human --max-results 50

# Download: downloads files to specified directory, returns summary JSON to stdout
conda run -n dataseek python scripts/download_geo.py --accession GSE123456 --output-dir downloads/GSE123456
```

## Slash Commands

### `/dataseek-search`

```
/dataseek-search --omic scRNAseq --disease glioblastoma --organism human
```

Parameters: `--omic` (required), `--disease`, `--organism` (default: human), `--tissue`, `--source`, `--min-samples`, `--date-from`, `--date-to`

### `/dataseek-download`

```
/dataseek-download GSE123456 SCP1234
```

Takes one or more accession IDs. Looks up source from cache or accession prefix.

## Pipeline Config Generation

Downloaded datasets get pipeline-ready starter configs in `downloads/{accession}/pipeline_config/`. Pipeline directories are siblings to dataseek (`../scRNAseq/`, `../bulkRNAseq/`, etc.).

## Error Handling

- **Source timeout during search:** Skip source, note in results
- **Source returns zero results:** Continue with other sources
- **Download failure:** Retry 3x with exponential backoff, then report
- **Malformed API response:** Log raw response, skip, report
- **Search cache miss during download:** Infer source from accession prefix

## Testing

```bash
conda run -n dataseek pytest tests/ -v
```

---
name: dataseek-search
description: Search public omics repositories for datasets matching user-specified omic type, disease, organism, and other filters
user_invocable: true
---

# /dataseek-search

Search public omics repositories (GEO, CCLE, Xena UCSC, DepMap, Single Cell Portal) for datasets matching user-specified filters.

## Usage

```
/dataseek-search --omic scRNAseq --disease glioblastoma --organism human
```

## Parameters

Parse the following from user args:

- `--omic` (REQUIRED): scRNAseq | bulkRNAseq | scMultiomicseq | scGenomicseq | bulkGenomicseq | ATACseq | ChIPseq | CRISPR
- `--disease`: Free text disease/condition filter
- `--organism`: human (default) | mouse | other
- `--tissue`: Tissue/cell type filter
- `--source`: Override auto source routing (comma-separated, e.g., GEO,SCP)
- `--min-samples`: Minimum sample count (post-filter)
- `--date-from` / `--date-to`: Publication date range (YYYY-MM-DD) — only passed to GEO

If `--omic` is missing, ask the user for it before proceeding.

## Source Routing

Unless `--source` is specified, auto-select sources based on omic type:

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

## Workflow

### Step 1: Dispatch Source Agents in Parallel

Dispatch all relevant source agents in a **single message** using the Agent tool. Each agent runs one search script:

**GEO Agent prompt:**
```
Search GEO for {omic} datasets. Run:
conda run -n dataseek python scripts/search_geo.py --omic {omic} --organism {organism} [--disease {disease}] [--tissue {tissue}] [--date-from {date_from}] [--date-to {date_to}] --max-results 50
Return the JSON output exactly as printed to stdout.
```

**CCLE Agent prompt:**
```
Search CCLE for {omic} datasets. Run:
conda run -n dataseek python scripts/search_ccle.py --omic {omic} [--disease {disease}] --max-results 50
Return the JSON output exactly as printed to stdout.
```

**Xena Agent prompt:**
```
Search Xena UCSC for {omic} datasets. Run:
conda run -n dataseek python scripts/search_xena.py --omic {omic} --organism {organism} [--disease {disease}] [--tissue {tissue}] --max-results 50
Return the JSON output exactly as printed to stdout.
```

**DepMap Agent prompt:**
```
Search DepMap for {omic} datasets. Run:
conda run -n dataseek python scripts/search_depmap.py --omic {omic} [--disease {disease}] --max-results 50
Return the JSON output exactly as printed to stdout.
```

**SCP Agent prompt:**
```
Search Single Cell Portal for {omic} datasets. Run:
conda run -n dataseek python scripts/search_scp.py --omic {omic} --organism {organism} [--disease {disease}] [--tissue {tissue}] --max-results 50
Return the JSON output exactly as printed to stdout.
```

### Step 2: Merge and Deduplicate

1. Parse JSON from each agent's response
2. Combine all results into one list
3. Deduplicate by accession (exact match) and title similarity (>90% match keeps the one with more metadata)
4. If `--min-samples` specified, filter out results where `sample_count < min_samples`

### Step 3: Rank Results

Sort by composite score (higher = better):
- `sample_count > 0`: +2 points
- `metadata_quality == "good"`: +2 points, `"partial"`: +1 point
- `paper.doi != ""`: +1 point
- `date_submitted` within last 2 years: +1 point
- Disease/tissue match in title: +1 point

### Step 4: Cache Results

Run this Python snippet to cache results:
```bash
conda run -n dataseek python -c "
import json, sys
sys.path.insert(0, '.')
from scripts.utils import save_search_cache
results = json.loads('''RESULTS_JSON''')
params = {PARAMS_DICT}
save_search_cache(results, '{omic}', params, 'results/search_cache')
"
```

### Step 5: Present Tiered Output

**Top 10 — Full Summary Cards:**

For each of the top 10 results, present:

```
### {rank}. {title}
**Accession:** {accession} | **Source:** {source} | **Samples:** {sample_count}
**Organism:** {organism} | **Platform:** {platform}
**Disease:** {disease} | **Tissue:** {tissue}
**Data Files:** {data_files summary}
**Paper:** {paper.title} ({paper.journal}, {paper.date}) DOI: {paper.doi}
**Quality:** {metadata_quality}
```

**Remaining Matches — Brief Table:**

```
| # | Accession | Title | Source | Samples |
|---|-----------|-------|--------|---------|
| 11 | GSE... | ... | GEO | 8 |
```

End with: "Use `/dataseek-download {accession}` to download any of these datasets."

## Error Handling

- If a source agent times out or fails: note it in output ("GEO was unreachable"), continue with other sources
- If all sources return 0 results: suggest broadening search terms
- If `--date-from`/`--date-to` specified: note that date filtering only applies to GEO results

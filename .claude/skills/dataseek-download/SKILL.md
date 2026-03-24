---
name: dataseek-download
description: Download selected omics datasets with metadata, paper info, and pipeline-ready starter configs
user_invocable: true
---

# /dataseek-download

Download one or more datasets from a previous `/dataseek-search`, including data files, metadata, paper info, and pipeline-ready starter configs.

## Usage

```
/dataseek-download GSE123456 SCP1234
```

Takes one or more accession IDs as arguments.

## Source Detection

For each accession, determine the source:

1. **Check search cache first:** Look in `results/search_cache/` for the most recent cache file containing this accession
2. **Fallback to prefix:** If not in cache, infer from accession format:
   - `GSE*` → GEO
   - `CCLE_*` → CCLE
   - `XENA_*` → Xena
   - `DepMap_*` → DepMap
   - `SCP*` → SCP

## Download Workflow

For each accession:

### Step 1: Run Download Script

```bash
conda run -n dataseek python scripts/download_{source}.py --accession {accession} --output-dir downloads/{accession} --omic-type {omic_type}
```

Where `{source}` is one of: geo, ccle, xena, depmap, scp.

### Step 2: Generate Summary Report

After download completes, read the metadata files and generate a summary report at `results/reports/{accession}_summary.md` containing:

1. **Dataset metadata** — accession, title, organism, tissue, disease, sample count
2. **Experimental design** — omic type, platform, conditions/groups
3. **Data availability** — downloaded files with formats and sizes
4. **Original paper** — title, authors, journal, DOI, date, abstract
5. **Quality indicators** — sample sizes per group, metadata completeness

Use `generate_summary_report()` from `scripts/download_geo.py` if the source is GEO. For other sources, generate equivalent markdown.

### Step 3: Generate Pipeline Config

Determine the target pipeline from omic type and generate starter configs in `downloads/{accession}/pipeline_config/`.

Pipeline directories are siblings to dataseek (`../scRNAseq/`, `../bulkRNAseq/`, etc.).

**Config generation rules:**

| Omic Type | Target Pipeline | Files Generated |
|-----------|----------------|-----------------|
| scRNAseq | `../scRNAseq/` | `sample_metadata.csv` (sample_id, condition, batch) |
| bulkRNAseq | `../bulkRNAseq/` | `sample_metadata.csv` + `config.yaml` |
| ATACseq | `../ATACseq/` | `sample_metadata.csv` |
| ChIPseq | `../ChIPseq/` | `sample_metadata.csv` |
| scMultiomicseq | `../scMultiomicseq/` | `sample_metadata.csv` |
| scGenomicseq | `../scGenomicseq/` | `sample_metadata.csv` |
| bulkGenomicseq | `../bulkGenomicseq/` | `sample_metadata.csv` |
| CRISPR | N/A | No pipeline config |

**sample_metadata.csv format:**
```csv
sample_id,condition,batch
GSM100001,TODO,TODO
GSM100002,TODO,TODO
```

Populate from source metadata where available. Mark unknown fields as `TODO` for user to fill in.

**bulkRNAseq config.yaml starter:**
```yaml
qc:
  min_count: 10
  min_samples: "auto"
de:
  padj_threshold: 0.05
  lfc_threshold: 1.0
pathway:
  consensus_min: 2
  ora_databases: ["GO:BP", "GO:MF", "GO:CC", "KEGG"]
  gsea_databases: ["MSigDB_Hallmark", "KEGG", "Reactome"]
organism: "human"
```

### Step 4: Report Status

Present download results to the user:

```
## Download Complete

### {accession} ({source})
- **Files downloaded:** {count} ({total_size_mb} MB)
- **Location:** downloads/{accession}/
- **Summary report:** results/reports/{accession}_summary.md
- **Pipeline config:** downloads/{accession}/pipeline_config/
- **Note:** Review sample_metadata.csv and fill in TODO fields before running the pipeline
```

## Error Handling

- **Download failure:** Retry 3x with exponential backoff, then report failure with details
- **Unknown accession prefix:** Ask user which source to use
- **Cache miss + unknown prefix:** Ask user for source
- **Multiple accessions:** Process sequentially to avoid overwhelming network

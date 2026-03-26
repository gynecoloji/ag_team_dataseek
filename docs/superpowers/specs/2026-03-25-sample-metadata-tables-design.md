# Sample-Level Metadata Tables for Dataseek Search

**Date:** 2026-03-25
**Status:** Approved

## Overview

Add per-sample clinical metadata extraction to all search scripts, producing structured tables with detailed sample information (tissue site, primary/metastasis, treatment status, etc.) during search — before any download. Supplement repository metadata with publication supplementary tables for richer clinical annotations.

## Decisions

- **When:** During search (not download)
- **Approach:** Per-script `fetch_samples()` in each search script, shared utilities in `utils.py`
- **Missing values:** Explicit `"N/A"`
- **Sample cap:** 500 per dataset
- **Cell line sources:** Unified standard columns + source-specific extras appended
- **Supplementary data:** Fetch from PMC when DOI/PMID available, takes priority over repository metadata

## Unified Sample Schema

### Standard Columns (all sources)

| Column | Description | Tier |
|--------|-------------|------|
| `sample_id` | Source-native sample ID (GSM*, cell line name, etc.) | 1 |
| `tissue_site` | Tissue/organ of origin | 1 |
| `sample_type` | primary / metastasis / cell_line / organoid / N/A | 1 |
| `treatment_status` | untreated / chemonaive / chemoresistant / chemosensitive / treated / N/A | 1 |
| `disease` | Disease or diagnosis | 2 |
| `cell_type` | Cell type annotation (esp. single-cell) | 2 |
| `age` | Patient age | 2 |
| `sex` | Patient sex | 2 |
| `stage` | Tumor stage/grade | 2 |

### Source-Specific Extra Columns

| Source | Extra Columns |
|--------|--------------|
| CCLE/DepMap | `lineage`, `sublineage`, `cosmic_id` |
| Xena | `cohort_name` |
| SCP | `cell_count`, `library_prep` |
| GEO | `platform_id` |

## Per-Source Fetch Strategy

### GEO

- After search, for each result call GEO SOFT endpoint: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={GSE}&form=text`
- Parse `!Sample_characteristics_ch1` lines to extract key-value pairs (submitter-formatted, e.g. `tissue: brain`, `treatment: cisplatin`)
- Map common keys to standard columns via case-insensitive keyword matching:
  - `tissue_site`: keys containing `tissue`, `organ`, `body site`, `anatomic site`
  - `sample_type`: keys containing `tumor type`, `sample type`, `primary`, `metastasis`, `metastatic`
  - `treatment_status`: keys containing `treatment`, `therapy`, `chemo`, `drug`
  - `disease`: keys containing `disease`, `diagnosis`, `pathology`, `histology`
  - `cell_type`: keys containing `cell type`, `cell line`, `sorted population`
  - `age`: keys containing `age`
  - `sex`: keys containing `sex`, `gender`
  - `stage`: keys containing `stage`, `grade`, `tnm`
- If >500 samples in SOFT record, parse first 500, note remainder in output
- Extra column: `platform_id`

### CCLE

- Fetch DepMap model metadata file (cell line annotations CSV) once per search session, cache it
- For each matched file, cross-reference cell lines in the data against the model metadata
- Map: `lineage` -> `tissue_site`, `"cell_line"` -> `sample_type`, `disease` -> `disease`
- Extra columns: `lineage`, `sublineage`, `cosmic_id`

### Xena

- After search, call `xenaPython.dataset_phenotypes(hub, dataset_name)` for each result
- Parse phenotype matrix for clinical columns (TCGA datasets have standardized `_primary_disease`, `sample_type`, `_PATIENT`, etc.)
- Map TCGA sample type codes (01=Primary, 06=Metastatic, etc.) to `sample_type`
- Extra column: `cohort_name`

### DepMap

- Same model metadata file as CCLE (shared cache)
- Same mapping strategy as CCLE
- Extra columns: `lineage`, `sublineage`, `cosmic_id`

### SCP

- Hit study detail endpoint `api/v1/studies/{accession}` for study metadata
- Extract available `cell_count`, `organ`, `disease`, `library_preparation_protocol`
- SCP doesn't expose per-cell metadata via API at search time — populate study-level attributes across all samples, mark per-cell fields `N/A`
- Extra columns: `cell_count`, `library_prep`

## Publication Supplementary Data Extraction

### Strategy

After resolving the paper DOI (already done via `resolve_doi()` in utils.py), attempt to fetch supplementary file listings from the publisher. Look for clinical/sample annotation tables and parse them to enrich the sample table.

### Publisher APIs

| Publisher | Method |
|-----------|--------|
| PubMed Central (PMC) | OA API: `https://www.ncbi.nlm.nih.gov/pmc/utils/oa.cgi?id={PMCID}` returns FTP links to supplementary files |
| CrossRef | DOI metadata includes `link` field with supplementary URLs |
| Direct journal sites | Fallback: scrape supplementary section from paper landing page |

### Flow

1. From the paper's DOI/PMID (already in SearchResult), resolve PMCID via Entrez
2. Fetch PMC supplementary file list via OA API
3. Filter for tabular files (`.csv`, `.tsv`, `.xlsx`, `.xls`) — skip images, PDFs, large archives
4. Download candidate tables (cap at 5 files, 10MB each)
5. Scan headers for sample ID columns and clinical attribute keywords (`patient`, `sample`, `primary`, `metastasis`, `treatment`, `stage`, `grade`, `age`, `sex`)
6. If a match is found, merge supplementary clinical data into the sample table, overwriting `N/A` values

### Priority Order

1. Publication supplementary tables (most curated)
2. Source API metadata (GEO SOFT, Xena phenotypes, etc.)
3. Mark remainder as `N/A`

### Limits

- Only attempt for datasets that have a resolved DOI/PMID
- Timeout per publisher: 15 seconds
- Cap at 5 supplementary files, 10MB each
- If no PMC access or no tabular supplements found, skip silently
- Log which datasets had supplementary enrichment in the search output

## Output

### CSV Files

- Each dataset's sample table saved to `results/search_cache/{accession}_samples.csv`
- Standard columns first, then source-specific extras
- First row is header, subsequent rows are samples (up to 500)
- If capped, last line of CSV is a comment: `# {N} additional samples not shown`

### Inline Summary

After each dataset's metadata block in search output, display:

```
Samples (12): 8 tumor primary, 4 normal | treated: 3, untreated: 9
```

Followed by a truncated markdown table (first 10 rows) with a note like `... and 2 more rows (see {accession}_samples.csv)`.

### SearchResult Schema Addition

New optional `sample_table` field:

```json
{
  "sample_table": {
    "total_samples": 12,
    "shown_samples": 12,
    "capped": false,
    "columns": ["sample_id", "tissue_site", "sample_type", "treatment_status", "disease", "cell_type", "age", "sex", "stage"],
    "rows": [
      ["GSM123", "brain", "primary", "untreated", "glioblastoma", "N/A", "45", "M", "IV"]
    ],
    "summary": {"primary": 8, "normal": 4, "treated": 3, "untreated": 9}
  }
}
```

If sample fetch fails or times out, omit `sample_table` and note in `metadata_quality`.

## Skill Orchestration

- After each search agent returns results, the coordinator writes CSV files and renders inline summaries
- No change to agent dispatch — sample fetching happens inside each search script

## utils.py Additions

- `normalize_sample_table(rows, source)` — maps source-specific columns to standard schema, fills N/A
- `write_sample_csv(accession, table_data, output_dir)` — writes CSV to cache directory
- `summarize_samples(table_data)` — generates the condensed summary string
- `fetch_supplementary_tables(doi, pmid)` — resolves PMCID, fetches and parses tabular supplements
- `match_supplement_to_samples(supplement_df, sample_ids)` — cross-references supplement with sample IDs
- `parse_tabular_file(filepath)` — reads CSV/TSV/XLSX into a standard format

## Testing Strategy

### Unit Tests

- Mock each source's sample metadata API response (SOFT record, Xena phenotypes, DepMap model CSV, SCP study detail)
- Verify column mapping: source-specific fields -> standard schema
- Verify N/A filling for missing attributes
- Verify 500-sample cap and "not shown" note

### Supplementary Extraction Tests

- Mock PMC OA API response with supplementary file list
- Mock tabular supplement download (CSV/TSV/XLSX)
- Verify header keyword matching and sample ID cross-referencing
- Verify merge priority (supplementary > source API > N/A)

### Integration Tests

- End-to-end search with sample table generation (mocked APIs)
- Verify CSV file written to `results/search_cache/{accession}_samples.csv`
- Verify `sample_table` field in SearchResult JSON
- Verify inline summary string format

### Utils Tests

- `normalize_sample_table()` with various source inputs
- `write_sample_csv()` file output correctness
- `summarize_samples()` output format
- Supplementary file parsing for CSV, TSV, XLSX formats

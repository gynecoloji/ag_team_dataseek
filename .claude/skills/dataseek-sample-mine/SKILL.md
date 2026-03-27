---
name: dataseek-sample-mine
description: Extract detailed sample-level information from a dataset's linked publication by scanning full text, tables, figures, and supplementary materials. Outputs a flexible-schema CSV.
user_invocable: true
---

# /dataseek-sample-mine

Extract detailed sample information from a dataset's publication. Scans the full paper — main text, tables, figure captions, and supplementary tables — to build a comprehensive sample metadata CSV.

## Usage

```
/dataseek-sample-mine GSE123456
/dataseek-sample-mine 10.1038/s41586-024-xxxx
/dataseek-sample-mine PMID:38654321
```

Takes one identifier: an accession ID (GSE/SCP), DOI, or PubMed ID.

## Workflow

### Step 1: Parse input and resolve identifiers

Determine the input type:
- **Accession ID** (starts with `GSE`, `SCP`, `CCLE_`, `XENA_`, `DepMap_`): resolve to DOI
- **DOI** (contains `/`): resolve to PMID
- **PMID** (starts with `PMID:` or is purely numeric): use directly

**Accession → DOI resolution:**
1. Search `results/search_cache/` for JSON files containing the accession. Use the Glob tool to find cache files, then Read each to check for the accession. Extract `paper.doi` from the matching result.
2. If not in cache, run the appropriate search script to get metadata:
   ```bash
   conda run -n dataseek python scripts/search_geo.py --omic bulkRNAseq --max-results 5 --accession {accession}
   ```
   Extract DOI from the result.
3. If no DOI found, report: "No linked publication found for {accession}. Try providing a DOI or PMID directly."

**DOI → PMID resolution:**
Use the PubMed MCP tool `convert_article_ids` to convert DOI to PMID. If that fails, try `lookup_article_by_citation` with the DOI.

**PMID → DOI resolution:**
Use the PubMed MCP tool `convert_article_ids` to get the DOI from the PMID.

Store both `pmid` and `doi` for use in subsequent steps. Either one being missing is OK — proceed with what's available.

### Step 2: Fetch full-text content (three-tier fallback)

**Tier 1 — PMC (PubMed Central):**
If PMID is available, call the PubMed MCP tool `get_full_text_article` with the PMID.
- If it returns full text: proceed with this content. Note `access_tier = "PMC"`.
- If it fails (paper not in PMC): fall through to Tier 2.

**Tier 2 — bioRxiv/medRxiv:**
If DOI is available, call the bioRxiv MCP tool `get_preprint` with the DOI.
- If it returns content: proceed. Note `access_tier = "bioRxiv"`.
- If it fails (not a preprint): fall through to Tier 3.

**Tier 3 — Abstract-only fallback:**
Call the PubMed MCP tool `get_article_metadata` with the PMID (or `search_articles` with the DOI if no PMID).
- Extract the abstract text.
- Also attempt supplementary table fetch (Step 5) — supplements may be open even if the paper is paywalled.
- Note `access_tier = "abstract-only"`.

### Step 3: Scan main text for sample descriptions

If full text was obtained (Tier 1 or 2), scan all sections. Focus on:

**Priority sections:** Methods, Materials and Methods, Patient Cohort, Study Design, Study Population, Sample Collection, Experimental Design, Results.

Extract from the text:
- Total cohort size and composition (e.g., "45 patients: 30 tumor, 15 normal")
- Sample collection details (e.g., "fresh-frozen tissue from surgical resection")
- Inclusion/exclusion criteria
- Any per-sample characteristics mentioned in running text
- Treatment protocols and response categories

Build a list of **text-extracted sample facts**. These supplement (not replace) tabular data found in Steps 4-5.

### Step 4: Scan tables — main article tables

If full text was obtained, identify all tables in the paper. For each table:

1. Read the table caption/title
2. Check if it contains sample-level data by looking for:
   - Column headers with keywords: sample, patient, subject, barcode, specimen, case, ID
   - Row counts suggesting per-sample data (>3 rows, <1000 rows)
   - Clinical/demographic column keywords: age, sex, gender, stage, grade, tissue, treatment, diagnosis, histology, response
3. If it looks like a sample table, extract ALL columns and rows as-is
4. Track the table number (e.g., "Table 1", "Table 2") for the `source_location` column

### Step 5: Scan tables — supplementary tables

Run the existing `supplement_fetch.py` to download supplementary tabular files:

```bash
conda run -n dataseek python -c "
import json
from scripts.supplement_fetch import fetch_supplementary_tables
rows = fetch_supplementary_tables(doi='{doi}', pmid='{pmid}')
print(json.dumps(rows) if rows else '[]')
"
```

If this returns data, it will be pre-matched sample rows with standardized schema fields. Include these rows and mark their `source_location` as "Supplementary tables (via supplement_fetch.py)".

Additionally, if the `supplement_fetch.py` approach misses tables (it only extracts tables matching its fixed schema), also check if any downloaded supplement files contain sample-level data with columns NOT in the standard schema. Read any supplementary CSV/TSV/Excel files in the temp directory and extract all columns from tables that have a sample ID column.

### Step 6: Scan figure captions

If full text was obtained, read all figure captions and legends. Look for per-sample annotations such as:
- "Patient 1 — Stage IV, treated with pembrolizumab"
- "Sample A: tumor, Sample B: adjacent normal"
- Color/symbol legends mapping to patient or sample identifiers

Extract any sample-level metadata from captions. Mark `source_location` as "Figure N caption".

### Step 7: Build flexible CSV

Merge all extracted sample data into a single table:

1. **Collect all sources:**
   - Main text tables (Step 4)
   - Supplementary tables (Step 5)
   - Figure caption extractions (Step 6)
   - Text-based descriptions (Step 3) — only if they provide per-sample facts not already in tables

2. **Merge strategy:**
   - If multiple sources share a `sample_id` column, perform an outer join on `sample_id`
   - If sources have different column names for the same concept (e.g., "Patient ID" vs "sample_id"), normalize to the most descriptive name
   - Keep ALL columns from ALL sources — do not enforce a fixed schema
   - Add a `source_location` column indicating the origin of each row's data

3. **First column:** `sample_id` (if identifiable). If no sample ID column exists across any source, use sequential row numbers and add a note.

4. **Sort** rows by `sample_id` (or row number).

5. **Write CSV** to `results/sample_mining/{identifier}_samples.csv` using the Write tool. Where `{identifier}` is:
   - The accession ID if input was an accession
   - A DOI-derived slug if input was DOI (replace `/` and `.` with `_`)
   - `PMID_{id}` if input was a PMID

### Step 8: Write extraction report

Write a markdown report to `results/sample_mining/{identifier}_extraction_report.md` with:

```markdown
# Sample Mining Report: {identifier}

## Paper
- **Title:** {title}
- **DOI:** {doi}
- **PMID:** {pmid}
- **Journal:** {journal}
- **Date:** {date}

## Access
- **Tier used:** {PMC | bioRxiv | abstract-only}
- **Full text available:** {yes/no}

## Extraction Summary
- **Sections scanned:** {list of sections found in paper}
- **Main tables found:** {count} ({count with sample data} contained sample data)
- **Supplementary files processed:** {count}
- **Figure captions scanned:** {count}

## Results
- **Total samples extracted:** {count}
- **Total columns/fields:** {count}
- **Column names:** {comma-separated list}

## Data Sources
| Source | Rows | Columns Added |
|--------|------|---------------|
| Table 1 | N | col1, col2, ... |
| Supp Table S2 | N | col3, col4, ... |
| ... | ... | ... |

## Gaps and Limitations
- {any gaps, e.g., "Supplementary Table S3 could not be parsed", "Paper is paywalled — only abstract available"}
```

### Step 9: Present results

Output a summary to the user:

```
## Sample Mining Complete: {identifier}

**Paper:** {title} ({journal}, {date})
**Access:** {tier}

**Extracted:** {N} samples × {M} fields
**CSV:** results/sample_mining/{identifier}_samples.csv
**Report:** results/sample_mining/{identifier}_extraction_report.md

**Fields found:** {comma-separated column names}

**Preview (first 5 rows):**
{formatted table of first 5 rows}
```

If no sample data was found:
```
## Sample Mining: {identifier}

**Paper:** {title} ({journal}, {date})
**Access:** {tier}

No structured sample data could be extracted.

**Text-based findings:** {summary of any cohort descriptions found in text}

**Suggestion:** The sample metadata for this dataset may need to be extracted manually from the paper's PDF tables.
```

## Error Handling

- **Identifier not resolvable:** "Could not resolve {input} to a publication. Please provide a DOI or PMID directly."
- **Paper not found in any tier:** "No publication content accessible for {identifier}. The paper may be paywalled with no open supplements."
- **supplement_fetch.py fails:** Continue with main text tables only. Note in report: "Supplementary table fetch failed: {error}"
- **No tables contain sample data:** Save text-based descriptions if any. Report "No structured sample tables found."
- **CSV write fails:** Report the error and print the extracted data to stdout as a fallback.

## Constraints

- This skill is **read-only** — it does not modify any existing dataseek data, search cache, or downloaded datasets.
- **Flexible schema** — never enforce a fixed set of columns. Extract whatever the paper provides.
- **Best-effort extraction** — papers vary widely in how they report sample info. Extract what's available, report gaps honestly.
- **No OCR** — does not attempt to read figure images, only captions and legends.
- **No web scraping** — only accesses papers through MCP tools and the existing supplement_fetch.py script.

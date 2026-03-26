# Sample-Level Metadata Tables Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add per-sample clinical metadata extraction to all search scripts, producing structured tables during search with tissue site, primary/metastasis, treatment status, and other clinical attributes.

**Architecture:** Each search script gets a `fetch_samples()` function that makes source-specific API calls after the initial search. Shared utilities in `scripts/sample_utils.py` handle schema normalization, CSV writing, summary generation, and supplementary table extraction. A new `scripts/supplement_fetch.py` handles PMC supplementary file discovery and parsing.

**Tech Stack:** Python 3.10+, requests, pandas (for tabular parsing), openpyxl (for XLSX), xenaPython, biopython, responses (testing)

---

### Task 1: Add openpyxl dependency

**Files:**
- Modify: `environment.yml:5-14`

- [ ] **Step 1: Add openpyxl to environment.yml**

Add `openpyxl` to conda dependencies for XLSX parsing of supplementary tables:

```yaml
dependencies:
  - python>=3.10
  - biopython
  - requests
  - pandas
  - openpyxl
  - pip
  - pip:
    - xenaPython
  - pytest
  - pytest-mock
  - responses
```

- [ ] **Step 2: Update conda environment**

Run: `conda env update -n dataseek -f environment.yml`
Expected: openpyxl installed successfully

- [ ] **Step 3: Commit**

```bash
git add environment.yml
git commit -m "chore: add openpyxl dependency for XLSX supplementary table parsing"
```

---

### Task 2: Create sample_utils.py with schema constants and normalize/write/summarize utilities

**Files:**
- Create: `scripts/sample_utils.py`
- Test: `tests/test_sample_utils.py`

- [ ] **Step 1: Write failing tests for normalize_sample_table**

```python
# tests/test_sample_utils.py
import csv
import os
import json

import pytest

from scripts.sample_utils import (
    STANDARD_COLUMNS,
    SOURCE_EXTRA_COLUMNS,
    SAMPLE_CAP,
    normalize_sample_table,
    write_sample_csv,
    summarize_samples,
    build_sample_table,
)


class TestNormalizeSampleTable:
    def test_fills_missing_columns_with_na(self):
        rows = [{"sample_id": "GSM001", "tissue_site": "brain"}]
        normalized = normalize_sample_table(rows, source="GEO")
        assert normalized[0]["sample_type"] == "N/A"
        assert normalized[0]["treatment_status"] == "N/A"
        assert normalized[0]["disease"] == "N/A"
        assert normalized[0]["cell_type"] == "N/A"
        assert normalized[0]["age"] == "N/A"
        assert normalized[0]["sex"] == "N/A"
        assert normalized[0]["stage"] == "N/A"

    def test_preserves_existing_values(self):
        rows = [{"sample_id": "GSM001", "tissue_site": "brain", "sample_type": "primary", "treatment_status": "untreated"}]
        normalized = normalize_sample_table(rows, source="GEO")
        assert normalized[0]["tissue_site"] == "brain"
        assert normalized[0]["sample_type"] == "primary"
        assert normalized[0]["treatment_status"] == "untreated"

    def test_appends_source_specific_columns_for_geo(self):
        rows = [{"sample_id": "GSM001", "platform_id": "GPL570"}]
        normalized = normalize_sample_table(rows, source="GEO")
        assert normalized[0]["platform_id"] == "GPL570"

    def test_appends_source_specific_columns_for_ccle(self):
        rows = [{"sample_id": "ACH-000001", "lineage": "lung", "sublineage": "NSCLC", "cosmic_id": "905949"}]
        normalized = normalize_sample_table(rows, source="CCLE")
        assert normalized[0]["lineage"] == "lung"
        assert normalized[0]["sublineage"] == "NSCLC"
        assert normalized[0]["cosmic_id"] == "905949"

    def test_caps_at_500_samples(self):
        rows = [{"sample_id": f"GSM{i:06d}"} for i in range(600)]
        normalized = normalize_sample_table(rows, source="GEO")
        assert len(normalized) == 500
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `conda run -n dataseek pytest tests/test_sample_utils.py -v`
Expected: FAIL with ModuleNotFoundError (scripts.sample_utils not found)

- [ ] **Step 3: Write sample_utils.py with constants and normalize_sample_table**

```python
# scripts/sample_utils.py
"""Shared utilities for sample-level metadata: schema, normalization, CSV output, summaries."""

import csv
import json
import os
from collections import Counter

STANDARD_COLUMNS = [
    "sample_id", "tissue_site", "sample_type", "treatment_status",
    "disease", "cell_type", "age", "sex", "stage",
]

SOURCE_EXTRA_COLUMNS = {
    "GEO": ["platform_id"],
    "CCLE": ["lineage", "sublineage", "cosmic_id"],
    "DepMap": ["lineage", "sublineage", "cosmic_id"],
    "Xena": ["cohort_name"],
    "SCP": ["cell_count", "library_prep"],
}

SAMPLE_CAP = 500


def normalize_sample_table(rows: list[dict], source: str) -> list[dict]:
    """Normalize sample rows to standard schema + source extras, fill N/A, cap at SAMPLE_CAP."""
    extras = SOURCE_EXTRA_COLUMNS.get(source, [])
    all_columns = STANDARD_COLUMNS + extras
    capped = rows[:SAMPLE_CAP]
    normalized = []
    for row in capped:
        norm_row = {}
        for col in all_columns:
            norm_row[col] = row.get(col, "N/A")
        normalized.append(norm_row)
    return normalized


def write_sample_csv(accession: str, table_data: dict, output_dir: str) -> str:
    """Write sample table to CSV. Returns file path."""
    os.makedirs(output_dir, exist_ok=True)
    filepath = os.path.join(output_dir, f"{accession}_samples.csv")
    columns = table_data["columns"]
    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        for row in table_data["rows"]:
            writer.writerow(row)
        if table_data.get("capped", False):
            remaining = table_data["total_samples"] - table_data["shown_samples"]
            f.write(f"# {remaining} additional samples not shown\n")
    return filepath


def summarize_samples(table_data: dict) -> str:
    """Generate a condensed summary string from sample table data."""
    total = table_data["total_samples"]
    rows = table_data["rows"]
    columns = table_data["columns"]

    parts = [f"Samples ({total})"]

    type_idx = columns.index("sample_type") if "sample_type" in columns else None
    treat_idx = columns.index("treatment_status") if "treatment_status" in columns else None

    if type_idx is not None:
        type_counts = Counter(r[type_idx] for r in rows if r[type_idx] != "N/A")
        if type_counts:
            type_parts = [f"{count} {label}" for label, count in type_counts.most_common()]
            parts.append(", ".join(type_parts))

    if treat_idx is not None:
        treat_counts = Counter(r[treat_idx] for r in rows if r[treat_idx] != "N/A")
        if treat_counts:
            treat_parts = [f"{label}: {count}" for label, count in treat_counts.most_common()]
            parts.append(" | ".join(treat_parts) if len(parts) <= 1 else "| " + ", ".join(treat_parts))

    return " | ".join(parts) if len(parts) > 1 else parts[0]


def build_sample_table(rows: list[dict], source: str, total_count: int | None = None) -> dict:
    """Build a sample_table dict from normalized rows."""
    normalized = normalize_sample_table(rows, source)
    extras = SOURCE_EXTRA_COLUMNS.get(source, [])
    columns = STANDARD_COLUMNS + extras
    total = total_count if total_count is not None else len(rows)
    shown = len(normalized)
    row_lists = [[r[col] for col in columns] for r in normalized]
    return {
        "total_samples": total,
        "shown_samples": shown,
        "capped": total > shown,
        "columns": columns,
        "rows": row_lists,
        "summary": _build_summary_dict(normalized),
    }


def _build_summary_dict(normalized: list[dict]) -> dict:
    """Build summary counts from sample_type and treatment_status."""
    summary = {}
    for field in ["sample_type", "treatment_status"]:
        counts = Counter(r[field] for r in normalized if r[field] != "N/A")
        summary.update(counts)
    return summary
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `conda run -n dataseek pytest tests/test_sample_utils.py::TestNormalizeSampleTable -v`
Expected: All 5 tests PASS

- [ ] **Step 5: Write failing tests for write_sample_csv**

Add to `tests/test_sample_utils.py`:

```python
class TestWriteSampleCsv:
    def test_writes_csv_with_header_and_rows(self, tmp_path):
        table_data = {
            "total_samples": 2,
            "shown_samples": 2,
            "capped": False,
            "columns": ["sample_id", "tissue_site", "sample_type"],
            "rows": [
                ["GSM001", "brain", "primary"],
                ["GSM002", "brain", "metastasis"],
            ],
        }
        filepath = write_sample_csv("GSE123", table_data, str(tmp_path))
        assert os.path.exists(filepath)
        with open(filepath) as f:
            reader = csv.reader(f)
            header = next(reader)
            assert header == ["sample_id", "tissue_site", "sample_type"]
            rows = list(reader)
            assert len(rows) == 2

    def test_writes_cap_comment_when_capped(self, tmp_path):
        table_data = {
            "total_samples": 600,
            "shown_samples": 500,
            "capped": True,
            "columns": ["sample_id"],
            "rows": [["GSM001"]],
        }
        filepath = write_sample_csv("GSE999", table_data, str(tmp_path))
        with open(filepath) as f:
            content = f.read()
            assert "# 100 additional samples not shown" in content
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `conda run -n dataseek pytest tests/test_sample_utils.py::TestWriteSampleCsv -v`
Expected: All 2 tests PASS

- [ ] **Step 7: Write failing tests for summarize_samples**

Add to `tests/test_sample_utils.py`:

```python
class TestSummarizeSamples:
    def test_generates_summary_with_type_and_treatment(self):
        table_data = {
            "total_samples": 4,
            "columns": ["sample_id", "tissue_site", "sample_type", "treatment_status"],
            "rows": [
                ["GSM001", "brain", "primary", "untreated"],
                ["GSM002", "brain", "primary", "treated"],
                ["GSM003", "brain", "metastasis", "untreated"],
                ["GSM004", "brain", "primary", "untreated"],
            ],
        }
        summary = summarize_samples(table_data)
        assert "Samples (4)" in summary
        assert "primary" in summary
        assert "untreated" in summary

    def test_handles_all_na_values(self):
        table_data = {
            "total_samples": 2,
            "columns": ["sample_id", "sample_type", "treatment_status"],
            "rows": [
                ["GSM001", "N/A", "N/A"],
                ["GSM002", "N/A", "N/A"],
            ],
        }
        summary = summarize_samples(table_data)
        assert "Samples (2)" in summary
```

- [ ] **Step 8: Run tests to verify they pass**

Run: `conda run -n dataseek pytest tests/test_sample_utils.py::TestSummarizeSamples -v`
Expected: All 2 tests PASS

- [ ] **Step 9: Write failing tests for build_sample_table**

Add to `tests/test_sample_utils.py`:

```python
class TestBuildSampleTable:
    def test_builds_complete_table_dict(self):
        rows = [
            {"sample_id": "GSM001", "tissue_site": "brain", "sample_type": "primary"},
            {"sample_id": "GSM002", "tissue_site": "lung", "sample_type": "metastasis"},
        ]
        table = build_sample_table(rows, source="GEO")
        assert table["total_samples"] == 2
        assert table["shown_samples"] == 2
        assert table["capped"] is False
        assert "sample_id" in table["columns"]
        assert "platform_id" in table["columns"]  # GEO extra
        assert len(table["rows"]) == 2
        assert table["summary"]["primary"] == 1
        assert table["summary"]["metastasis"] == 1

    def test_marks_capped_when_exceeds_limit(self):
        rows = [{"sample_id": f"GSM{i:06d}"} for i in range(600)]
        table = build_sample_table(rows, source="GEO", total_count=600)
        assert table["total_samples"] == 600
        assert table["shown_samples"] == 500
        assert table["capped"] is True
```

- [ ] **Step 10: Run tests to verify they pass**

Run: `conda run -n dataseek pytest tests/test_sample_utils.py::TestBuildSampleTable -v`
Expected: All 2 tests PASS

- [ ] **Step 11: Run full test suite to verify no regressions**

Run: `conda run -n dataseek pytest tests/test_sample_utils.py -v`
Expected: All 11 tests PASS

- [ ] **Step 12: Commit**

```bash
git add scripts/sample_utils.py tests/test_sample_utils.py
git commit -m "feat: add sample_utils with schema constants, normalization, CSV writer, and summary generator"
```

---

### Task 3: Create supplement_fetch.py for PMC supplementary table extraction

**Files:**
- Create: `scripts/supplement_fetch.py`
- Test: `tests/test_supplement_fetch.py`

- [ ] **Step 1: Write failing tests for resolve_pmcid**

```python
# tests/test_supplement_fetch.py
import os
import json
import tempfile

import pytest
import responses

from scripts.supplement_fetch import (
    resolve_pmcid,
    fetch_pmc_supplement_list,
    parse_tabular_file,
    match_supplement_to_samples,
    fetch_supplementary_tables,
    CLINICAL_KEYWORDS,
    SAMPLE_ID_KEYWORDS,
)


class TestResolvePmcid:
    @responses.activate
    def test_resolves_pmid_to_pmcid(self):
        responses.add(
            responses.GET,
            "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/",
            json={"records": [{"pmid": "38000001", "pmcid": "PMC10500001"}]},
            status=200,
        )
        result = resolve_pmcid(pmid="38000001")
        assert result == "PMC10500001"

    @responses.activate
    def test_returns_none_on_no_match(self):
        responses.add(
            responses.GET,
            "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/",
            json={"records": [{"pmid": "99999999"}]},
            status=200,
        )
        result = resolve_pmcid(pmid="99999999")
        assert result is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `conda run -n dataseek pytest tests/test_supplement_fetch.py::TestResolvePmcid -v`
Expected: FAIL with ModuleNotFoundError

- [ ] **Step 3: Write supplement_fetch.py with resolve_pmcid**

```python
# scripts/supplement_fetch.py
"""Fetch and parse publication supplementary tables for sample-level clinical metadata."""

import csv
import io
import os
import tempfile

import pandas as pd
import requests

from scripts.utils import fetch_with_retry, setup_logger

logger = setup_logger("supplement_fetch")

IDCONV_URL = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
PMC_OA_URL = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa.cgi"

TABULAR_EXTENSIONS = {".csv", ".tsv", ".xlsx", ".xls"}
MAX_SUPPLEMENT_FILES = 5
MAX_FILE_SIZE_MB = 10

SAMPLE_ID_KEYWORDS = {"sample", "patient", "subject", "barcode", "sample_id", "patient_id", "case_id", "specimen"}
CLINICAL_KEYWORDS = {
    "tissue_site": {"tissue", "organ", "body site", "anatomic site", "site"},
    "sample_type": {"tumor type", "sample type", "primary", "metastasis", "metastatic", "tumor"},
    "treatment_status": {"treatment", "therapy", "chemo", "drug", "regimen", "response"},
    "disease": {"disease", "diagnosis", "pathology", "histology", "cancer type"},
    "cell_type": {"cell type", "cell line", "sorted population", "cluster"},
    "age": {"age"},
    "sex": {"sex", "gender"},
    "stage": {"stage", "grade", "tnm"},
}


def resolve_pmcid(pmid: str = None, doi: str = None) -> str | None:
    """Resolve a PMID or DOI to a PMCID via NCBI ID Converter API."""
    params = {"format": "json", "tool": "dataseek", "email": "dataseek@example.com"}
    if pmid:
        params["ids"] = pmid
    elif doi:
        params["ids"] = doi
    else:
        return None
    try:
        resp = fetch_with_retry(IDCONV_URL, params=params, timeout=15, max_retries=2)
        if resp.status_code != 200:
            return None
        records = resp.json().get("records", [])
        if records and "pmcid" in records[0]:
            return records[0]["pmcid"]
    except Exception as e:
        logger.warning(f"Failed to resolve PMCID: {e}")
    return None


def fetch_pmc_supplement_list(pmcid: str) -> list[dict]:
    """Fetch list of supplementary files from PMC OA API. Returns list of {url, filename} dicts."""
    try:
        resp = fetch_with_retry(
            PMC_OA_URL,
            params={"id": pmcid, "format": "json"},
            timeout=15,
            max_retries=2,
        )
        if resp.status_code != 200:
            return []
        data = resp.json()
        records = data.get("records", [])
        if not records:
            return []
        record = records[0]
        supplements = []
        for key in ["supplementary_files", "supp_files"]:
            if key in record:
                supplements.extend(record[key])
        # Also check OA file list
        if "oa_pdf" in record or "oa_xml" in record:
            for file_info in record.get("files", []):
                href = file_info.get("href", "")
                ext = os.path.splitext(href)[1].lower()
                if ext in TABULAR_EXTENSIONS:
                    supplements.append({"url": href, "filename": os.path.basename(href)})
        return supplements[:MAX_SUPPLEMENT_FILES]
    except Exception as e:
        logger.warning(f"Failed to fetch PMC supplements for {pmcid}: {e}")
        return []


def parse_tabular_file(filepath: str) -> pd.DataFrame | None:
    """Read a CSV, TSV, or XLSX file into a DataFrame. Returns None on failure."""
    ext = os.path.splitext(filepath)[1].lower()
    try:
        if ext == ".csv":
            return pd.read_csv(filepath, nrows=500)
        elif ext == ".tsv":
            return pd.read_csv(filepath, sep="\t", nrows=500)
        elif ext in (".xlsx", ".xls"):
            return pd.read_excel(filepath, nrows=500, engine="openpyxl")
    except Exception as e:
        logger.warning(f"Failed to parse {filepath}: {e}")
    return None


def _column_matches_keywords(col_name: str, keywords: set[str]) -> bool:
    """Check if a column name matches any of the keywords (case-insensitive)."""
    col_lower = col_name.lower().strip()
    return any(kw in col_lower for kw in keywords)


def _find_sample_id_column(df: pd.DataFrame) -> str | None:
    """Find the column most likely to contain sample IDs."""
    for col in df.columns:
        if _column_matches_keywords(col, SAMPLE_ID_KEYWORDS):
            return col
    return None


def _map_columns_to_schema(df: pd.DataFrame) -> dict[str, str]:
    """Map DataFrame columns to standard schema columns. Returns {schema_col: df_col}."""
    mapping = {}
    for schema_col, keywords in CLINICAL_KEYWORDS.items():
        for df_col in df.columns:
            if _column_matches_keywords(df_col, keywords):
                mapping[schema_col] = df_col
                break
    return mapping


def match_supplement_to_samples(df: pd.DataFrame, sample_ids: list[str]) -> list[dict] | None:
    """Cross-reference a supplement DataFrame with sample IDs. Returns enriched rows or None."""
    id_col = _find_sample_id_column(df)
    if id_col is None:
        return None
    col_mapping = _map_columns_to_schema(df)
    if not col_mapping:
        return None
    # Normalize sample IDs for matching
    df_ids = df[id_col].astype(str).str.strip()
    sample_set = {s.strip() for s in sample_ids}
    matched_rows = []
    for idx, df_id in df_ids.items():
        if df_id in sample_set:
            row = {"sample_id": df_id}
            for schema_col, df_col in col_mapping.items():
                val = df.at[idx, df_col]
                row[schema_col] = str(val).strip() if pd.notna(val) else "N/A"
            matched_rows.append(row)
    return matched_rows if matched_rows else None


def fetch_supplementary_tables(doi: str = None, pmid: str = None, sample_ids: list[str] = None) -> list[dict] | None:
    """Full pipeline: resolve PMCID -> fetch supplements -> parse -> match to samples.
    Returns enriched sample rows or None."""
    if not doi and not pmid:
        return None
    pmcid = resolve_pmcid(pmid=pmid, doi=doi)
    if not pmcid:
        return None
    supplements = fetch_pmc_supplement_list(pmcid)
    if not supplements:
        return None
    for supp in supplements:
        url = supp.get("url", "")
        filename = supp.get("filename", "")
        ext = os.path.splitext(filename)[1].lower()
        if ext not in TABULAR_EXTENSIONS:
            continue
        try:
            resp = requests.get(url, timeout=15)
            if resp.status_code != 200:
                continue
            if len(resp.content) > MAX_FILE_SIZE_MB * 1024 * 1024:
                continue
            with tempfile.NamedTemporaryFile(suffix=ext, delete=False) as tmp:
                tmp.write(resp.content)
                tmp_path = tmp.name
            df = parse_tabular_file(tmp_path)
            os.unlink(tmp_path)
            if df is None:
                continue
            if sample_ids:
                matched = match_supplement_to_samples(df, sample_ids)
                if matched:
                    return matched
        except Exception as e:
            logger.warning(f"Failed to process supplement {filename}: {e}")
    return None
```

- [ ] **Step 4: Run resolve_pmcid tests to verify they pass**

Run: `conda run -n dataseek pytest tests/test_supplement_fetch.py::TestResolvePmcid -v`
Expected: All 2 tests PASS

- [ ] **Step 5: Write failing tests for parse_tabular_file**

Add to `tests/test_supplement_fetch.py`:

```python
class TestParseTabularFile:
    def test_parses_csv(self, tmp_path):
        csv_path = tmp_path / "test.csv"
        csv_path.write_text("sample_id,tissue,treatment\nS001,brain,untreated\nS002,lung,treated\n")
        df = parse_tabular_file(str(csv_path))
        assert df is not None
        assert len(df) == 2
        assert "sample_id" in df.columns

    def test_parses_tsv(self, tmp_path):
        tsv_path = tmp_path / "test.tsv"
        tsv_path.write_text("sample_id\ttissue\ttreatment\nS001\tbrain\tuntreated\n")
        df = parse_tabular_file(str(tsv_path))
        assert df is not None
        assert len(df) == 1

    def test_returns_none_for_invalid_file(self, tmp_path):
        bad_path = tmp_path / "test.csv"
        bad_path.write_text("not,valid\x00csv\x00content")
        df = parse_tabular_file(str(bad_path))
        # Should return a DataFrame or None, not crash
        assert df is None or isinstance(df, pd.DataFrame)
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `conda run -n dataseek pytest tests/test_supplement_fetch.py::TestParseTabularFile -v`
Expected: All 3 tests PASS

- [ ] **Step 7: Write failing tests for match_supplement_to_samples**

Add to `tests/test_supplement_fetch.py`:

```python
class TestMatchSupplementToSamples:
    def test_matches_sample_ids_and_maps_columns(self):
        df = pd.DataFrame({
            "Patient ID": ["S001", "S002", "S003"],
            "Tissue Site": ["brain", "lung", "liver"],
            "Treatment": ["untreated", "cisplatin", "untreated"],
            "Stage": ["IV", "III", "II"],
        })
        matched = match_supplement_to_samples(df, ["S001", "S003"])
        assert matched is not None
        assert len(matched) == 2
        assert matched[0]["sample_id"] == "S001"
        assert matched[0]["tissue_site"] == "brain"
        assert matched[0]["treatment_status"] == "untreated"

    def test_returns_none_when_no_sample_id_column(self):
        df = pd.DataFrame({
            "gene": ["TP53", "BRCA1"],
            "expression": [5.2, 3.1],
        })
        matched = match_supplement_to_samples(df, ["TP53"])
        assert matched is None

    def test_returns_none_when_no_clinical_columns(self):
        df = pd.DataFrame({
            "sample_id": ["S001", "S002"],
            "gene_count": [20000, 18000],
            "umi_count": [50000, 45000],
        })
        matched = match_supplement_to_samples(df, ["S001"])
        assert matched is None
```

- [ ] **Step 8: Run tests to verify they pass**

Run: `conda run -n dataseek pytest tests/test_supplement_fetch.py::TestMatchSupplementToSamples -v`
Expected: All 3 tests PASS

- [ ] **Step 9: Write failing test for full fetch_supplementary_tables pipeline**

Add to `tests/test_supplement_fetch.py`:

```python
class TestFetchSupplementaryTables:
    def test_returns_none_when_no_doi_or_pmid(self):
        result = fetch_supplementary_tables()
        assert result is None

    @responses.activate
    def test_returns_none_when_pmcid_not_found(self):
        responses.add(
            responses.GET,
            "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/",
            json={"records": [{"pmid": "99999999"}]},
            status=200,
        )
        result = fetch_supplementary_tables(pmid="99999999", sample_ids=["S001"])
        assert result is None
```

- [ ] **Step 10: Run tests to verify they pass**

Run: `conda run -n dataseek pytest tests/test_supplement_fetch.py::TestFetchSupplementaryTables -v`
Expected: All 2 tests PASS

- [ ] **Step 11: Run full supplement_fetch test suite**

Run: `conda run -n dataseek pytest tests/test_supplement_fetch.py -v`
Expected: All 10 tests PASS

- [ ] **Step 12: Commit**

```bash
git add scripts/supplement_fetch.py tests/test_supplement_fetch.py
git commit -m "feat: add supplement_fetch for PMC supplementary table extraction and sample matching"
```

---

### Task 4: Add fetch_samples to search_geo.py (GEO SOFT parsing)

**Files:**
- Modify: `scripts/search_geo.py:1-148`
- Test: `tests/test_search_geo.py`

- [ ] **Step 1: Write failing tests for GEO SOFT sample parsing**

Add to `tests/test_search_geo.py`:

```python
from scripts.search_geo import fetch_geo_samples, parse_soft_samples

SAMPLE_SOFT_RESPONSE = """\
^SERIES = GSE123456
!Series_title = Single-cell RNA-seq of glioblastoma
!Series_sample_id = GSM000001
!Series_sample_id = GSM000002
!Series_sample_id = GSM000003
^SAMPLE = GSM000001
!Sample_title = Tumor_1_primary
!Sample_characteristics_ch1 = tissue: brain
!Sample_characteristics_ch1 = disease state: glioblastoma
!Sample_characteristics_ch1 = tumor type: primary
!Sample_characteristics_ch1 = treatment: untreated
!Sample_characteristics_ch1 = age: 55
!Sample_characteristics_ch1 = Sex: Male
!Sample_characteristics_ch1 = stage: IV
!Sample_platform_id = GPL24676
^SAMPLE = GSM000002
!Sample_title = Tumor_2_metastasis
!Sample_characteristics_ch1 = tissue: spine
!Sample_characteristics_ch1 = disease state: glioblastoma
!Sample_characteristics_ch1 = tumor type: metastasis
!Sample_characteristics_ch1 = treatment: cisplatin
!Sample_characteristics_ch1 = age: 42
!Sample_characteristics_ch1 = Sex: Female
!Sample_characteristics_ch1 = stage: IV
!Sample_platform_id = GPL24676
^SAMPLE = GSM000003
!Sample_title = Normal_1
!Sample_characteristics_ch1 = tissue: brain
!Sample_characteristics_ch1 = disease state: normal
!Sample_characteristics_ch1 = cell type: astrocyte
!Sample_characteristics_ch1 = treatment: untreated
!Sample_characteristics_ch1 = age: 60
!Sample_characteristics_ch1 = Sex: Male
!Sample_platform_id = GPL24676
"""


class TestParseSoftSamples:
    def test_parses_sample_characteristics(self):
        rows = parse_soft_samples(SAMPLE_SOFT_RESPONSE)
        assert len(rows) == 3
        assert rows[0]["sample_id"] == "GSM000001"
        assert rows[0]["tissue_site"] == "brain"
        assert rows[0]["sample_type"] == "primary"
        assert rows[0]["treatment_status"] == "untreated"
        assert rows[0]["disease"] == "glioblastoma"
        assert rows[0]["age"] == "55"
        assert rows[0]["sex"] == "Male"
        assert rows[0]["stage"] == "IV"
        assert rows[0]["platform_id"] == "GPL24676"

    def test_maps_metastasis_sample(self):
        rows = parse_soft_samples(SAMPLE_SOFT_RESPONSE)
        assert rows[1]["sample_type"] == "metastasis"
        assert rows[1]["treatment_status"] == "cisplatin"
        assert rows[1]["tissue_site"] == "spine"

    def test_maps_cell_type(self):
        rows = parse_soft_samples(SAMPLE_SOFT_RESPONSE)
        assert rows[2]["cell_type"] == "astrocyte"

    def test_caps_at_500_samples(self):
        lines = ["^SERIES = GSE999\n"]
        for i in range(600):
            lines.append(f"^SAMPLE = GSM{i:06d}\n")
            lines.append(f"!Sample_characteristics_ch1 = tissue: brain\n")
        soft = "".join(lines)
        rows = parse_soft_samples(soft)
        assert len(rows) == 500
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `conda run -n dataseek pytest tests/test_search_geo.py::TestParseSoftSamples -v`
Expected: FAIL with ImportError (parse_soft_samples not found)

- [ ] **Step 3: Implement parse_soft_samples and fetch_geo_samples in search_geo.py**

Add after the existing imports at line 10 of `scripts/search_geo.py`:

```python
from scripts.sample_utils import build_sample_table, SAMPLE_CAP
from scripts.supplement_fetch import fetch_supplementary_tables
```

Add after `parse_geo_record` (after line 107) in `scripts/search_geo.py`:

```python
GEO_SOFT_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

CHARACTERISTIC_KEY_MAP = {
    "tissue_site": ["tissue", "organ", "body site", "anatomic site"],
    "sample_type": ["tumor type", "sample type", "primary", "metastasis", "metastatic"],
    "treatment_status": ["treatment", "therapy", "chemo", "drug"],
    "disease": ["disease", "diagnosis", "pathology", "histology", "disease state"],
    "cell_type": ["cell type", "cell line", "sorted population"],
    "age": ["age"],
    "sex": ["sex", "gender"],
    "stage": ["stage", "grade", "tnm"],
}


def _map_characteristic(key: str, value: str) -> tuple[str, str] | None:
    """Map a GEO SOFT characteristic key to a standard schema column."""
    key_lower = key.lower().strip()
    for schema_col, keywords in CHARACTERISTIC_KEY_MAP.items():
        if any(kw in key_lower for kw in keywords):
            return (schema_col, value.strip())
    return None


def parse_soft_samples(soft_text: str) -> list[dict]:
    """Parse GEO SOFT format text into sample metadata rows."""
    samples = []
    current_sample = None
    sample_count = 0

    for line in soft_text.split("\n"):
        line = line.strip()
        if line.startswith("^SAMPLE"):
            if current_sample is not None:
                samples.append(current_sample)
                sample_count += 1
                if sample_count >= SAMPLE_CAP:
                    break
            sample_id = line.split("=", 1)[1].strip() if "=" in line else ""
            current_sample = {"sample_id": sample_id}
        elif current_sample is not None:
            if line.startswith("!Sample_characteristics_ch1"):
                _, _, char_value = line.partition("=")
                char_value = char_value.strip()
                if ":" in char_value:
                    char_key, _, char_val = char_value.partition(":")
                    mapped = _map_characteristic(char_key, char_val)
                    if mapped:
                        current_sample[mapped[0]] = mapped[1]
            elif line.startswith("!Sample_platform_id"):
                _, _, plat = line.partition("=")
                current_sample["platform_id"] = plat.strip()

    if current_sample is not None and sample_count < SAMPLE_CAP:
        samples.append(current_sample)

    return samples


def fetch_geo_samples(accession: str, paper: dict = None) -> dict | None:
    """Fetch sample-level metadata for a GEO series. Returns sample_table dict or None."""
    try:
        resp = fetch_with_retry(
            GEO_SOFT_URL,
            params={"acc": accession, "form": "text", "view": "brief"},
            timeout=30,
            max_retries=2,
        )
        if resp.status_code != 200:
            logger.warning(f"GEO SOFT fetch failed for {accession}: {resp.status_code}")
            return None
        rows = parse_soft_samples(resp.text)
        if not rows:
            return None

        # Try supplementary enrichment
        if paper:
            pmid = paper.get("pubmed_id", "")
            doi = paper.get("doi", "")
            sample_ids = [r["sample_id"] for r in rows]
            supp_rows = fetch_supplementary_tables(doi=doi, pmid=pmid, sample_ids=sample_ids)
            if supp_rows:
                supp_map = {r["sample_id"]: r for r in supp_rows}
                for row in rows:
                    if row["sample_id"] in supp_map:
                        for key, val in supp_map[row["sample_id"]].items():
                            if key != "sample_id" and row.get(key, "N/A") == "N/A":
                                row[key] = val

        return build_sample_table(rows, source="GEO")
    except Exception as e:
        logger.warning(f"Failed to fetch samples for {accession}: {e}")
        return None
```

Modify `search_geo()` (around line 127-131) to call `fetch_geo_samples` for each result:

Replace:
```python
    results = []
    for doc in root.findall("DocSum"):
        accession = _get_item_text(doc, "Accession")
        if accession and accession.startswith("GSE"):
            results.append(parse_geo_record(doc, omic))
    logger.info(f"Parsed {len(results)} GEO series records")
    return results
```

With:
```python
    results = []
    for doc in root.findall("DocSum"):
        accession = _get_item_text(doc, "Accession")
        if accession and accession.startswith("GSE"):
            record = parse_geo_record(doc, omic)
            sample_table = fetch_geo_samples(accession, paper=record.get("paper"))
            if sample_table:
                record["sample_table"] = sample_table
            results.append(record)
    logger.info(f"Parsed {len(results)} GEO series records")
    return results
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `conda run -n dataseek pytest tests/test_search_geo.py::TestParseSoftSamples -v`
Expected: All 4 tests PASS

- [ ] **Step 5: Run full GEO test suite to verify no regressions**

Run: `conda run -n dataseek pytest tests/test_search_geo.py -v`
Expected: All existing tests still PASS

- [ ] **Step 6: Commit**

```bash
git add scripts/search_geo.py tests/test_search_geo.py
git commit -m "feat: add GEO SOFT sample metadata parsing with supplementary enrichment"
```

---

### Task 5: Add fetch_samples to search_ccle.py (DepMap model metadata)

**Files:**
- Modify: `scripts/search_ccle.py:1-78`
- Test: `tests/test_search_ccle.py`

- [ ] **Step 1: Write failing tests for CCLE sample fetching**

Add to `tests/test_search_ccle.py`:

```python
from scripts.search_ccle import fetch_ccle_samples

SAMPLE_MODEL_METADATA = [
    {
        "ModelID": "ACH-000001",
        "CellLineName": "A549",
        "OncotreeLineage": "Lung",
        "OncotreeSubtype": "Non-Small Cell Lung Cancer",
        "OncotreePrimaryDisease": "Lung Cancer",
        "COSMICID": "905949",
    },
    {
        "ModelID": "ACH-000002",
        "CellLineName": "MCF7",
        "OncotreeLineage": "Breast",
        "OncotreeSubtype": "Breast Invasive Ductal Carcinoma",
        "OncotreePrimaryDisease": "Breast Cancer",
        "COSMICID": "905946",
    },
]


class TestFetchCcleSamples:
    @responses.activate
    def test_fetches_and_maps_model_metadata(self):
        # Mock the model metadata endpoint
        responses.add(
            responses.GET,
            "https://depmap.org/portal/api/download/all",
            json=[
                {"fileName": "Model.csv", "fileDescription": "Model metadata", "tagsUrl": "model", "size": 1000, "releaseName": "24Q4", "downloadUrl": "https://depmap.org/portal/api/download/Model.csv"},
            ],
            status=200,
        )
        responses.add(
            responses.GET,
            "https://depmap.org/portal/api/download/Model.csv",
            body="ModelID,CellLineName,OncotreeLineage,OncotreeSubtype,OncotreePrimaryDisease,COSMICID\nACH-000001,A549,Lung,Non-Small Cell Lung Cancer,Lung Cancer,905949\nACH-000002,MCF7,Breast,Breast Invasive Ductal Carcinoma,Breast Cancer,905946\n",
            status=200,
        )
        table = fetch_ccle_samples()
        assert table is not None
        assert table["total_samples"] >= 1
        # Check first row has cell_line as sample_type
        cols = table["columns"]
        type_idx = cols.index("sample_type")
        assert table["rows"][0][type_idx] == "cell_line"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `conda run -n dataseek pytest tests/test_search_ccle.py::TestFetchCcleSamples -v`
Expected: FAIL with ImportError

- [ ] **Step 3: Implement fetch_ccle_samples in search_ccle.py**

Add after the existing imports at line 8 of `scripts/search_ccle.py`:

```python
import csv
import io

from scripts.sample_utils import build_sample_table
```

Add after `search_ccle` function (after line 65) in `scripts/search_ccle.py`:

```python
_model_metadata_cache = None


def _fetch_model_metadata() -> list[dict]:
    """Fetch DepMap model metadata CSV. Cached for session."""
    global _model_metadata_cache
    if _model_metadata_cache is not None:
        return _model_metadata_cache
    try:
        resp = fetch_with_retry(DEPMAP_DOWNLOAD_API)
        if resp.status_code != 200:
            return []
        entries = resp.json()
        model_entry = None
        for entry in entries:
            if entry.get("fileName", "").lower() == "model.csv":
                model_entry = entry
                break
        if not model_entry:
            return []
        download_url = model_entry.get("downloadUrl", "")
        if not download_url:
            return []
        resp = fetch_with_retry(download_url, timeout=60)
        if resp.status_code != 200:
            return []
        reader = csv.DictReader(io.StringIO(resp.text))
        _model_metadata_cache = list(reader)
        return _model_metadata_cache
    except Exception as e:
        logger.warning(f"Failed to fetch DepMap model metadata: {e}")
        return []


def fetch_ccle_samples() -> dict | None:
    """Fetch cell line sample metadata from DepMap model metadata. Returns sample_table dict."""
    models = _fetch_model_metadata()
    if not models:
        return None
    rows = []
    for model in models:
        rows.append({
            "sample_id": model.get("ModelID", ""),
            "tissue_site": model.get("OncotreeLineage", "N/A"),
            "sample_type": "cell_line",
            "treatment_status": "N/A",
            "disease": model.get("OncotreePrimaryDisease", "N/A"),
            "cell_type": model.get("CellLineName", "N/A"),
            "age": "N/A",
            "sex": "N/A",
            "stage": "N/A",
            "lineage": model.get("OncotreeLineage", "N/A"),
            "sublineage": model.get("OncotreeSubtype", "N/A"),
            "cosmic_id": model.get("COSMICID", "N/A"),
        })
    if not rows:
        return None
    return build_sample_table(rows, source="CCLE")
```

Modify `search_ccle()` to attach sample_table to each result. Replace lines 57-64:

Replace:
```python
    results = []
    for entry in all_entries:
        entry_tags = entry.get("tagsUrl", "")
        if any(tag in entry_tags for tag in tags):
            result = parse_ccle_dataset(entry, omic)
            if disease and disease.lower() not in result["title"].lower():
                continue
            results.append(result)
    logger.info(f"Found {len(results)} CCLE datasets")
    return results[:max_results]
```

With:
```python
    results = []
    sample_table = fetch_ccle_samples()
    for entry in all_entries:
        entry_tags = entry.get("tagsUrl", "")
        if any(tag in entry_tags for tag in tags):
            result = parse_ccle_dataset(entry, omic)
            if disease and disease.lower() not in result["title"].lower():
                continue
            if sample_table:
                result["sample_table"] = sample_table
            results.append(result)
    logger.info(f"Found {len(results)} CCLE datasets")
    return results[:max_results]
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `conda run -n dataseek pytest tests/test_search_ccle.py -v`
Expected: All tests PASS

- [ ] **Step 5: Commit**

```bash
git add scripts/search_ccle.py tests/test_search_ccle.py
git commit -m "feat: add CCLE cell line sample metadata via DepMap model metadata"
```

---

### Task 6: Add fetch_samples to search_depmap.py (shared model metadata)

**Files:**
- Modify: `scripts/search_depmap.py:1-61`
- Test: `tests/test_search_depmap.py`

- [ ] **Step 1: Write failing tests for DepMap sample fetching**

Add to `tests/test_search_depmap.py`:

```python
from scripts.search_depmap import fetch_depmap_samples


class TestFetchDepmapSamples:
    @responses.activate
    def test_fetches_model_metadata(self):
        responses.add(
            responses.GET,
            "https://depmap.org/portal/api/download/all",
            json=[
                {"fileName": "Model.csv", "fileDescription": "Model metadata", "tagsUrl": "model", "size": 1000, "releaseName": "24Q4", "downloadUrl": "https://depmap.org/portal/api/download/Model.csv"},
            ],
            status=200,
        )
        responses.add(
            responses.GET,
            "https://depmap.org/portal/api/download/Model.csv",
            body="ModelID,CellLineName,OncotreeLineage,OncotreeSubtype,OncotreePrimaryDisease,COSMICID\nACH-000001,A549,Lung,NSCLC,Lung Cancer,905949\n",
            status=200,
        )
        table = fetch_depmap_samples()
        assert table is not None
        assert table["total_samples"] >= 1
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `conda run -n dataseek pytest tests/test_search_depmap.py::TestFetchDepmapSamples -v`
Expected: FAIL with ImportError

- [ ] **Step 3: Implement fetch_depmap_samples in search_depmap.py**

Add after existing imports at line 8 of `scripts/search_depmap.py`:

```python
import csv
import io

from scripts.sample_utils import build_sample_table
```

Add after `search_depmap` function (after line 50) in `scripts/search_depmap.py`:

```python
_model_metadata_cache = None


def _fetch_model_metadata() -> list[dict]:
    """Fetch DepMap model metadata CSV. Cached for session."""
    global _model_metadata_cache
    if _model_metadata_cache is not None:
        return _model_metadata_cache
    try:
        resp = fetch_with_retry(DEPMAP_DOWNLOAD_API)
        if resp.status_code != 200:
            return []
        entries = resp.json()
        model_entry = None
        for entry in entries:
            if entry.get("fileName", "").lower() == "model.csv":
                model_entry = entry
                break
        if not model_entry:
            return []
        download_url = model_entry.get("downloadUrl", "")
        if not download_url:
            return []
        resp = fetch_with_retry(download_url, timeout=60)
        if resp.status_code != 200:
            return []
        reader = csv.DictReader(io.StringIO(resp.text))
        _model_metadata_cache = list(reader)
        return _model_metadata_cache
    except Exception as e:
        logger.warning(f"Failed to fetch DepMap model metadata: {e}")
        return []


def fetch_depmap_samples() -> dict | None:
    """Fetch cell line sample metadata from DepMap model metadata. Returns sample_table dict."""
    models = _fetch_model_metadata()
    if not models:
        return None
    rows = []
    for model in models:
        rows.append({
            "sample_id": model.get("ModelID", ""),
            "tissue_site": model.get("OncotreeLineage", "N/A"),
            "sample_type": "cell_line",
            "treatment_status": "N/A",
            "disease": model.get("OncotreePrimaryDisease", "N/A"),
            "cell_type": model.get("CellLineName", "N/A"),
            "age": "N/A",
            "sex": "N/A",
            "stage": "N/A",
            "lineage": model.get("OncotreeLineage", "N/A"),
            "sublineage": model.get("OncotreeSubtype", "N/A"),
            "cosmic_id": model.get("COSMICID", "N/A"),
        })
    if not rows:
        return None
    return build_sample_table(rows, source="DepMap")
```

Modify `search_depmap()` to attach sample_table. Replace lines 45-50:

Replace:
```python
    results = []
    for entry in all_entries:
        entry_tags = entry.get("tagsUrl", "")
        if any(tag in entry_tags for tag in tags):
            results.append(parse_depmap_dataset(entry, omic))
    logger.info(f"Found {len(results)} DepMap datasets")
    return results[:max_results]
```

With:
```python
    results = []
    sample_table = fetch_depmap_samples()
    for entry in all_entries:
        entry_tags = entry.get("tagsUrl", "")
        if any(tag in entry_tags for tag in tags):
            result = parse_depmap_dataset(entry, omic)
            if sample_table:
                result["sample_table"] = sample_table
            results.append(result)
    logger.info(f"Found {len(results)} DepMap datasets")
    return results[:max_results]
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `conda run -n dataseek pytest tests/test_search_depmap.py -v`
Expected: All tests PASS

- [ ] **Step 5: Commit**

```bash
git add scripts/search_depmap.py tests/test_search_depmap.py
git commit -m "feat: add DepMap cell line sample metadata via model metadata"
```

---

### Task 7: Add fetch_samples to search_xena.py (phenotype matrices)

**Files:**
- Modify: `scripts/search_xena.py:1-105`
- Test: `tests/test_search_xena.py`

- [ ] **Step 1: Write failing tests for Xena sample fetching**

Add to `tests/test_search_xena.py`:

```python
from unittest.mock import patch, MagicMock
from scripts.search_xena import fetch_xena_samples, parse_xena_phenotypes

SAMPLE_PHENOTYPE_DATA = [
    {"sampleID": "TCGA-06-0125-01", "_primary_disease": "GBM", "sample_type": "Primary Tumor", "age_at_initial_pathologic_diagnosis": "55", "gender": "MALE", "_PATIENT": "TCGA-06-0125"},
    {"sampleID": "TCGA-06-0126-01", "_primary_disease": "GBM", "sample_type": "Primary Tumor", "age_at_initial_pathologic_diagnosis": "42", "gender": "FEMALE", "_PATIENT": "TCGA-06-0126"},
    {"sampleID": "TCGA-06-0127-06", "_primary_disease": "GBM", "sample_type": "Metastatic", "age_at_initial_pathologic_diagnosis": "60", "gender": "MALE", "_PATIENT": "TCGA-06-0127"},
]


class TestParseXenaPhenotypes:
    def test_maps_tcga_phenotype_fields(self):
        rows = parse_xena_phenotypes(SAMPLE_PHENOTYPE_DATA, cohort="TCGA GBM")
        assert len(rows) == 3
        assert rows[0]["sample_id"] == "TCGA-06-0125-01"
        assert rows[0]["disease"] == "GBM"
        assert rows[0]["sample_type"] == "primary"
        assert rows[0]["age"] == "55"
        assert rows[0]["sex"] == "MALE"
        assert rows[0]["cohort_name"] == "TCGA GBM"

    def test_maps_metastatic_sample_type(self):
        rows = parse_xena_phenotypes(SAMPLE_PHENOTYPE_DATA, cohort="TCGA GBM")
        assert rows[2]["sample_type"] == "metastasis"
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `conda run -n dataseek pytest tests/test_search_xena.py::TestParseXenaPhenotypes -v`
Expected: FAIL with ImportError

- [ ] **Step 3: Implement parse_xena_phenotypes and fetch_xena_samples**

Add after existing imports at line 12 of `scripts/search_xena.py`:

```python
from scripts.sample_utils import build_sample_table, SAMPLE_CAP
from scripts.supplement_fetch import fetch_supplementary_tables
```

Add after `_search_xena_api` function (after line 92) in `scripts/search_xena.py`:

```python
TCGA_SAMPLE_TYPE_MAP = {
    "primary tumor": "primary",
    "primary solid tumor": "primary",
    "recurrent tumor": "primary",
    "recurrent solid tumor": "primary",
    "metastatic": "metastasis",
    "solid tissue normal": "normal",
    "blood derived normal": "normal",
}


def parse_xena_phenotypes(phenotype_data: list[dict], cohort: str = "") -> list[dict]:
    """Parse Xena phenotype data into sample metadata rows."""
    rows = []
    for entry in phenotype_data[:SAMPLE_CAP]:
        sample_id = entry.get("sampleID", entry.get("sample", ""))
        if not sample_id:
            continue
        raw_type = entry.get("sample_type", entry.get("_sample_type", ""))
        sample_type = TCGA_SAMPLE_TYPE_MAP.get(raw_type.lower().strip(), raw_type) if raw_type else "N/A"
        rows.append({
            "sample_id": sample_id,
            "tissue_site": entry.get("_primary_site", entry.get("primary_site", "N/A")),
            "sample_type": sample_type,
            "treatment_status": entry.get("treatment_outcome_first_course", "N/A"),
            "disease": entry.get("_primary_disease", entry.get("disease", "N/A")),
            "cell_type": "N/A",
            "age": entry.get("age_at_initial_pathologic_diagnosis", entry.get("age", "N/A")),
            "sex": entry.get("gender", entry.get("sex", "N/A")),
            "stage": entry.get("pathologic_stage", entry.get("clinical_stage", "N/A")),
            "cohort_name": cohort,
        })
    return rows


def fetch_xena_samples(hub: str, dataset_name: str, cohort: str = "", paper: dict = None) -> dict | None:
    """Fetch sample phenotype data for a Xena dataset. Returns sample_table dict or None."""
    if xena is None:
        return None
    try:
        samples = xena.dataset_samples(hub, dataset_name, None)
        if not samples:
            return None
        # Fetch phenotype data for these samples
        phenotype_datasets = xena.dataset_phenotypes(hub, dataset_name)
        if not phenotype_datasets:
            return None
        # phenotype_datasets is a list of phenotype dataset names
        pheno_data = []
        for pheno_ds in phenotype_datasets[:1]:  # Use first phenotype dataset
            try:
                pheno_samples = xena.dataset_samples(hub, pheno_ds, None)
                if not pheno_samples:
                    continue
                fields = xena.dataset_field(hub, pheno_ds)
                if not fields:
                    continue
                # Fetch field values for each sample
                for i, sample_id in enumerate(pheno_samples[:SAMPLE_CAP]):
                    entry = {"sampleID": sample_id}
                    for field in fields:
                        values = xena.dataset_probe_values(hub, pheno_ds, [sample_id], [field])
                        if values and values[0]:
                            entry[field] = str(values[0][0]) if values[0][0] is not None else "N/A"
                    pheno_data.append(entry)
            except Exception as e:
                logger.warning(f"Failed to fetch Xena phenotype {pheno_ds}: {e}")
        if not pheno_data:
            return None
        rows = parse_xena_phenotypes(pheno_data, cohort=cohort)
        if not rows:
            return None
        return build_sample_table(rows, source="Xena")
    except Exception as e:
        logger.warning(f"Failed to fetch Xena samples for {dataset_name}: {e}")
        return None
```

Modify `search_xena()` to attach sample_table. In the loop at lines 68-72, replace:

```python
              for ds in filtered:
                  ds["host"] = hub
                  if disease and disease.lower() not in ds.get("cohort", "").lower():
                      continue
                  results.append(parse_xena_dataset(ds, omic))
```

With:
```python
              for ds in filtered:
                  ds["host"] = hub
                  if disease and disease.lower() not in ds.get("cohort", "").lower():
                      continue
                  record = parse_xena_dataset(ds, omic)
                  sample_table = fetch_xena_samples(hub, ds.get("name", ""), cohort=ds.get("cohort", ""))
                  if sample_table:
                      record["sample_table"] = sample_table
                  results.append(record)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `conda run -n dataseek pytest tests/test_search_xena.py -v`
Expected: All tests PASS

- [ ] **Step 5: Commit**

```bash
git add scripts/search_xena.py tests/test_search_xena.py
git commit -m "feat: add Xena phenotype-based sample metadata with TCGA sample type mapping"
```

---

### Task 8: Add fetch_samples to search_scp.py (study-level metadata)

**Files:**
- Modify: `scripts/search_scp.py:1-69`
- Test: `tests/test_search_scp.py`

- [ ] **Step 1: Write failing tests for SCP sample fetching**

Add to `tests/test_search_scp.py`:

```python
from scripts.search_scp import fetch_scp_samples

SAMPLE_SCP_STUDY_DETAIL = {
    "accession": "SCP1234",
    "name": "Single-cell atlas of glioblastoma",
    "description": "Comprehensive single-cell atlas",
    "cell_count": 50000,
    "organ": ["brain"],
    "disease": ["glioblastoma"],
    "library_preparation_protocol": ["10x 3' v3"],
    "species": ["Homo sapiens"],
}


class TestFetchScpSamples:
    @responses.activate
    def test_builds_study_level_sample_table(self):
        responses.add(
            responses.GET,
            "https://singlecell.broadinstitute.org/single_cell/api/v1/studies/SCP1234",
            json=SAMPLE_SCP_STUDY_DETAIL,
            status=200,
        )
        table = fetch_scp_samples("SCP1234")
        assert table is not None
        assert table["total_samples"] == 1
        cols = table["columns"]
        assert "cell_count" in cols
        assert "library_prep" in cols
        # Check tissue_site mapped from organ
        site_idx = cols.index("tissue_site")
        assert table["rows"][0][site_idx] == "brain"

    @responses.activate
    def test_returns_none_on_api_failure(self):
        responses.add(
            responses.GET,
            "https://singlecell.broadinstitute.org/single_cell/api/v1/studies/SCP9999",
            status=404,
        )
        table = fetch_scp_samples("SCP9999")
        assert table is None
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `conda run -n dataseek pytest tests/test_search_scp.py::TestFetchScpSamples -v`
Expected: FAIL with ImportError

- [ ] **Step 3: Implement fetch_scp_samples in search_scp.py**

Add after the existing imports at line 9 of `scripts/search_scp.py`:

```python
from scripts.sample_utils import build_sample_table
```

Add after `search_scp` function (after line 56) in `scripts/search_scp.py`:

```python
def fetch_scp_samples(accession: str) -> dict | None:
    """Fetch study-level sample metadata from SCP. Returns sample_table dict or None.
    SCP doesn't expose per-cell metadata via API, so we create a single study-level row."""
    token = os.environ.get("SCP_TOKEN", "")
    headers = {}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    try:
        resp = fetch_with_retry(
            f"{SCP_API_BASE}/studies/{accession}",
            headers=headers,
            timeout=15,
            max_retries=2,
        )
        if resp.status_code != 200:
            return None
        study = resp.json()
        organs = study.get("organ", [])
        diseases = study.get("disease", [])
        protocols = study.get("library_preparation_protocol", [])
        cell_count = study.get("cell_count", 0)
        rows = [{
            "sample_id": accession,
            "tissue_site": ", ".join(organs) if organs else "N/A",
            "sample_type": "N/A",
            "treatment_status": "N/A",
            "disease": ", ".join(diseases) if diseases else "N/A",
            "cell_type": "N/A",
            "age": "N/A",
            "sex": "N/A",
            "stage": "N/A",
            "cell_count": str(cell_count) if cell_count else "N/A",
            "library_prep": ", ".join(protocols) if protocols else "N/A",
        }]
        return build_sample_table(rows, source="SCP")
    except Exception as e:
        logger.warning(f"Failed to fetch SCP samples for {accession}: {e}")
        return None
```

Modify `search_scp()` to attach sample_table. Replace lines 53-55:

Replace:
```python
    for study in studies:
        results.append(parse_scp_study(study, omic))
```

With:
```python
    for study in studies:
        record = parse_scp_study(study, omic)
        sample_table = fetch_scp_samples(study.get("accession", ""))
        if sample_table:
            record["sample_table"] = sample_table
        results.append(record)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `conda run -n dataseek pytest tests/test_search_scp.py -v`
Expected: All tests PASS

- [ ] **Step 5: Commit**

```bash
git add scripts/search_scp.py tests/test_search_scp.py
git commit -m "feat: add SCP study-level sample metadata"
```

---

### Task 9: Update conftest.py and CLAUDE.md with new schema

**Files:**
- Modify: `tests/conftest.py:41-44`
- Modify: `CLAUDE.md`

- [ ] **Step 1: Update SEARCH_RESULT_REQUIRED_KEYS in conftest.py**

The `sample_table` field is optional (not required), so `SEARCH_RESULT_REQUIRED_KEYS` stays the same. But add a new constant for sample_table validation:

Add after `PAPER_REQUIRED_KEYS` (after line 47) in `tests/conftest.py`:

```python
SAMPLE_TABLE_KEYS = {
    "total_samples", "shown_samples", "capped", "columns", "rows", "summary"
}
```

- [ ] **Step 2: Update sample_search_result fixture to include sample_table**

Add `sample_table` to the fixture in `tests/conftest.py` after line 38 (`"date_submitted": "2025-05-01"`):

```python
        "sample_table": {
            "total_samples": 12,
            "shown_samples": 12,
            "capped": False,
            "columns": ["sample_id", "tissue_site", "sample_type", "treatment_status", "disease", "cell_type", "age", "sex", "stage", "platform_id"],
            "rows": [
                ["GSM000001", "brain", "primary", "untreated", "glioblastoma", "N/A", "55", "M", "IV", "GPL24676"],
            ],
            "summary": {"primary": 1},
        }
```

- [ ] **Step 3: Update CLAUDE.md SearchResult Schema**

Add `sample_table` to the SearchResult schema in CLAUDE.md after `"date_submitted"`:

```json
  "sample_table": {
    "total_samples": 12,
    "shown_samples": 12,
    "capped": false,
    "columns": ["sample_id", "tissue_site", "sample_type", "treatment_status", "disease", "cell_type", "age", "sex", "stage"],
    "rows": [["GSM123456", "brain", "primary", "untreated", "glioblastoma", "N/A", "45", "M", "IV"]],
    "summary": {"primary": 8, "normal": 4, "treated": 3, "untreated": 9}
  }
```

Also add `sample_utils.py` and `supplement_fetch.py` to the Directory Structure section under `scripts/`.

- [ ] **Step 4: Run full test suite**

Run: `conda run -n dataseek pytest tests/ -v`
Expected: All tests PASS

- [ ] **Step 5: Commit**

```bash
git add tests/conftest.py CLAUDE.md
git commit -m "docs: update schema and test fixtures with sample_table field"
```

---

### Task 10: Integration test — end-to-end search with sample tables

**Files:**
- Create: `tests/test_sample_integration.py`

- [ ] **Step 1: Write integration test for GEO search with sample table**

```python
# tests/test_sample_integration.py
"""Integration tests for sample metadata table generation during search."""

import csv
import json
import os

import pytest
import responses

from scripts.search_geo import search_geo
from scripts.sample_utils import write_sample_csv, summarize_samples
from tests.conftest import SAMPLE_TABLE_KEYS


SAMPLE_GEO_ESEARCH = '<?xml version="1.0"?><eSearchResult><Count>1</Count><RetMax>1</RetMax><IdList><Id>200123456</Id></IdList></eSearchResult>'

SAMPLE_GEO_ESUMMARY = """<?xml version="1.0"?>
<eSummaryResult>
    <DocSum>
        <Id>200123456</Id>
        <Item Name="Accession" Type="String">GSE123456</Item>
        <Item Name="title" Type="String">Test study</Item>
        <Item Name="summary" Type="String">Test summary</Item>
        <Item Name="GPL" Type="String">GPL24676</Item>
        <Item Name="taxon" Type="String">Homo sapiens</Item>
        <Item Name="n_samples" Type="Integer">2</Item>
        <Item Name="PDAT" Type="String">2025/05/01</Item>
        <Item Name="PubMedIds" Type="List"></Item>
        <Item Name="suppFile" Type="String">H5</Item>
    </DocSum>
</eSummaryResult>"""

SAMPLE_SOFT = """\
^SERIES = GSE123456
^SAMPLE = GSM000001
!Sample_characteristics_ch1 = tissue: brain
!Sample_characteristics_ch1 = tumor type: primary
!Sample_characteristics_ch1 = treatment: untreated
!Sample_platform_id = GPL24676
^SAMPLE = GSM000002
!Sample_characteristics_ch1 = tissue: brain
!Sample_characteristics_ch1 = tumor type: metastasis
!Sample_characteristics_ch1 = treatment: cisplatin
!Sample_platform_id = GPL24676
"""


class TestGeoSearchWithSampleTable:
    @responses.activate
    def test_search_result_includes_sample_table(self):
        responses.add(responses.GET, "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", body=SAMPLE_GEO_ESEARCH, status=200)
        responses.add(responses.GET, "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", body=SAMPLE_GEO_ESUMMARY, status=200)
        responses.add(responses.GET, "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi", body=SAMPLE_SOFT, status=200)
        # Block supplementary fetch (no PMID)
        results = search_geo(omic="scRNAseq", organism="human", max_results=5)
        assert len(results) >= 1
        result = results[0]
        assert "sample_table" in result
        table = result["sample_table"]
        for key in SAMPLE_TABLE_KEYS:
            assert key in table
        assert table["total_samples"] == 2
        assert len(table["rows"]) == 2

    @responses.activate
    def test_sample_csv_write_from_search_result(self, tmp_path):
        responses.add(responses.GET, "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", body=SAMPLE_GEO_ESEARCH, status=200)
        responses.add(responses.GET, "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", body=SAMPLE_GEO_ESUMMARY, status=200)
        responses.add(responses.GET, "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi", body=SAMPLE_SOFT, status=200)
        results = search_geo(omic="scRNAseq", organism="human", max_results=5)
        table = results[0]["sample_table"]
        filepath = write_sample_csv("GSE123456", table, str(tmp_path))
        assert os.path.exists(filepath)
        with open(filepath) as f:
            reader = csv.reader(f)
            header = next(reader)
            assert "sample_id" in header
            assert "tissue_site" in header
            rows = list(reader)
            assert len(rows) == 2

    @responses.activate
    def test_sample_summary_from_search_result(self):
        responses.add(responses.GET, "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", body=SAMPLE_GEO_ESEARCH, status=200)
        responses.add(responses.GET, "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", body=SAMPLE_GEO_ESUMMARY, status=200)
        responses.add(responses.GET, "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi", body=SAMPLE_SOFT, status=200)
        results = search_geo(omic="scRNAseq", organism="human", max_results=5)
        table = results[0]["sample_table"]
        summary = summarize_samples(table)
        assert "Samples (2)" in summary
```

- [ ] **Step 2: Run integration tests**

Run: `conda run -n dataseek pytest tests/test_sample_integration.py -v`
Expected: All 3 tests PASS

- [ ] **Step 3: Run full test suite**

Run: `conda run -n dataseek pytest tests/ -v`
Expected: All tests PASS

- [ ] **Step 4: Commit**

```bash
git add tests/test_sample_integration.py
git commit -m "test: add integration tests for sample metadata table generation during search"
```

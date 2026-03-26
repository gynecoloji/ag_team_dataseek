"""Fetch and parse publication supplementary tables for sample-level clinical metadata."""

import os
import tempfile
import xml.etree.ElementTree as ET

import pandas as pd
import requests

from scripts.utils import fetch_with_retry, setup_logger

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

IDCONV_URL = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
PMC_OA_URL = "https://www.ncbi.nlm.nih.gov/pmc/utils/oa.cgi"
TABULAR_EXTENSIONS = {".csv", ".tsv", ".xlsx", ".xls"}
MAX_SUPPLEMENT_FILES = 5
MAX_FILE_SIZE_MB = 10

SAMPLE_ID_KEYWORDS = {
    "sample",
    "patient",
    "subject",
    "barcode",
    "sample_id",
    "patient_id",
    "case_id",
    "specimen",
}

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

logger = setup_logger("supplement_fetch")


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------


def resolve_pmcid(pmid: str | None = None, doi: str | None = None) -> str | None:
    """Convert a PMID or DOI to a PMCID using the NCBI ID Converter API.

    Returns the PMCID string (e.g. "PMC10500001") or None if not found.
    """
    if pmid is None and doi is None:
        return None

    params: dict = {"tool": "dataseek", "email": "dataseek@example.com", "format": "xml"}
    if pmid is not None:
        params["ids"] = str(pmid)
    else:
        params["ids"] = str(doi)

    try:
        resp = fetch_with_retry(IDCONV_URL, params=params, max_retries=2, base_delay=0.5)
        if resp.status_code != 200:
            logger.warning("IDCONV returned status %s", resp.status_code)
            return None
        root = ET.fromstring(resp.text)
        record = root.find("record")
        if record is None:
            return None
        pmcid = record.get("pmcid")
        return pmcid if pmcid else None
    except Exception as exc:
        logger.warning("resolve_pmcid failed: %s", exc)
        return None


def fetch_pmc_supplement_list(pmcid: str) -> list[dict]:
    """Query the PMC OA API for supplementary file links.

    Returns a list of dicts with keys ``url`` and ``filename``, capped at
    MAX_SUPPLEMENT_FILES entries.
    """
    params = {"id": pmcid}
    try:
        resp = fetch_with_retry(PMC_OA_URL, params=params, max_retries=2, base_delay=0.5)
        if resp.status_code != 200:
            logger.warning("PMC OA returned status %s for %s", resp.status_code, pmcid)
            return []
        root = ET.fromstring(resp.text)
        supplements = []
        for link in root.iter("link"):
            href = link.get("href", "")
            format_ = link.get("format", "").lower()
            # Collect supplementary files only
            if not href:
                continue
            filename = href.split("/")[-1]
            supplements.append({"url": href, "filename": filename})
            if len(supplements) >= MAX_SUPPLEMENT_FILES:
                break
        return supplements
    except Exception as exc:
        logger.warning("fetch_pmc_supplement_list failed for %s: %s", pmcid, exc)
        return []


def parse_tabular_file(filepath: str) -> pd.DataFrame | None:
    """Read a CSV, TSV, or Excel file into a DataFrame (max 500 rows).

    Returns None on any failure (file missing, parse error, unsupported format).
    """
    if not os.path.exists(filepath):
        return None
    ext = os.path.splitext(filepath)[1].lower()
    if ext not in TABULAR_EXTENSIONS:
        return None
    try:
        if ext == ".csv":
            return pd.read_csv(filepath, nrows=500)
        elif ext == ".tsv":
            return pd.read_csv(filepath, sep="\t", nrows=500)
        elif ext in {".xlsx", ".xls"}:
            return pd.read_excel(filepath, nrows=500)
    except Exception as exc:
        logger.warning("parse_tabular_file failed for %s: %s", filepath, exc)
    return None


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _column_matches_keywords(col_name: str, keywords: set) -> bool:
    """Return True if col_name (lowercased) contains any keyword from the set."""
    col_lower = col_name.lower()
    return any(kw in col_lower for kw in keywords)


def _find_sample_id_column(df: pd.DataFrame) -> str | None:
    """Return the name of the first column that matches SAMPLE_ID_KEYWORDS, or None."""
    for col in df.columns:
        if _column_matches_keywords(col, SAMPLE_ID_KEYWORDS):
            return col
    return None


def _map_columns_to_schema(df: pd.DataFrame) -> dict[str, str]:
    """Map DataFrame column names to standard schema field names.

    Returns a dict of {schema_field: df_column_name} for matched columns.
    """
    mapping: dict[str, str] = {}
    for schema_field, keywords in CLINICAL_KEYWORDS.items():
        for col in df.columns:
            if _column_matches_keywords(col, keywords):
                mapping[schema_field] = col
                break
    return mapping


def match_supplement_to_samples(
    df: pd.DataFrame, sample_ids: list[str] | None
) -> list[dict] | None:
    """Cross-reference a DataFrame with a list of sample IDs.

    Steps:
    1. Find the sample-ID column; return None if absent.
    2. Map columns to schema; return None if no clinical columns found.
    3. Filter rows matching sample_ids (or all rows if sample_ids is None/empty).
    4. Return enriched rows as a list of dicts with standardised keys.
    """
    sample_id_col = _find_sample_id_column(df)
    if sample_id_col is None:
        return None

    schema_mapping = _map_columns_to_schema(df)
    if not schema_mapping:
        return None

    # Filter to matching sample IDs when provided
    if sample_ids:
        mask = df[sample_id_col].astype(str).isin([str(s) for s in sample_ids])
        subset = df[mask]
    else:
        subset = df

    rows = []
    for _, row in subset.iterrows():
        record: dict = {"sample_id": str(row[sample_id_col])}
        for schema_field, col in schema_mapping.items():
            val = row.get(col)
            if val is not None and str(val).strip() != "":
                record[schema_field] = str(val)
        rows.append(record)

    return rows if rows else None


def fetch_supplementary_tables(
    doi: str | None = None,
    pmid: str | None = None,
    sample_ids: list[str] | None = None,
) -> list[dict] | None:
    """Full pipeline: resolve PMCID -> fetch supplements -> parse -> match.

    Returns enriched sample rows or None on any failure.
    """
    if doi is None and pmid is None:
        return None

    pmcid = resolve_pmcid(pmid=pmid, doi=doi)
    if pmcid is None:
        logger.info("Could not resolve PMCID for doi=%s pmid=%s", doi, pmid)
        return None

    supplement_list = fetch_pmc_supplement_list(pmcid)
    if not supplement_list:
        logger.info("No tabular supplements found for %s", pmcid)
        return None

    all_rows: list[dict] = []
    for supp in supplement_list:
        url = supp["url"]
        filename = supp["filename"]
        ext = os.path.splitext(filename)[1].lower()
        if ext not in TABULAR_EXTENSIONS:
            continue

        # Download to a temp file
        try:
            resp = fetch_with_retry(url, timeout=15)
            resp.raise_for_status()
            size_mb = len(resp.content) / (1024 * 1024)
            if size_mb > MAX_FILE_SIZE_MB:
                logger.info("Skipping %s: size %.1f MB exceeds limit", filename, size_mb)
                continue

            with tempfile.NamedTemporaryFile(suffix=ext, delete=False) as tmp:
                tmp.write(resp.content)
                tmp_path = tmp.name

            df = parse_tabular_file(tmp_path)
            os.unlink(tmp_path)

            if df is None:
                continue

            matched = match_supplement_to_samples(df, sample_ids)
            if matched:
                all_rows.extend(matched)

        except Exception as exc:
            logger.warning("Failed to fetch/parse %s: %s", filename, exc)
            continue

    return all_rows if all_rows else None

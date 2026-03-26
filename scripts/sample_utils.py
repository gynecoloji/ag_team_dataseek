"""Shared utilities for sample-level metadata: schema, normalization, CSV output, summaries."""

import csv
import os
from collections import Counter

# ---------------------------------------------------------------------------
# Schema constants
# ---------------------------------------------------------------------------

STANDARD_COLUMNS = [
    "sample_id",
    "tissue_site",
    "sample_type",
    "treatment_status",
    "disease",
    "cell_type",
    "age",
    "sex",
    "stage",
]

SOURCE_EXTRA_COLUMNS = {
    "GEO": ["platform_id"],
    "CCLE": ["lineage", "sublineage", "cosmic_id"],
    "DepMap": ["lineage", "sublineage", "cosmic_id"],
    "Xena": ["cohort_name"],
    "SCP": ["cell_count", "library_prep"],
}

SAMPLE_CAP = 500


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _build_summary_dict(normalized: list[dict]) -> dict:
    """Build summary counts from sample_type and treatment_status fields.

    Returns a dict with two keys:
        ``{"sample_type": {value: count, ...}, "treatment_status": {value: count, ...}}``

    N/A values are excluded from counts.  Empty input returns the two-key
    structure with empty inner dicts.
    """
    sample_type_counts: dict[str, int] = {}
    treatment_status_counts: dict[str, int] = {}

    for row in normalized:
        st = row.get("sample_type", "N/A")
        if st and st != "N/A":
            sample_type_counts[st] = sample_type_counts.get(st, 0) + 1

        ts = row.get("treatment_status", "N/A")
        if ts and ts != "N/A":
            treatment_status_counts[ts] = treatment_status_counts.get(ts, 0) + 1

    return {
        "sample_type": sample_type_counts,
        "treatment_status": treatment_status_counts,
    }


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------

def normalize_sample_table(rows: list[dict], source: str) -> list[dict]:
    """Normalize rows to standard schema + source extras, fill missing with 'N/A', cap at SAMPLE_CAP.

    Args:
        rows: Raw sample dicts from a data source.
        source: Source name (e.g. "GEO", "CCLE").

    Returns:
        Normalized list of dicts, at most SAMPLE_CAP entries.
    """
    extra_cols = SOURCE_EXTRA_COLUMNS.get(source, [])
    all_cols = STANDARD_COLUMNS + extra_cols

    normalized = []
    for row in rows[:SAMPLE_CAP]:
        norm = {col: row.get(col, "N/A") or "N/A" for col in all_cols}
        normalized.append(norm)

    return normalized


def write_sample_csv(accession: str, table_data: dict, output_dir: str) -> str:
    """Write sample table to CSV. If capped, appends a comment at end.

    Args:
        accession: Dataset accession ID used to name the file.
        table_data: Dict as produced by build_sample_table.
        output_dir: Directory to write the CSV into (created if missing).

    Returns:
        Absolute path to the written CSV file.
    """
    os.makedirs(output_dir, exist_ok=True)
    filepath = os.path.join(output_dir, f"{accession}_samples.csv")

    columns = table_data["columns"]
    rows = table_data["rows"]
    capped = table_data.get("capped", False)
    total_samples = table_data.get("total_samples", len(rows))
    shown_samples = table_data.get("shown_samples", len(rows))

    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        for row in rows:
            writer.writerow(row)
        if capped:
            hidden = total_samples - shown_samples
            f.write(f"\n# {hidden} additional samples not shown\n")

    return filepath


def summarize_samples(table_data: dict) -> str:
    """Generate a condensed one-line summary of sample metadata.

    Format: ``Samples (N): A primary, B normal | untreated: X, treated: Y``

    Args:
        table_data: Dict as produced by build_sample_table.

    Returns:
        Human-readable summary string.
    """
    shown = table_data.get("shown_samples", 0)
    summary_dict = table_data.get("summary", {})

    sample_type_counts = summary_dict.get("sample_type", {})
    treatment_status_counts = summary_dict.get("treatment_status", {})

    type_parts = [f"{count} {key}" for key, count in sample_type_counts.items()]
    treatment_parts = [f"{key}: {count}" for key, count in treatment_status_counts.items()]

    type_str = ", ".join(type_parts) if type_parts else "N/A"
    treatment_str = ", ".join(treatment_parts) if treatment_parts else "N/A"

    return f"Samples ({shown}): {type_str} | {treatment_str}"


def build_sample_table(rows: list[dict], source: str, total_count: int | None = None) -> dict:
    """Build sample_table dict with normalized rows, columns, and summary.

    Args:
        rows: Raw sample dicts from a data source.
        source: Source name (e.g. "GEO", "CCLE").
        total_count: True total sample count from the source (may exceed len(rows) when
                     the caller already pre-sliced). If None, defaults to len(rows).

    Returns:
        Dict with keys: total_samples, shown_samples, capped, columns, rows, summary.
    """
    normalized = normalize_sample_table(rows, source)

    extra_cols = SOURCE_EXTRA_COLUMNS.get(source, [])
    columns = STANDARD_COLUMNS + extra_cols

    csv_rows = [[r[c] for c in columns] for r in normalized]

    raw_total = total_count if total_count is not None else len(rows)
    shown = len(normalized)
    capped = raw_total > SAMPLE_CAP or len(rows) > SAMPLE_CAP

    summary = _build_summary_dict(normalized)

    return {
        "total_samples": raw_total,
        "shown_samples": shown,
        "capped": capped,
        "columns": columns,
        "rows": csv_rows,
        "summary": summary,
    }

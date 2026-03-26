#!/usr/bin/env python3
"""Search Cancer Cell Line Encyclopedia (CCLE) via DepMap Portal API."""

import argparse
import csv
import io
import json
import sys

from scripts.sample_utils import build_sample_table
from scripts.utils import fetch_with_retry, setup_logger

DEPMAP_DOWNLOAD_API = "https://depmap.org/portal/api/download/all"

OMIC_TO_CCLE_TAGS = {
    "bulkRNAseq": ["expression"],
    "bulkGenomicseq": ["copy_number", "mutation", "wgs", "wes"],
}

logger = setup_logger("search_ccle")

_model_metadata_cache = None


def _fetch_model_metadata() -> list[dict]:
    """Fetch DepMap model metadata CSV and return parsed rows.

    Results are cached in the module-level ``_model_metadata_cache`` variable
    so subsequent calls within a session do not repeat the HTTP requests.
    """
    global _model_metadata_cache
    if _model_metadata_cache is not None:
        return _model_metadata_cache

    resp = fetch_with_retry(DEPMAP_DOWNLOAD_API)
    if resp.status_code != 200:
        logger.warning(f"DepMap API returned status {resp.status_code} while fetching model metadata")
        _model_metadata_cache = []
        return _model_metadata_cache

    entries = resp.json()
    model_entry = None
    for entry in entries:
        if entry.get("fileName", "").lower() == "model.csv":
            model_entry = entry
            break

    if model_entry is None:
        logger.warning("Model.csv not found in DepMap download list")
        _model_metadata_cache = []
        return _model_metadata_cache

    download_url = model_entry.get("downloadUrl", "")
    if not download_url:
        logger.warning("Model.csv entry has no downloadUrl")
        _model_metadata_cache = []
        return _model_metadata_cache

    csv_resp = fetch_with_retry(download_url)
    if csv_resp.status_code != 200:
        logger.warning(f"Failed to download Model.csv: status {csv_resp.status_code}")
        _model_metadata_cache = []
        return _model_metadata_cache

    reader = csv.DictReader(io.StringIO(csv_resp.text))
    _model_metadata_cache = list(reader)
    logger.info(f"Fetched {len(_model_metadata_cache)} model metadata rows from DepMap")
    return _model_metadata_cache


def fetch_ccle_samples() -> dict | None:
    """Map DepMap model metadata to a sample table.

    Returns:
        Sample table dict as produced by ``build_sample_table``, or ``None``
        if model metadata could not be retrieved.
    """
    models = _fetch_model_metadata()
    if not models:
        return None

    rows = []
    for model in models:
        rows.append({
            "sample_id": model.get("ModelID", "N/A"),
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

    return build_sample_table(rows, source="CCLE")


def parse_ccle_dataset(entry, omic):
    size_mb = round(entry.get("size", 0) / 1_000_000, 1)
    filename = entry.get("fileName", "")
    release = entry.get("releaseName", "")
    return {
        "accession": f"CCLE_{filename.replace('.csv', '').replace('.', '_')}",
        "source": "CCLE",
        "title": entry.get("fileDescription", filename),
        "organism": "Homo sapiens",
        "omic_type": omic,
        "platform": "multiple (cell line panel)",
        "disease": "cancer cell lines",
        "tissue": "multiple",
        "sample_count": 0,
        "condition_groups": [],
        "data_files": [{"type": "processed_matrix", "format": "csv", "size_mb": size_mb}],
        "paper": {
            "title": "Cancer Cell Line Encyclopedia",
            "authors": "Broad Institute",
            "journal": "DepMap Portal",
            "doi": "",
            "date": "",
            "abstract": f"Release: {release}. {entry.get('fileDescription', '')}",
        },
        "metadata_quality": "good",
        "date_submitted": "",
    }


def search_ccle(omic, disease=None, max_results=50):
    tags = OMIC_TO_CCLE_TAGS.get(omic, [])
    if not tags:
        logger.info(f"No CCLE mapping for omic type: {omic}")
        return []
    logger.info(f"Querying DepMap download API for tags: {tags}")
    resp = fetch_with_retry(DEPMAP_DOWNLOAD_API)
    if resp.status_code != 200:
        logger.warning(f"DepMap API returned status {resp.status_code}")
        return []
    all_entries = resp.json()
    sample_table = fetch_ccle_samples()
    results = []
    for entry in all_entries:
        entry_tags = entry.get("tagsUrl", "")
        if any(tag in entry_tags for tag in tags):
            result = parse_ccle_dataset(entry, omic)
            if disease and disease.lower() not in result["title"].lower():
                continue
            result["sample_table"] = sample_table
            results.append(result)
    logger.info(f"Found {len(results)} CCLE datasets")
    return results[:max_results]


def main():
    parser = argparse.ArgumentParser(description="Search CCLE via DepMap Portal")
    parser.add_argument("--omic", required=True)
    parser.add_argument("--disease", default=None)
    parser.add_argument("--max-results", type=int, default=50)
    args = parser.parse_args()
    results = search_ccle(omic=args.omic, disease=args.disease, max_results=args.max_results)
    json.dump(results, sys.stdout, indent=2)


if __name__ == "__main__":
    main()

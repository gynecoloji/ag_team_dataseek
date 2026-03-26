#!/usr/bin/env python3
"""Search Single Cell Portal (Broad Institute) for single-cell datasets."""

import argparse
import json
import os
import sys

from scripts.sample_utils import build_sample_table
from scripts.utils import fetch_with_retry, setup_logger

SCP_API_BASE = "https://singlecell.broadinstitute.org/single_cell/api/v1"

logger = setup_logger("search_scp")


def parse_scp_study(study, omic):
    species = study.get("species", [])
    organism = species[0] if species else ""
    organs = study.get("organ", [])
    diseases = study.get("disease", [])
    protocols = study.get("library_preparation_protocol", [])
    return {
        "accession": study.get("accession", ""),
        "source": "SCP",
        "title": study.get("name", ""),
        "organism": organism,
        "omic_type": omic,
        "platform": ", ".join(protocols) if protocols else "unknown",
        "disease": ", ".join(diseases) if diseases else "",
        "tissue": ", ".join(organs) if organs else "",
        "sample_count": 0,
        "condition_groups": [],
        "data_files": [{"type": "study_files", "format": "mixed", "size_mb": 0}],
        "paper": {"title": "", "authors": "", "journal": "", "doi": "", "date": "", "abstract": study.get("description", "")},
        "metadata_quality": "good" if study.get("description") else "partial",
        "date_submitted": "",
    }


def fetch_scp_samples(accession: str) -> dict | None:
    """Fetch study-level metadata from SCP API and return a sample table.

    Args:
        accession: SCP study accession (e.g. "SCP1234").

    Returns:
        Sample table dict from build_sample_table, or None on HTTP error/exception.
    """
    token = os.environ.get("SCP_TOKEN", "")
    headers = {}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    try:
        resp = fetch_with_retry(f"{SCP_API_BASE}/studies/{accession}", headers=headers)
        if resp.status_code != 200:
            logger.warning(f"SCP studies API returned status {resp.status_code} for {accession}")
            return None
        study = resp.json()
    except Exception as exc:
        logger.warning(f"Failed to fetch SCP study detail for {accession}: {exc}")
        return None

    organs = study.get("organ", [])
    diseases = study.get("disease", [])
    protocols = study.get("library_preparation_protocol", [])
    cell_count = study.get("cell_count")

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


def search_scp(omic, organism="human", disease=None, tissue=None, max_results=50):
    token = os.environ.get("SCP_TOKEN", "")
    headers = {}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    terms = []
    organism_map = {"human": "Homo sapiens", "mouse": "Mus musculus"}
    org = organism_map.get(organism.lower(), organism)
    terms.append(org)
    if disease:
        terms.append(disease)
    if tissue:
        terms.append(tissue)
    params = {"terms": " ".join(terms), "limit": max_results}
    logger.info(f"Searching SCP: {params}")
    resp = fetch_with_retry(f"{SCP_API_BASE}/search", params=params, headers=headers)
    if resp.status_code != 200:
        logger.warning(f"SCP API returned status {resp.status_code}")
        return []
    data = resp.json()
    studies = data.get("studies", [])
    results = []
    for study in studies:
        record = parse_scp_study(study, omic)
        accession = record.get("accession", "")
        if accession:
            scp_sample_table = fetch_scp_samples(accession)
            if scp_sample_table:
                record["sample_table"] = scp_sample_table
        results.append(record)
    logger.info(f"Found {len(results)} SCP studies")
    return results


def main():
    parser = argparse.ArgumentParser(description="Search Single Cell Portal")
    parser.add_argument("--omic", required=True)
    parser.add_argument("--organism", default="human")
    parser.add_argument("--disease", default=None)
    parser.add_argument("--tissue", default=None)
    parser.add_argument("--max-results", type=int, default=50)
    args = parser.parse_args()
    results = search_scp(omic=args.omic, organism=args.organism, disease=args.disease, tissue=args.tissue, max_results=args.max_results)
    json.dump(results, sys.stdout, indent=2)


if __name__ == "__main__":
    main()

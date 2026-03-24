#!/usr/bin/env python3
"""Search DepMap portal for functional genomics datasets."""

import argparse
import json
import sys

from scripts.utils import fetch_with_retry, setup_logger

DEPMAP_DOWNLOAD_API = "https://depmap.org/portal/api/download/all"

OMIC_TO_DEPMAP_TAGS = {
    "CRISPR": ["crispr"],
    "bulkGenomicseq": ["copy_number", "mutation", "wgs", "wes", "fusion"],
}

logger = setup_logger("search_depmap")


def parse_depmap_dataset(entry, omic):
    size_mb = round(entry.get("size", 0) / 1_000_000, 1)
    filename = entry.get("fileName", "")
    release = entry.get("releaseName", "")
    return {
        "accession": f"DepMap_{filename.replace('.csv', '').replace('.', '_')}",
        "source": "DepMap",
        "title": entry.get("fileDescription", filename),
        "organism": "Homo sapiens",
        "omic_type": omic,
        "platform": "DepMap (Broad Institute)",
        "disease": "cancer cell lines",
        "tissue": "multiple",
        "sample_count": 0,
        "condition_groups": [],
        "data_files": [{"type": "processed_matrix", "format": "csv", "size_mb": size_mb}],
        "paper": {"title": "Cancer Dependency Map", "authors": "Broad Institute", "journal": "DepMap Portal", "doi": "", "date": "", "abstract": f"Release: {release}. {entry.get('fileDescription', '')}"},
        "metadata_quality": "good",
        "date_submitted": "",
    }


def search_depmap(omic, disease=None, max_results=50):
    tags = OMIC_TO_DEPMAP_TAGS.get(omic, [])
    if not tags:
        logger.info(f"No DepMap mapping for omic type: {omic}")
        return []
    logger.info(f"Querying DepMap API for tags: {tags}")
    resp = fetch_with_retry(DEPMAP_DOWNLOAD_API)
    if resp.status_code != 200:
        return []
    all_entries = resp.json()
    results = []
    for entry in all_entries:
        entry_tags = entry.get("tagsUrl", "")
        if any(tag in entry_tags for tag in tags):
            results.append(parse_depmap_dataset(entry, omic))
    logger.info(f"Found {len(results)} DepMap datasets")
    return results[:max_results]


def main():
    parser = argparse.ArgumentParser(description="Search DepMap for omics datasets")
    parser.add_argument("--omic", required=True)
    parser.add_argument("--disease", default=None)
    parser.add_argument("--max-results", type=int, default=50)
    args = parser.parse_args()
    results = search_depmap(omic=args.omic, disease=args.disease, max_results=args.max_results)
    json.dump(results, sys.stdout, indent=2)


if __name__ == "__main__":
    main()

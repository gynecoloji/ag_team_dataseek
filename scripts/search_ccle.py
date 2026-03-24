#!/usr/bin/env python3
"""Search Cancer Cell Line Encyclopedia (CCLE) via DepMap Portal API."""

import argparse
import json
import sys

from scripts.utils import fetch_with_retry, setup_logger

DEPMAP_DOWNLOAD_API = "https://depmap.org/portal/api/download/all"

OMIC_TO_CCLE_TAGS = {
    "bulkRNAseq": ["expression"],
    "bulkGenomicseq": ["copy_number", "mutation", "wgs", "wes"],
}

logger = setup_logger("search_ccle")


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

#!/usr/bin/env python3
"""Download DepMap datasets."""

import argparse
import json
import os
import sys

from scripts.utils import fetch_with_retry, download_file, setup_logger

DEPMAP_DOWNLOAD_API = "https://depmap.org/portal/api/download/all"

logger = setup_logger("download_depmap")


def download_depmap(accession, output_dir, omic_type="unknown"):
    os.makedirs(output_dir, exist_ok=True)
    data_dir = os.path.join(output_dir, "data")
    meta_dir = os.path.join(output_dir, "metadata")
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(meta_dir, exist_ok=True)

    filename_stem = accession.replace("DepMap_", "")
    logger.info(f"Downloading DepMap dataset: {accession}")

    resp = fetch_with_retry(DEPMAP_DOWNLOAD_API)
    entries = resp.json() if resp.status_code == 200 else []

    downloaded_files = []
    for entry in entries:
        entry_stem = entry.get("fileName", "").replace(".csv", "").replace(".", "_")
        if entry_stem == filename_stem:
            url = entry.get("downloadUrl", "")
            if url:
                dest = os.path.join(data_dir, entry["fileName"])
                try:
                    download_file(url, dest)
                    size_mb = round(os.path.getsize(dest) / 1_000_000, 1)
                    downloaded_files.append({"name": entry["fileName"], "size_mb": size_mb})
                except Exception as e:
                    logger.warning(f"Download failed: {e}")

    meta = {"accession": accession, "source": "DepMap", "organism": "Homo sapiens", "disease": "cancer cell lines"}
    with open(os.path.join(meta_dir, "study_metadata.json"), "w") as f:
        json.dump(meta, f, indent=2)
    with open(os.path.join(meta_dir, "paper.json"), "w") as f:
        json.dump({"title": "DepMap", "authors": "Broad Institute", "journal": "DepMap Portal", "doi": "", "date": "", "abstract": ""}, f, indent=2)

    return {"accession": accession, "output_dir": output_dir, "files_downloaded": len(downloaded_files), "data_files": downloaded_files}


def main():
    parser = argparse.ArgumentParser(description="Download DepMap dataset")
    parser.add_argument("--accession", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--omic-type", default="unknown")
    args = parser.parse_args()
    result = download_depmap(args.accession, args.output_dir, args.omic_type)
    json.dump(result, sys.stdout, indent=2)


if __name__ == "__main__":
    main()

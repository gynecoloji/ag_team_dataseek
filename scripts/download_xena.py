#!/usr/bin/env python3
"""Download UCSC Xena datasets."""

import argparse
import json
import os
import sys

from scripts.utils import download_file, setup_logger

logger = setup_logger("download_xena")


def build_xena_download_url(host, dataset_name):
    return f"{host}/download/{dataset_name}"


def download_xena(accession, output_dir, omic_type="unknown"):
    os.makedirs(output_dir, exist_ok=True)
    data_dir = os.path.join(output_dir, "data")
    meta_dir = os.path.join(output_dir, "metadata")
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(meta_dir, exist_ok=True)

    dataset_name = accession.replace("XENA_", "").replace("_", ".")
    logger.info(f"Downloading Xena dataset: {dataset_name}")

    hubs = ["https://tcga.xenahubs.net", "https://ucscpublic.xenahubs.net", "https://toil.xenahubs.net"]
    downloaded_files = []
    for hub in hubs:
        url = build_xena_download_url(hub, dataset_name)
        dest = os.path.join(data_dir, os.path.basename(dataset_name))
        try:
            download_file(url, dest)
            if os.path.exists(dest) and os.path.getsize(dest) > 0:
                size_mb = round(os.path.getsize(dest) / 1_000_000, 1)
                downloaded_files.append({"name": os.path.basename(dataset_name), "size_mb": size_mb})
                break
        except Exception:
            continue

    meta = {"accession": accession, "source": "Xena", "organism": "Homo sapiens"}
    with open(os.path.join(meta_dir, "study_metadata.json"), "w") as f:
        json.dump(meta, f, indent=2)
    with open(os.path.join(meta_dir, "paper.json"), "w") as f:
        json.dump({"title": "", "authors": "", "journal": "", "doi": "", "date": "", "abstract": ""}, f, indent=2)

    return {"accession": accession, "output_dir": output_dir, "files_downloaded": len(downloaded_files), "data_files": downloaded_files}


def main():
    parser = argparse.ArgumentParser(description="Download Xena dataset")
    parser.add_argument("--accession", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--omic-type", default="unknown")
    args = parser.parse_args()
    result = download_xena(args.accession, args.output_dir, args.omic_type)
    json.dump(result, sys.stdout, indent=2)


if __name__ == "__main__":
    main()

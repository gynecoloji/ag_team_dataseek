#!/usr/bin/env python3
"""Download Single Cell Portal study files."""

import argparse
import json
import os
import sys

from scripts.utils import fetch_with_retry, download_file, setup_logger

SCP_API_BASE = "https://singlecell.broadinstitute.org/single_cell/api/v1"

logger = setup_logger("download_scp")


def download_scp(accession, output_dir, omic_type="unknown"):
    os.makedirs(output_dir, exist_ok=True)
    data_dir = os.path.join(output_dir, "data")
    meta_dir = os.path.join(output_dir, "metadata")
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(meta_dir, exist_ok=True)

    logger.info(f"Downloading SCP study: {accession}")

    token = os.environ.get("SCP_TOKEN", "")
    headers = {"Authorization": f"Bearer {token}"} if token else {}

    resp = fetch_with_retry(f"{SCP_API_BASE}/studies/{accession}", headers=headers)
    if resp.status_code != 200:
        logger.warning(f"SCP API returned {resp.status_code} for {accession}")
        return {"accession": accession, "output_dir": output_dir, "files_downloaded": 0, "data_files": []}

    study = resp.json()

    meta = {"accession": accession, "source": "SCP", "name": study.get("name", ""), "description": study.get("description", "")}
    with open(os.path.join(meta_dir, "study_metadata.json"), "w") as f:
        json.dump(meta, f, indent=2)
    with open(os.path.join(meta_dir, "paper.json"), "w") as f:
        json.dump({"title": "", "authors": "", "journal": "", "doi": "", "date": "", "abstract": study.get("description", "")}, f, indent=2)

    downloaded_files = []
    for sf in study.get("study_files", []):
        url = sf.get("download_url", "")
        name = sf.get("name", "unknown_file")
        if url:
            dest = os.path.join(data_dir, name)
            try:
                download_file(url, dest)
                size_mb = round(os.path.getsize(dest) / 1_000_000, 1)
                downloaded_files.append({"name": name, "size_mb": size_mb})
            except Exception as e:
                logger.warning(f"Failed to download {name}: {e}")

    with open(os.path.join(meta_dir, "sample_metadata.csv"), "w") as f:
        f.write("sample_id,condition,batch\n")

    return {"accession": accession, "output_dir": output_dir, "files_downloaded": len(downloaded_files), "data_files": downloaded_files}


def main():
    parser = argparse.ArgumentParser(description="Download SCP study")
    parser.add_argument("--accession", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--omic-type", default="unknown")
    args = parser.parse_args()
    result = download_scp(args.accession, args.output_dir, args.omic_type)
    json.dump(result, sys.stdout, indent=2)


if __name__ == "__main__":
    main()

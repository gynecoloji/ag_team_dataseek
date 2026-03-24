#!/usr/bin/env python3
"""Download GEO datasets with metadata and supplementary files."""

import argparse
import json
import os
import sys

from scripts.utils import fetch_with_retry, download_file, resolve_doi, setup_logger

GEO_SOFT_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

logger = setup_logger("download_geo")


def _accession_range(accession):
    numeric = accession.replace("GSE", "")
    return f"GSE{numeric[:-3]}nnn"


def build_geo_ftp_urls(accession):
    acc_range = _accession_range(accession)
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{acc_range}/{accession}/suppl/"
    return [base, f"{base}{accession}_RAW.tar"]


def _fetch_geo_metadata(accession):
    params = {"acc": accession, "targ": "self", "form": "text", "view": "brief"}
    resp = fetch_with_retry(GEO_SOFT_URL, params=params)
    if resp.status_code != 200:
        return {}
    metadata = {"title": "", "summary": "", "overall_design": "", "pubmed_ids": [], "samples": [], "platform": "", "organism": ""}
    for line in resp.text.split("\n"):
        if line.startswith("!Series_title"):
            metadata["title"] = line.split("=", 1)[1].strip() if "=" in line else ""
        elif line.startswith("!Series_summary"):
            metadata["summary"] += line.split("=", 1)[1].strip() + " " if "=" in line else ""
        elif line.startswith("!Series_overall_design"):
            metadata["overall_design"] = line.split("=", 1)[1].strip() if "=" in line else ""
        elif line.startswith("!Series_pubmed_id"):
            val = line.split("=", 1)[1].strip() if "=" in line else ""
            if val:
                metadata["pubmed_ids"].append(val)
        elif line.startswith("!Series_sample_id"):
            val = line.split("=", 1)[1].strip() if "=" in line else ""
            if val:
                metadata["samples"].append(val)
        elif line.startswith("!Series_platform_id"):
            metadata["platform"] = line.split("=", 1)[1].strip() if "=" in line else ""
        elif line.startswith("!Platform_organism"):
            metadata["organism"] = line.split("=", 1)[1].strip() if "=" in line else ""
    return metadata


def _fetch_supplementary_file_list(accession):
    acc_range = _accession_range(accession)
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{acc_range}/{accession}/suppl/"
    try:
        resp = fetch_with_retry(url)
        if resp.status_code != 200:
            return []
        files = []
        for line in resp.text.split("\n"):
            if 'href="' in line and accession in line:
                href_start = line.index('href="') + 6
                href_end = line.index('"', href_start)
                filename = line[href_start:href_end]
                if filename and not filename.endswith("/"):
                    files.append({"name": filename, "url": url + filename, "size_mb": 0})
        return files
    except Exception as e:
        logger.warning(f"Failed to list supplementary files: {e}")
        return []


def generate_summary_report(accession, title, organism, omic_type, platform, disease, tissue, sample_count, condition_groups, data_files, paper, metadata_quality):
    conditions = ", ".join(condition_groups) if condition_groups else "unknown"
    files_list = "\n".join(f"  - {f.get('name', 'unknown')} ({f.get('size_mb', 0)} MB)" for f in data_files)
    return f"""# {accession}

## Dataset Metadata
- **Accession:** {accession}
- **Title:** {title}
- **Organism:** {organism}
- **Tissue:** {tissue or 'not specified'}
- **Disease/Condition:** {disease or 'not specified'}
- **Sample Count:** {sample_count}

## Experimental Design
- **Omic Type:** {omic_type}
- **Platform:** {platform}
- **Conditions/Groups:** {conditions}

## Data Availability
{files_list if files_list.strip() else '  - No supplementary files listed'}

## Original Paper
- **Title:** {paper.get('title', 'N/A')}
- **Authors:** {paper.get('authors', 'N/A')}
- **Journal:** {paper.get('journal', 'N/A')}
- **DOI:** {paper.get('doi', 'N/A')}
- **Date:** {paper.get('date', 'N/A')}
- **Abstract:** {paper.get('abstract', 'N/A')}

## Quality Indicators
- **Metadata Quality:** {metadata_quality}
- **Samples per Group:** {conditions}
"""


def download_geo(accession, output_dir, omic_type="unknown"):
    os.makedirs(output_dir, exist_ok=True)
    data_dir = os.path.join(output_dir, "data")
    meta_dir = os.path.join(output_dir, "metadata")
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(meta_dir, exist_ok=True)

    logger.info(f"Downloading GEO dataset: {accession}")

    metadata = _fetch_geo_metadata(accession)
    with open(os.path.join(meta_dir, "study_metadata.json"), "w") as f:
        json.dump(metadata, f, indent=2)

    paper = {"title": "", "authors": "", "journal": "", "doi": "", "date": "", "abstract": metadata.get("summary", "")}
    if metadata.get("pubmed_ids"):
        pmid = metadata["pubmed_ids"][0]
        try:
            resp = fetch_with_retry("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", params={"db": "pubmed", "id": pmid, "retmode": "json"})
            if resp.status_code == 200:
                pub_data = resp.json().get("result", {}).get(pmid, {})
                doi = ""
                for aid in pub_data.get("articleids", []):
                    if aid.get("idtype") == "doi":
                        doi = aid.get("value", "")
                        break
                if doi:
                    paper = resolve_doi(doi)
                else:
                    paper["title"] = pub_data.get("title", "")
                    paper["authors"] = ", ".join(a.get("name", "") for a in pub_data.get("authors", []))
                    paper["journal"] = pub_data.get("fulljournalname", "")
                    paper["date"] = pub_data.get("pubdate", "")
        except Exception as e:
            logger.warning(f"Failed to resolve paper for PMID {pmid}: {e}")

    with open(os.path.join(meta_dir, "paper.json"), "w") as f:
        json.dump(paper, f, indent=2)

    supp_files = _fetch_supplementary_file_list(accession)
    downloaded_files = []
    for sf in supp_files:
        dest = os.path.join(data_dir, sf["name"])
        try:
            logger.info(f"Downloading {sf['name']}...")
            download_file(sf["url"], dest)
            size_mb = round(os.path.getsize(dest) / 1_000_000, 1)
            downloaded_files.append({"name": sf["name"], "size_mb": size_mb})
        except Exception as e:
            logger.warning(f"Failed to download {sf['name']}: {e}")

    sample_metadata_path = os.path.join(meta_dir, "sample_metadata.csv")
    with open(sample_metadata_path, "w") as f:
        f.write("sample_id,condition,batch\n")
        for sample_id in metadata.get("samples", []):
            f.write(f"{sample_id},TODO,TODO\n")

    return {"accession": accession, "output_dir": output_dir, "files_downloaded": len(downloaded_files), "data_files": downloaded_files, "paper": paper, "metadata": metadata}


def main():
    parser = argparse.ArgumentParser(description="Download GEO dataset")
    parser.add_argument("--accession", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--omic-type", default="unknown")
    args = parser.parse_args()
    result = download_geo(args.accession, args.output_dir, args.omic_type)
    json.dump(result, sys.stdout, indent=2)


if __name__ == "__main__":
    main()

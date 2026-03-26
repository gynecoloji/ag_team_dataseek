#!/usr/bin/env python3
"""Search UCSC Xena for omics datasets via xenaPython."""

import argparse
import json
import sys

try:
    import xenaPython as xena
except ImportError:
    xena = None

from scripts.utils import fetch_with_retry, setup_logger
from scripts.sample_utils import build_sample_table, SAMPLE_CAP
from scripts.supplement_fetch import fetch_supplementary_tables

XENA_HUB_URL = "https://ucscpublic.xenahubs.net"
TCGA_HUB_URL = "https://tcga.xenahubs.net"
TOIL_HUB_URL = "https://toil.xenahubs.net"
ALL_HUBS = [TCGA_HUB_URL, XENA_HUB_URL, TOIL_HUB_URL]

TCGA_SAMPLE_TYPE_MAP = {
    "primary tumor": "primary",
    "primary solid tumor": "primary",
    "recurrent tumor": "primary",
    "recurrent solid tumor": "primary",
    "metastatic": "metastasis",
    "solid tissue normal": "normal",
    "blood derived normal": "normal",
}

OMIC_TO_XENA_SUBTYPES = {
    "bulkRNAseq": ["gene expression rnaseq", "gene expression", "rnaseq"],
    "bulkGenomicseq": ["copy number", "somatic mutation", "mutation", "cnv"],
}

logger = setup_logger("search_xena")


def filter_datasets_by_omic(datasets, omic):
    subtypes = OMIC_TO_XENA_SUBTYPES.get(omic, [])
    if not subtypes:
        return []
    return [ds for ds in datasets if any(st in ds.get("dataSubType", "").lower() for st in subtypes)]


def parse_xena_dataset(ds, omic):
    cohort = ds.get("cohort", "")
    name = ds.get("name", "")
    label = ds.get("label", "")
    return {
        "accession": f"XENA_{name.replace('/', '_').replace('.', '_')}",
        "source": "Xena",
        "title": f"{cohort} - {label}",
        "organism": "Homo sapiens",
        "omic_type": omic,
        "platform": ds.get("dataSubType", "unknown"),
        "disease": cohort,
        "tissue": "",
        "sample_count": ds.get("sampleCount", 0),
        "condition_groups": [],
        "data_files": [{"type": "genomic_matrix", "format": "tsv", "size_mb": 0}],
        "paper": {
            "title": cohort,
            "authors": "",
            "journal": "",
            "doi": "",
            "date": "",
            "abstract": f"Xena dataset: {label}. Cohort: {cohort}.",
        },
        "metadata_quality": "good" if ds.get("sampleCount", 0) > 0 else "partial",
        "date_submitted": "",
    }


def parse_xena_phenotypes(phenotype_data: list[dict], cohort: str = "") -> list[dict]:
    """Parse Xena phenotype data into sample rows.

    Args:
        phenotype_data: List of phenotype dicts from Xena dataset_phenotypes.
        cohort: Cohort name to attach as extra column.

    Returns:
        List of sample row dicts (capped at SAMPLE_CAP).
    """
    rows = []
    for entry in phenotype_data[:SAMPLE_CAP]:
        sample_id = entry.get("sampleID", "")
        disease = entry.get("_primary_disease", "")
        raw_sample_type = entry.get("sample_type", "")
        mapped_sample_type = TCGA_SAMPLE_TYPE_MAP.get(raw_sample_type.lower(), "")
        age = entry.get("age_at_initial_pathologic_diagnosis", "")
        sex = entry.get("gender", entry.get("sex", ""))
        stage = entry.get("pathologic_stage", entry.get("stage", ""))
        primary_site = entry.get("primary_site", "")
        rows.append({
            "sample_id": sample_id,
            "disease": disease,
            "sample_type": mapped_sample_type,
            "age": age,
            "sex": sex,
            "stage": stage,
            "tissue_site": primary_site,
            "cohort_name": cohort,
        })
    return rows


def fetch_xena_samples(hub: str, dataset_name: str, cohort: str = "", paper: dict | None = None) -> dict | None:
    """Fetch sample phenotype data from a Xena dataset.

    Args:
        hub: Xena hub URL.
        dataset_name: Dataset identifier on the hub.
        cohort: Cohort name for annotation.
        paper: Optional paper dict (unused, for API compatibility).

    Returns:
        Sample table dict from build_sample_table, or None on failure.
    """
    if xena is None:
        return None
    try:
        phenotype_data = xena.dataset_phenotypes(hub, dataset_name)
        if not phenotype_data:
            return None
        rows = parse_xena_phenotypes(phenotype_data, cohort=cohort)
        if not rows:
            return None
        return build_sample_table(rows, source="Xena")
    except Exception as e:
        logger.warning(f"fetch_xena_samples failed for {dataset_name} on {hub}: {e}")
        return None


def search_xena(omic, organism="human", disease=None, tissue=None, max_results=50):
    if xena is None:
        logger.warning("xenaPython not installed, falling back to API")
        return _search_xena_api(omic, disease, max_results)
    results = []
    for hub in ALL_HUBS:
        try:
            datasets = xena.all_datasets(hub)
            if not datasets:
                continue
            filtered = filter_datasets_by_omic(datasets, omic)
            for ds in filtered:
                ds["host"] = hub
                if disease and disease.lower() not in ds.get("cohort", "").lower():
                    continue
                result = parse_xena_dataset(ds, omic)
                dataset_name = ds.get("name", "")
                cohort = ds.get("cohort", "")
                sample_table = fetch_xena_samples(hub, dataset_name, cohort=cohort)
                if sample_table is not None:
                    result["sample_table"] = sample_table
                results.append(result)
        except Exception as e:
            logger.warning(f"Failed to query Xena hub {hub}: {e}")
    logger.info(f"Found {len(results)} Xena datasets")
    return results[:max_results]


def _search_xena_api(omic, disease, max_results):
    results = []
    for hub in ALL_HUBS:
        try:
            resp = fetch_with_retry(f"{hub}/data/", params={"format": "json"})
            if resp.status_code != 200:
                continue
            datasets = resp.json()
            filtered = filter_datasets_by_omic(datasets, omic)
            for ds in filtered:
                if disease and disease.lower() not in ds.get("cohort", "").lower():
                    continue
                results.append(parse_xena_dataset(ds, omic))
        except Exception as e:
            logger.warning(f"Xena API fallback failed for {hub}: {e}")
    return results[:max_results]


def main():
    parser = argparse.ArgumentParser(description="Search UCSC Xena for omics datasets")
    parser.add_argument("--omic", required=True)
    parser.add_argument("--organism", default="human")
    parser.add_argument("--disease", default=None)
    parser.add_argument("--tissue", default=None)
    parser.add_argument("--max-results", type=int, default=50)
    args = parser.parse_args()
    results = search_xena(omic=args.omic, organism=args.organism, disease=args.disease, tissue=args.tissue, max_results=args.max_results)
    json.dump(results, sys.stdout, indent=2)


if __name__ == "__main__":
    main()

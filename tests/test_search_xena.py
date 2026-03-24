import json

import pytest

from scripts.search_xena import search_xena, parse_xena_dataset, filter_datasets_by_omic
from tests.conftest import SEARCH_RESULT_REQUIRED_KEYS


SAMPLE_XENA_DATASET = {
    "name": "TCGA-GBM.htseq_fpkm.tsv",
    "label": "Gene Expression RNAseq - HTSeq - FPKM",
    "cohort": "TCGA Glioblastoma (GBM)",
    "type": "genomicMatrix",
    "dataSubType": "gene expression RNAseq",
    "unit": "fpkm",
    "host": "https://tcga.xenahubs.net",
    "sampleCount": 174,
}


class TestParsXenaDataset:
    def test_parses_to_search_result(self):
        result = parse_xena_dataset(SAMPLE_XENA_DATASET, omic="bulkRNAseq")
        assert result["source"] == "Xena"
        assert result["sample_count"] == 174
        assert "TCGA" in result["title"] or "TCGA" in result["accession"]
        for key in SEARCH_RESULT_REQUIRED_KEYS:
            assert key in result


class TestFilterDatasetsByOmic:
    def test_filters_expression_for_bulkRNAseq(self):
        datasets = [
            {"dataSubType": "gene expression RNAseq", "label": "RNA"},
            {"dataSubType": "copy number", "label": "CN"},
            {"dataSubType": "somatic mutation", "label": "Mut"},
        ]
        filtered = filter_datasets_by_omic(datasets, "bulkRNAseq")
        assert len(filtered) == 1
        assert filtered[0]["label"] == "RNA"

    def test_filters_cnv_for_bulkGenomicseq(self):
        datasets = [
            {"dataSubType": "gene expression RNAseq", "label": "RNA"},
            {"dataSubType": "copy number", "label": "CN"},
        ]
        filtered = filter_datasets_by_omic(datasets, "bulkGenomicseq")
        assert len(filtered) == 1
        assert filtered[0]["label"] == "CN"

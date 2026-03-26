import json

import pytest
from unittest.mock import patch, MagicMock

from scripts.search_xena import search_xena, parse_xena_dataset, filter_datasets_by_omic, parse_xena_phenotypes, TCGA_SAMPLE_TYPE_MAP
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

SAMPLE_PHENOTYPE_DATA = [
    {"sampleID": "TCGA-06-0125-01", "_primary_disease": "GBM", "sample_type": "Primary Tumor", "age_at_initial_pathologic_diagnosis": "55", "gender": "MALE", "_PATIENT": "TCGA-06-0125"},
    {"sampleID": "TCGA-06-0126-01", "_primary_disease": "GBM", "sample_type": "Primary Tumor", "age_at_initial_pathologic_diagnosis": "42", "gender": "FEMALE", "_PATIENT": "TCGA-06-0126"},
    {"sampleID": "TCGA-06-0127-06", "_primary_disease": "GBM", "sample_type": "Metastatic", "age_at_initial_pathologic_diagnosis": "60", "gender": "MALE", "_PATIENT": "TCGA-06-0127"},
]


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


class TestParseXenaPhenotypes:
    def test_maps_tcga_phenotype_fields(self):
        rows = parse_xena_phenotypes(SAMPLE_PHENOTYPE_DATA, cohort="TCGA GBM")
        assert len(rows) == 3
        assert rows[0]["sample_id"] == "TCGA-06-0125-01"
        assert rows[0]["disease"] == "GBM"
        assert rows[0]["sample_type"] == "primary"
        assert rows[0]["age"] == "55"
        assert rows[0]["sex"] == "MALE"
        assert rows[0]["cohort_name"] == "TCGA GBM"

    def test_maps_metastatic_sample_type(self):
        rows = parse_xena_phenotypes(SAMPLE_PHENOTYPE_DATA, cohort="TCGA GBM")
        assert rows[2]["sample_type"] == "metastasis"

    def test_caps_at_500(self):
        big_data = [{"sampleID": f"TCGA-{i:04d}"} for i in range(600)]
        rows = parse_xena_phenotypes(big_data)
        assert len(rows) == 500

    def test_empty_input_returns_empty_list(self):
        rows = parse_xena_phenotypes([])
        assert rows == []

    def test_unknown_sample_type_preserved_or_empty(self):
        data = [{"sampleID": "S1", "sample_type": "Unknown Type XYZ"}]
        rows = parse_xena_phenotypes(data)
        # Unknown types not in map should result in empty string or original value
        assert "sample_type" in rows[0]

    def test_cohort_name_defaults_to_empty_string(self):
        rows = parse_xena_phenotypes(SAMPLE_PHENOTYPE_DATA)
        assert rows[0]["cohort_name"] == ""


class TestTcgaSampleTypeMap:
    def test_primary_tumor_maps_to_primary(self):
        assert TCGA_SAMPLE_TYPE_MAP["primary tumor"] == "primary"

    def test_metastatic_maps_to_metastasis(self):
        assert TCGA_SAMPLE_TYPE_MAP["metastatic"] == "metastasis"

    def test_solid_tissue_normal_maps_to_normal(self):
        assert TCGA_SAMPLE_TYPE_MAP["solid tissue normal"] == "normal"

    def test_blood_derived_normal_maps_to_normal(self):
        assert TCGA_SAMPLE_TYPE_MAP["blood derived normal"] == "normal"


class TestFetchXenaSamples:
    def test_returns_none_when_xena_is_none(self):
        """When xenaPython is not available, fetch_xena_samples returns None."""
        from scripts.search_xena import fetch_xena_samples
        with patch("scripts.search_xena.xena", None):
            result = fetch_xena_samples("https://tcga.xenahubs.net", "TCGA-GBM.htseq_fpkm.tsv", cohort="TCGA GBM")
        assert result is None

    def test_returns_sample_table_on_success(self):
        """When xenaPython returns phenotype data, build_sample_table is called."""
        from scripts.search_xena import fetch_xena_samples
        mock_xena = MagicMock()
        mock_xena.dataset_samples.return_value = ["TCGA-06-0125-01", "TCGA-06-0126-01"]
        mock_xena.dataset_phenotypes.return_value = [
            {"sampleID": "TCGA-06-0125-01", "_primary_disease": "GBM", "sample_type": "Primary Tumor",
             "age_at_initial_pathologic_diagnosis": "55", "gender": "MALE"},
            {"sampleID": "TCGA-06-0126-01", "_primary_disease": "GBM", "sample_type": "Primary Tumor",
             "age_at_initial_pathologic_diagnosis": "42", "gender": "FEMALE"},
        ]
        with patch("scripts.search_xena.xena", mock_xena):
            result = fetch_xena_samples("https://tcga.xenahubs.net", "TCGA-GBM.htseq_fpkm.tsv", cohort="TCGA GBM")
        assert result is not None
        assert "rows" in result
        assert "columns" in result
        assert result["shown_samples"] == 2

    def test_returns_none_on_exception(self):
        """If xenaPython raises an exception, return None gracefully."""
        from scripts.search_xena import fetch_xena_samples
        mock_xena = MagicMock()
        mock_xena.dataset_samples.side_effect = Exception("Network error")
        with patch("scripts.search_xena.xena", mock_xena):
            result = fetch_xena_samples("https://tcga.xenahubs.net", "TCGA-GBM.htseq_fpkm.tsv")
        assert result is None

    def test_returns_none_on_empty_phenotypes(self):
        """If dataset_phenotypes returns empty, return None gracefully."""
        from scripts.search_xena import fetch_xena_samples
        mock_xena = MagicMock()
        mock_xena.dataset_samples.return_value = []
        mock_xena.dataset_phenotypes.return_value = []
        with patch("scripts.search_xena.xena", mock_xena):
            result = fetch_xena_samples("https://tcga.xenahubs.net", "TCGA-GBM.htseq_fpkm.tsv")
        assert result is None

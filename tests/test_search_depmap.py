import json

import pytest
import responses

import scripts.search_depmap as search_depmap_module
from scripts.search_depmap import search_depmap, parse_depmap_dataset, fetch_depmap_samples
from tests.conftest import SEARCH_RESULT_REQUIRED_KEYS


SAMPLE_DEPMAP_DATASETS = [
    {"fileName": "CRISPRGeneEffect.csv", "fileDescription": "CRISPR knockout gene effect scores (Chronos)", "releaseName": "DepMap Public 24Q4", "tagsUrl": "crispr", "size": 300000000},
    {"fileName": "CRISPRGeneDependency.csv", "fileDescription": "CRISPR gene dependency probability scores", "releaseName": "DepMap Public 24Q4", "tagsUrl": "crispr", "size": 250000000},
    {"fileName": "OmicsExpressionProteinCodingGenesTPMLogp1.csv", "fileDescription": "Gene expression TPM values", "releaseName": "DepMap Public 24Q4", "tagsUrl": "expression", "size": 450000000},
]


class TestParseDepmapDataset:
    def test_parses_to_search_result(self):
        result = parse_depmap_dataset(SAMPLE_DEPMAP_DATASETS[0], omic="CRISPR")
        assert result["source"] == "DepMap"
        for key in SEARCH_RESULT_REQUIRED_KEYS:
            assert key in result


class TestSearchDepmap:
    @responses.activate
    def test_filters_crispr_datasets(self):
        responses.add(responses.GET, "https://depmap.org/portal/api/download/all", json=SAMPLE_DEPMAP_DATASETS, status=200)
        results = search_depmap(omic="CRISPR")
        assert len(results) >= 1
        assert all(r["source"] == "DepMap" for r in results)

    @responses.activate
    def test_returns_empty_for_unmapped_omic(self):
        responses.add(responses.GET, "https://depmap.org/portal/api/download/all", json=SAMPLE_DEPMAP_DATASETS, status=200)
        results = search_depmap(omic="ChIPseq")
        assert results == []


@pytest.fixture(autouse=False)
def reset_model_metadata_cache():
    """Reset module-level _model_metadata_cache before and after each test."""
    search_depmap_module._model_metadata_cache = None
    yield
    search_depmap_module._model_metadata_cache = None


MOCK_DOWNLOAD_API_WITH_MODEL = [
    {
        "fileName": "Model.csv",
        "fileDescription": "Model metadata",
        "tagsUrl": "model",
        "size": 1000,
        "releaseName": "24Q4",
        "downloadUrl": "https://depmap.org/portal/api/download/Model.csv",
    }
]

MODEL_CSV_BODY = (
    "ModelID,CellLineName,OncotreeLineage,OncotreeSubtype,OncotreePrimaryDisease,COSMICID\n"
    "ACH-000001,A549,Lung,NSCLC,Lung Cancer,905949\n"
    "ACH-000002,MCF7,Breast,Breast Invasive Ductal Carcinoma,Breast Cancer,905946\n"
)


class TestFetchDepmapSamples:
    @responses.activate
    def test_fetches_model_metadata(self, reset_model_metadata_cache):
        responses.add(responses.GET, "https://depmap.org/portal/api/download/all",
            json=[{"fileName": "Model.csv", "fileDescription": "Model metadata", "tagsUrl": "model", "size": 1000, "releaseName": "24Q4", "downloadUrl": "https://depmap.org/portal/api/download/Model.csv"}],
            status=200)
        responses.add(responses.GET, "https://depmap.org/portal/api/download/Model.csv",
            body="ModelID,CellLineName,OncotreeLineage,OncotreeSubtype,OncotreePrimaryDisease,COSMICID\nACH-000001,A549,Lung,NSCLC,Lung Cancer,905949\n",
            status=200)
        table = fetch_depmap_samples()
        assert table is not None
        assert table["total_samples"] >= 1

    @responses.activate
    def test_returns_none_when_model_csv_missing(self, reset_model_metadata_cache):
        responses.add(
            responses.GET,
            "https://depmap.org/portal/api/download/all",
            json=[{"fileName": "OtherFile.csv", "fileDescription": "Other", "tagsUrl": "other", "size": 0, "releaseName": "24Q4", "downloadUrl": "https://depmap.org/portal/api/download/OtherFile.csv"}],
            status=200,
        )
        table = fetch_depmap_samples()
        assert table is None

    @responses.activate
    def test_uses_cache_on_second_call(self, reset_model_metadata_cache):
        responses.add(
            responses.GET,
            "https://depmap.org/portal/api/download/all",
            json=MOCK_DOWNLOAD_API_WITH_MODEL,
            status=200,
        )
        responses.add(
            responses.GET,
            "https://depmap.org/portal/api/download/Model.csv",
            body=MODEL_CSV_BODY,
            status=200,
        )
        fetch_depmap_samples()
        table = fetch_depmap_samples()
        assert table is not None
        assert len(responses.calls) == 2

    @responses.activate
    def test_search_depmap_attaches_sample_table(self, reset_model_metadata_cache):
        responses.add(
            responses.GET,
            "https://depmap.org/portal/api/download/all",
            json=SAMPLE_DEPMAP_DATASETS + MOCK_DOWNLOAD_API_WITH_MODEL,
            status=200,
        )
        responses.add(
            responses.GET,
            "https://depmap.org/portal/api/download/Model.csv",
            body=MODEL_CSV_BODY,
            status=200,
        )
        results = search_depmap(omic="CRISPR")
        assert len(results) >= 1
        for r in results:
            assert "sample_table" in r
            assert r["sample_table"] is not None
            assert "columns" in r["sample_table"]

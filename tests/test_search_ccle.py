import json

import pytest
import responses

import scripts.search_ccle as search_ccle_module
from scripts.search_ccle import search_ccle, parse_ccle_dataset, fetch_ccle_samples
from tests.conftest import SEARCH_RESULT_REQUIRED_KEYS


SAMPLE_DEPMAP_DOWNLOAD_LIST = [
    {
        "fileName": "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
        "fileDescription": "Gene expression TPM values (log2(TPM+1)) for protein coding genes",
        "releaseName": "DepMap Public 24Q4",
        "tagsUrl": "expression",
        "size": 450000000,
    },
    {
        "fileName": "OmicsCNGene.csv",
        "fileDescription": "Gene-level copy number data",
        "releaseName": "DepMap Public 24Q4",
        "tagsUrl": "copy_number",
        "size": 120000000,
    },
    {
        "fileName": "OmicsSomaticMutations.csv",
        "fileDescription": "Somatic mutation calls",
        "releaseName": "DepMap Public 24Q4",
        "tagsUrl": "mutation",
        "size": 80000000,
    },
]


class TestParseCcleDataset:
    def test_parses_download_entry(self):
        entry = SAMPLE_DEPMAP_DOWNLOAD_LIST[0]
        result = parse_ccle_dataset(entry, omic="bulkRNAseq")
        assert result["source"] == "CCLE"
        assert result["organism"] == "Homo sapiens"
        for key in SEARCH_RESULT_REQUIRED_KEYS:
            assert key in result


class TestSearchCcle:
    @responses.activate
    def test_returns_matching_datasets(self):
        responses.add(responses.GET, "https://depmap.org/portal/api/download/all", json=SAMPLE_DEPMAP_DOWNLOAD_LIST, status=200)
        results = search_ccle(omic="bulkRNAseq")
        assert isinstance(results, list)
        assert len(results) >= 1
        assert all(r["source"] == "CCLE" for r in results)

    @responses.activate
    def test_filters_by_omic_type(self):
        responses.add(responses.GET, "https://depmap.org/portal/api/download/all", json=SAMPLE_DEPMAP_DOWNLOAD_LIST, status=200)
        results = search_ccle(omic="bulkRNAseq")
        assert any("expression" in r["title"].lower() or "expression" in r.get("accession", "").lower() for r in results)


@pytest.fixture(autouse=False)
def reset_model_metadata_cache():
    """Reset module-level _model_metadata_cache before and after each test."""
    search_ccle_module._model_metadata_cache = None
    yield
    search_ccle_module._model_metadata_cache = None


MODEL_CSV_BODY = (
    "ModelID,CellLineName,OncotreeLineage,OncotreeSubtype,OncotreePrimaryDisease,COSMICID\n"
    "ACH-000001,A549,Lung,Non-Small Cell Lung Cancer,Lung Cancer,905949\n"
    "ACH-000002,MCF7,Breast,Breast Invasive Ductal Carcinoma,Breast Cancer,905946\n"
)

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


class TestFetchCcleSamples:
    @responses.activate
    def test_fetches_and_maps_model_metadata(self, reset_model_metadata_cache):
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
        table = fetch_ccle_samples()
        assert table is not None
        assert table["total_samples"] >= 1
        cols = table["columns"]
        type_idx = cols.index("sample_type")
        assert table["rows"][0][type_idx] == "cell_line"
        # Verify CCLE extras
        assert "lineage" in cols
        assert "sublineage" in cols
        assert "cosmic_id" in cols

    @responses.activate
    def test_returns_none_when_model_csv_missing(self, reset_model_metadata_cache):
        responses.add(
            responses.GET,
            "https://depmap.org/portal/api/download/all",
            json=[{"fileName": "OtherFile.csv", "fileDescription": "Other", "tagsUrl": "other", "size": 0, "releaseName": "24Q4", "downloadUrl": "https://depmap.org/portal/api/download/OtherFile.csv"}],
            status=200,
        )
        table = fetch_ccle_samples()
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
        # First call populates cache
        fetch_ccle_samples()
        # Second call should use cache — no new HTTP calls needed
        table = fetch_ccle_samples()
        assert table is not None
        # Only 2 HTTP calls total (API list + CSV download), not 4
        assert len(responses.calls) == 2

    @responses.activate
    def test_search_ccle_attaches_sample_table(self, reset_model_metadata_cache):
        responses.add(
            responses.GET,
            "https://depmap.org/portal/api/download/all",
            json=SAMPLE_DEPMAP_DOWNLOAD_LIST + MOCK_DOWNLOAD_API_WITH_MODEL,
            status=200,
        )
        responses.add(
            responses.GET,
            "https://depmap.org/portal/api/download/Model.csv",
            body=MODEL_CSV_BODY,
            status=200,
        )
        results = search_ccle(omic="bulkRNAseq")
        assert len(results) >= 1
        for r in results:
            assert "sample_table" in r
            assert r["sample_table"] is not None
            assert "columns" in r["sample_table"]

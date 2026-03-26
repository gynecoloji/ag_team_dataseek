import json

import pytest
import responses

from scripts.search_scp import search_scp, parse_scp_study, fetch_scp_samples
from tests.conftest import SEARCH_RESULT_REQUIRED_KEYS


SAMPLE_SCP_RESPONSE = {
    "studies": [
        {
            "accession": "SCP1234",
            "name": "Single-cell atlas of glioblastoma",
            "description": "Comprehensive single-cell transcriptomic atlas of glioblastoma tumors",
            "cell_count": 50000,
            "gene_count": 20000,
            "disease": ["glioblastoma"],
            "organ": ["brain"],
            "species": ["Homo sapiens"],
            "library_preparation_protocol": ["10x 3' v3"],
            "study_url": "https://singlecell.broadinstitute.org/single_cell/study/SCP1234",
        }
    ],
    "total_studies": 1,
}


class TestParsScpStudy:
    def test_parses_to_search_result(self):
        study = SAMPLE_SCP_RESPONSE["studies"][0]
        result = parse_scp_study(study, omic="scRNAseq")
        assert result["accession"] == "SCP1234"
        assert result["source"] == "SCP"
        assert result["organism"] == "Homo sapiens"
        for key in SEARCH_RESULT_REQUIRED_KEYS:
            assert key in result


class TestSearchScp:
    @responses.activate
    def test_returns_matching_studies(self):
        responses.add(responses.GET, "https://singlecell.broadinstitute.org/single_cell/api/v1/search", json=SAMPLE_SCP_RESPONSE, status=200)
        results = search_scp(omic="scRNAseq", disease="glioblastoma", organism="human")
        assert isinstance(results, list)
        assert len(results) >= 1
        assert results[0]["source"] == "SCP"

    @responses.activate
    def test_returns_empty_on_no_results(self):
        responses.add(responses.GET, "https://singlecell.broadinstitute.org/single_cell/api/v1/search", json={"studies": [], "total_studies": 0}, status=200)
        results = search_scp(omic="scRNAseq", disease="nonexistent")
        assert results == []


SAMPLE_SCP_STUDY_DETAIL = {
    "accession": "SCP1234",
    "name": "Single-cell atlas of glioblastoma",
    "description": "Comprehensive single-cell atlas",
    "cell_count": 50000,
    "organ": ["brain"],
    "disease": ["glioblastoma"],
    "library_preparation_protocol": ["10x 3' v3"],
    "species": ["Homo sapiens"],
}


class TestFetchScpSamples:
    @responses.activate
    def test_builds_study_level_sample_table(self):
        responses.add(responses.GET, "https://singlecell.broadinstitute.org/single_cell/api/v1/studies/SCP1234", json=SAMPLE_SCP_STUDY_DETAIL, status=200)
        table = fetch_scp_samples("SCP1234")
        assert table is not None
        assert table["total_samples"] == 1
        cols = table["columns"]
        assert "cell_count" in cols
        assert "library_prep" in cols
        site_idx = cols.index("tissue_site")
        assert table["rows"][0][site_idx] == "brain"

    @responses.activate
    def test_returns_none_on_api_failure(self):
        responses.add(responses.GET, "https://singlecell.broadinstitute.org/single_cell/api/v1/studies/SCP9999", status=404)
        table = fetch_scp_samples("SCP9999")
        assert table is None

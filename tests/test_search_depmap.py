import json

import pytest
import responses

from scripts.search_depmap import search_depmap, parse_depmap_dataset
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

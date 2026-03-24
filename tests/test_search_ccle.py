import json

import pytest
import responses

from scripts.search_ccle import search_ccle, parse_ccle_dataset
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

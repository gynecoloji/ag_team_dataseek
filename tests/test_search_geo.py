import json
import xml.etree.ElementTree as ET

import pytest
import responses

from scripts.search_geo import search_geo, parse_geo_record, build_esearch_query
from tests.conftest import SEARCH_RESULT_REQUIRED_KEYS, PAPER_REQUIRED_KEYS


class TestBuildEsearchQuery:
    def test_basic_query(self):
        q = build_esearch_query(omic="scRNAseq", organism="human")
        assert "single cell" in q.lower() or "scRNA" in q
        assert "Homo sapiens" in q or "human" in q

    def test_disease_filter(self):
        q = build_esearch_query(omic="scRNAseq", organism="human", disease="glioblastoma")
        assert "glioblastoma" in q

    def test_tissue_filter(self):
        q = build_esearch_query(omic="bulkRNAseq", organism="mouse", tissue="liver")
        assert "liver" in q


SAMPLE_GEO_ESEARCH_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<eSearchResult>
    <Count>1</Count>
    <RetMax>1</RetMax>
    <IdList>
        <Id>200123456</Id>
    </IdList>
</eSearchResult>"""

SAMPLE_GEO_ESUMMARY_RESPONSE = """<?xml version="1.0" encoding="UTF-8"?>
<eSummaryResult>
    <DocSum>
        <Id>200123456</Id>
        <Item Name="Accession" Type="String">GSE123456</Item>
        <Item Name="title" Type="String">Single-cell RNA-seq of glioblastoma</Item>
        <Item Name="summary" Type="String">We performed scRNAseq analysis on glioblastoma tumors to characterize cellular heterogeneity.</Item>
        <Item Name="GPL" Type="String">GPL24676</Item>
        <Item Name="GSE" Type="String">GSE123456</Item>
        <Item Name="taxon" Type="String">Homo sapiens</Item>
        <Item Name="gdsType" Type="String">Expression profiling by high throughput sequencing</Item>
        <Item Name="n_samples" Type="Integer">12</Item>
        <Item Name="PDAT" Type="String">2025/05/01</Item>
        <Item Name="PubMedIds" Type="List">
            <Item Name="int" Type="Integer">38000001</Item>
        </Item>
        <Item Name="suppFile" Type="String">H5</Item>
    </DocSum>
</eSummaryResult>"""


class TestParseGeoRecord:
    def test_parses_esummary_to_search_result(self):
        root = ET.fromstring(SAMPLE_GEO_ESUMMARY_RESPONSE)
        doc = root.find("DocSum")
        result = parse_geo_record(doc, omic="scRNAseq")
        assert result["accession"] == "GSE123456"
        assert result["source"] == "GEO"
        assert result["organism"] == "Homo sapiens"
        assert result["sample_count"] == 12
        for key in SEARCH_RESULT_REQUIRED_KEYS:
            assert key in result


class TestSearchGeo:
    @responses.activate
    def test_returns_list_of_search_results(self):
        responses.add(responses.GET, "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", body=SAMPLE_GEO_ESEARCH_RESPONSE, status=200)
        responses.add(responses.GET, "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", body=SAMPLE_GEO_ESUMMARY_RESPONSE, status=200)
        results = search_geo(omic="scRNAseq", organism="human", disease="glioblastoma", max_results=10)
        assert isinstance(results, list)
        assert len(results) >= 1
        assert results[0]["accession"] == "GSE123456"

    @responses.activate
    def test_returns_empty_on_no_results(self):
        empty_response = '<?xml version="1.0"?><eSearchResult><Count>0</Count><RetMax>0</RetMax><IdList></IdList></eSearchResult>'
        responses.add(responses.GET, "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", body=empty_response, status=200)
        results = search_geo(omic="scRNAseq", organism="human", disease="nonexistent_disease_xyz")
        assert results == []

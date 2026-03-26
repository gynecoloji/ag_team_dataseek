import json
import xml.etree.ElementTree as ET

import pytest
import responses

from scripts.search_geo import (
    search_geo,
    parse_geo_record,
    build_esearch_query,
    parse_soft_samples,
    fetch_geo_samples,
)
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


SAMPLE_SOFT_RESPONSE = """\
^SAMPLE = GSM000001
!Sample_characteristics_ch1 = tissue: brain
!Sample_characteristics_ch1 = disease state: glioblastoma
!Sample_characteristics_ch1 = tumor type: primary
!Sample_characteristics_ch1 = treatment: untreated
!Sample_characteristics_ch1 = age: 55
!Sample_characteristics_ch1 = Sex: Male
!Sample_characteristics_ch1 = stage: IV
!Sample_platform_id = GPL24676

^SAMPLE = GSM000002
!Sample_characteristics_ch1 = tissue: spine
!Sample_characteristics_ch1 = disease state: glioblastoma
!Sample_characteristics_ch1 = tumor type: metastasis
!Sample_characteristics_ch1 = treatment: cisplatin
!Sample_characteristics_ch1 = age: 42
!Sample_characteristics_ch1 = Sex: Female
!Sample_characteristics_ch1 = stage: IV
!Sample_platform_id = GPL24676

^SAMPLE = GSM000003
!Sample_characteristics_ch1 = tissue: brain
!Sample_characteristics_ch1 = disease state: normal
!Sample_characteristics_ch1 = cell type: astrocyte
!Sample_characteristics_ch1 = treatment: untreated
!Sample_characteristics_ch1 = age: 60
!Sample_characteristics_ch1 = Sex: Male
!Sample_platform_id = GPL24676
"""


class TestParseSoftSamples:
    def test_parses_sample_characteristics(self):
        rows = parse_soft_samples(SAMPLE_SOFT_RESPONSE)
        gsm1 = next(r for r in rows if r["sample_id"] == "GSM000001")
        assert gsm1["tissue_site"] == "brain"
        assert gsm1["disease"] == "glioblastoma"
        assert gsm1["sample_type"] == "primary"
        assert gsm1["treatment_status"] == "untreated"
        assert gsm1["age"] == "55"
        assert gsm1["sex"] == "Male"
        assert gsm1["stage"] == "IV"
        assert gsm1["platform_id"] == "GPL24676"

    def test_maps_metastasis_sample(self):
        rows = parse_soft_samples(SAMPLE_SOFT_RESPONSE)
        gsm2 = next(r for r in rows if r["sample_id"] == "GSM000002")
        assert gsm2["tissue_site"] == "spine"
        assert gsm2["disease"] == "glioblastoma"
        assert gsm2["sample_type"] == "metastasis"
        assert gsm2["treatment_status"] == "cisplatin"
        assert gsm2["age"] == "42"
        assert gsm2["sex"] == "Female"
        assert gsm2["stage"] == "IV"
        assert gsm2["platform_id"] == "GPL24676"

    def test_maps_cell_type(self):
        rows = parse_soft_samples(SAMPLE_SOFT_RESPONSE)
        gsm3 = next(r for r in rows if r["sample_id"] == "GSM000003")
        assert gsm3["cell_type"] == "astrocyte"
        assert gsm3["tissue_site"] == "brain"

    def test_caps_at_500_samples(self):
        lines = []
        for i in range(600):
            lines.append(f"^SAMPLE = GSM{i:06d}")
            lines.append(f"!Sample_characteristics_ch1 = tissue: brain")
            lines.append(f"!Sample_platform_id = GPL24676")
            lines.append("")
        soft_text = "\n".join(lines)
        rows = parse_soft_samples(soft_text)
        assert len(rows) == 500


class TestFetchGeoSamples:
    @responses.activate
    def test_returns_sample_table_dict(self):
        responses.add(
            responses.GET,
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi",
            body=SAMPLE_SOFT_RESPONSE,
            status=200,
        )
        result = fetch_geo_samples("GSE000001")
        assert result is not None
        assert "rows" in result
        assert "columns" in result
        assert "total_samples" in result
        assert result["total_samples"] == 3

    @responses.activate
    def test_returns_none_on_http_error(self):
        responses.add(
            responses.GET,
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi",
            body=Exception("connection error"),
        )
        result = fetch_geo_samples("GSE_BAD")
        assert result is None


class TestSearchGeo:
    @responses.activate
    def test_returns_list_of_search_results(self):
        responses.add(responses.GET, "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", body=SAMPLE_GEO_ESEARCH_RESPONSE, status=200)
        responses.add(responses.GET, "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", body=SAMPLE_GEO_ESUMMARY_RESPONSE, status=200)
        responses.add(responses.GET, "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi", body=SAMPLE_SOFT_RESPONSE, status=200)
        results = search_geo(omic="scRNAseq", organism="human", disease="glioblastoma", max_results=10)
        assert isinstance(results, list)
        assert len(results) >= 1
        assert results[0]["accession"] == "GSE123456"
        assert "sample_table" in results[0]

    @responses.activate
    def test_returns_empty_on_no_results(self):
        empty_response = '<?xml version="1.0"?><eSearchResult><Count>0</Count><RetMax>0</RetMax><IdList></IdList></eSearchResult>'
        responses.add(responses.GET, "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", body=empty_response, status=200)
        results = search_geo(omic="scRNAseq", organism="human", disease="nonexistent_disease_xyz")
        assert results == []

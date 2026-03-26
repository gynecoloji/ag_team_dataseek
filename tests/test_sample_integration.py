"""Integration tests for sample metadata table generation during search."""

import csv
import json
import os

import pytest
import responses

from scripts.search_geo import search_geo
from scripts.sample_utils import write_sample_csv, summarize_samples
from tests.conftest import SAMPLE_TABLE_KEYS


# ---------------------------------------------------------------------------
# Test data
# ---------------------------------------------------------------------------

SAMPLE_GEO_ESEARCH = """<?xml version="1.0" encoding="UTF-8"?>
<eSearchResult>
    <Count>1</Count>
    <RetMax>1</RetMax>
    <IdList>
        <Id>200123456</Id>
    </IdList>
</eSearchResult>"""

SAMPLE_GEO_ESUMMARY = """<?xml version="1.0" encoding="UTF-8"?>
<eSummaryResult>
    <DocSum>
        <Id>200123456</Id>
        <Item Name="Accession" Type="String">GSE123456</Item>
        <Item Name="title" Type="String">Single-cell RNA-seq of glioblastoma</Item>
        <Item Name="summary" Type="String">We performed scRNAseq on glioblastoma tumors.</Item>
        <Item Name="GPL" Type="String">GPL24676</Item>
        <Item Name="GSE" Type="String">GSE123456</Item>
        <Item Name="taxon" Type="String">Homo sapiens</Item>
        <Item Name="gdsType" Type="String">Expression profiling by high throughput sequencing</Item>
        <Item Name="n_samples" Type="Integer">2</Item>
        <Item Name="PDAT" Type="String">2025/05/01</Item>
        <Item Name="PubMedIds" Type="List">
        </Item>
        <Item Name="suppFile" Type="String">H5</Item>
    </DocSum>
</eSummaryResult>"""

SAMPLE_SOFT = """\
^SAMPLE = GSM000001
!Sample_characteristics_ch1 = tissue: brain
!Sample_characteristics_ch1 = tumor type: primary
!Sample_characteristics_ch1 = treatment: untreated
!Sample_platform_id = GPL24676

^SAMPLE = GSM000002
!Sample_characteristics_ch1 = tissue: brain
!Sample_characteristics_ch1 = tumor type: metastasis
!Sample_characteristics_ch1 = treatment: cisplatin
!Sample_platform_id = GPL24676
"""


def _register_mocks():
    """Register all required HTTP mocks for a complete search_geo call."""
    responses.add(
        responses.GET,
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        body=SAMPLE_GEO_ESEARCH,
        status=200,
    )
    responses.add(
        responses.GET,
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
        body=SAMPLE_GEO_ESUMMARY,
        status=200,
    )
    responses.add(
        responses.GET,
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi",
        body=SAMPLE_SOFT,
        status=200,
    )
    # Mock NCBI ID converter — no PubMedIds in test data, prevents real HTTP calls
    responses.add(
        responses.GET,
        "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/",
        json={"records": []},
        status=200,
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestGeoSearchWithSampleTable:

    @responses.activate
    def test_search_result_includes_sample_table(self):
        _register_mocks()
        results = search_geo(omic="scRNAseq", organism="human", disease="glioblastoma", max_results=10)

        assert len(results) == 1
        result = results[0]

        assert "sample_table" in result
        sample_table = result["sample_table"]
        assert sample_table is not None

        # All required keys must be present
        for key in SAMPLE_TABLE_KEYS:
            assert key in sample_table, f"Missing key: {key}"

        assert sample_table["total_samples"] == 2
        assert len(sample_table["rows"]) == 2

    @responses.activate
    def test_sample_csv_write_from_search_result(self, tmp_path):
        _register_mocks()
        results = search_geo(omic="scRNAseq", organism="human", disease="glioblastoma", max_results=10)

        sample_table = results[0]["sample_table"]
        output_dir = str(tmp_path / "metadata")
        csv_path = write_sample_csv("GSE123456", sample_table, output_dir)

        assert os.path.exists(csv_path)

        with open(csv_path, newline="") as f:
            reader = csv.reader(f)
            rows = list(reader)

        # First row is header
        header = rows[0]
        assert "sample_id" in header
        assert "tissue_site" in header

        # Exactly 2 data rows
        data_rows = [r for r in rows[1:] if r and not r[0].startswith("#")]
        assert len(data_rows) == 2

    @responses.activate
    def test_sample_summary_from_search_result(self):
        _register_mocks()
        results = search_geo(omic="scRNAseq", organism="human", disease="glioblastoma", max_results=10)

        sample_table = results[0]["sample_table"]
        summary = summarize_samples(sample_table)

        assert "Samples (2)" in summary

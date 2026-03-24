import json
import os

import pytest
import responses

from scripts.download_geo import download_geo, build_geo_ftp_urls, generate_summary_report


class TestBuildGeoFtpUrls:
    def test_builds_correct_urls(self):
        urls = build_geo_ftp_urls("GSE123456")
        assert any("GSE123nnn" in url for url in urls)
        assert any("GSE123456" in url for url in urls)


class TestGenerateSummaryReport:
    def test_generates_markdown(self, tmp_path):
        report = generate_summary_report(
            accession="GSE123456",
            title="Test Dataset",
            organism="Homo sapiens",
            omic_type="scRNAseq",
            platform="10x Chromium",
            disease="glioblastoma",
            tissue="brain",
            sample_count=12,
            condition_groups=["tumor", "normal"],
            data_files=[{"name": "matrix.h5", "size_mb": 450}],
            paper={"title": "Paper", "authors": "Smith", "journal": "Nature", "doi": "10.1038/x", "date": "2025-01-01", "abstract": "Abstract"},
            metadata_quality="good",
        )
        assert "# GSE123456" in report
        assert "glioblastoma" in report
        assert "10x Chromium" in report

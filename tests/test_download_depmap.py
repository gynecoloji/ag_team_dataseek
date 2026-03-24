import json
import os
import pytest
import responses
from scripts.download_depmap import download_depmap

class TestDownloadDepmap:
    @responses.activate
    def test_downloads_to_output_dir(self, tmp_path):
        output_dir = str(tmp_path / "DepMap_test")
        responses.add(responses.GET, "https://depmap.org/portal/api/download/all", json=[{"fileName": "CRISPRGeneEffect.csv", "downloadUrl": "https://depmap.org/portal/api/download/file/CRISPRGeneEffect.csv", "size": 100}], status=200)
        responses.add(responses.GET, "https://depmap.org/portal/api/download/file/CRISPRGeneEffect.csv", body=b"gene,cell1\nTP53,-0.5", status=200)
        result = download_depmap("DepMap_CRISPRGeneEffect", output_dir)
        assert result["accession"] == "DepMap_CRISPRGeneEffect"

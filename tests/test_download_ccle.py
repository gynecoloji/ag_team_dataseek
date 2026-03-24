import json
import os
import pytest
import responses
from scripts.download_ccle import download_ccle

class TestDownloadCcle:
    @responses.activate
    def test_downloads_to_output_dir(self, tmp_path):
        output_dir = str(tmp_path / "CCLE_test")
        responses.add(responses.GET, "https://depmap.org/portal/api/download/all", json=[{"fileName": "OmicsExpression.csv", "downloadUrl": "https://depmap.org/portal/api/download/file/OmicsExpression.csv", "size": 100}], status=200)
        responses.add(responses.GET, "https://depmap.org/portal/api/download/file/OmicsExpression.csv", body=b"gene,sample1\nTP53,5.2", status=200)
        result = download_ccle("CCLE_OmicsExpression", output_dir)
        assert result["accession"] == "CCLE_OmicsExpression"
        assert os.path.isdir(output_dir)

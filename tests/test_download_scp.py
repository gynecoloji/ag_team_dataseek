import json
import os
import pytest
import responses
from scripts.download_scp import download_scp

class TestDownloadScp:
    @responses.activate
    def test_downloads_study_files(self, tmp_path):
        output_dir = str(tmp_path / "SCP_test")
        study_info = {"accession": "SCP1234", "name": "Test Study", "description": "A test study", "study_files": [{"name": "expression_matrix.tsv.gz", "download_url": "https://singlecell.broadinstitute.org/data/SCP1234/expression_matrix.tsv.gz"}]}
        responses.add(responses.GET, "https://singlecell.broadinstitute.org/single_cell/api/v1/studies/SCP1234", json=study_info, status=200)
        responses.add(responses.GET, "https://singlecell.broadinstitute.org/data/SCP1234/expression_matrix.tsv.gz", body=b"fake compressed data", status=200)
        result = download_scp("SCP1234", output_dir)
        assert result["accession"] == "SCP1234"

import json
import os
import pytest
from scripts.download_xena import build_xena_download_url

class TestDownloadXena:
    def test_builds_download_url(self):
        url = build_xena_download_url("https://tcga.xenahubs.net", "TCGA-GBM.htseq_fpkm.tsv")
        assert "tcga.xenahubs.net" in url
        assert "TCGA-GBM" in url

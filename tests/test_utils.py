# tests/test_utils.py
import hashlib
import json
import os
import time

import pytest
import responses

from scripts.utils import (
    build_cache_filename,
    fetch_with_retry,
    resolve_doi,
    save_search_cache,
    load_search_cache,
    download_file,
    setup_logger,
)


class TestFetchWithRetry:
    @responses.activate
    def test_successful_request(self):
        responses.add(responses.GET, "https://example.com/api", json={"ok": True}, status=200)
        result = fetch_with_retry("https://example.com/api")
        assert result.status_code == 200
        assert result.json() == {"ok": True}

    @responses.activate
    def test_retries_on_500(self):
        responses.add(responses.GET, "https://example.com/api", status=500)
        responses.add(responses.GET, "https://example.com/api", status=500)
        responses.add(responses.GET, "https://example.com/api", json={"ok": True}, status=200)
        result = fetch_with_retry("https://example.com/api", max_retries=3, base_delay=0.01)
        assert result.status_code == 200

    @responses.activate
    def test_raises_after_max_retries(self):
        responses.add(responses.GET, "https://example.com/api", status=500)
        responses.add(responses.GET, "https://example.com/api", status=500)
        responses.add(responses.GET, "https://example.com/api", status=500)
        with pytest.raises(Exception, match="failed after 3 retries"):
            fetch_with_retry("https://example.com/api", max_retries=3, base_delay=0.01)


class TestResolveDoi:
    @responses.activate
    def test_resolves_doi_to_paper_metadata(self):
        crossref_response = {
            "message": {
                "title": ["Test Paper Title"],
                "author": [{"given": "John", "family": "Doe"}],
                "container-title": ["Nature"],
                "DOI": "10.1038/test",
                "published-print": {"date-parts": [[2025, 6, 15]]},
                "abstract": "Test abstract."
            }
        }
        responses.add(
            responses.GET,
            "https://api.crossref.org/works/10.1038/test",
            json=crossref_response,
            status=200,
        )
        paper = resolve_doi("10.1038/test")
        assert paper["title"] == "Test Paper Title"
        assert paper["authors"] == "John Doe"
        assert paper["journal"] == "Nature"
        assert paper["doi"] == "10.1038/test"
        assert paper["date"] == "2025-06-15"
        assert paper["abstract"] == "Test abstract."

    @responses.activate
    def test_returns_empty_on_failure(self):
        responses.add(
            responses.GET,
            "https://api.crossref.org/works/10.1038/missing",
            status=404,
        )
        paper = resolve_doi("10.1038/missing")
        assert paper["title"] == ""
        assert paper["authors"] == ""


class TestCacheOperations:
    def test_build_cache_filename(self):
        fname = build_cache_filename("scRNAseq", {"disease": "glioblastoma", "organism": "human"})
        # Format: {date}_{omic}_{hash}.json
        assert "scRNAseq_" in fname
        assert fname.endswith(".json")
        # Same params produce same filename
        fname2 = build_cache_filename("scRNAseq", {"disease": "glioblastoma", "organism": "human"})
        assert fname == fname2
        # Different params produce different filename
        fname3 = build_cache_filename("scRNAseq", {"disease": "melanoma", "organism": "human"})
        assert fname != fname3

    def test_save_and_load_search_cache(self, tmp_cache_dir):
        results = [{"accession": "GSE123", "source": "GEO", "title": "Test"}]
        params = {"disease": "glioblastoma", "organism": "human"}
        save_search_cache(results, "scRNAseq", params, tmp_cache_dir)
        loaded = load_search_cache("scRNAseq", params, tmp_cache_dir)
        assert loaded == results

    def test_load_missing_cache_returns_none(self, tmp_cache_dir):
        result = load_search_cache("scRNAseq", {"disease": "nothing"}, tmp_cache_dir)
        assert result is None


class TestDownloadFile:
    @responses.activate
    def test_downloads_file_to_disk(self, tmp_path):
        responses.add(
            responses.GET,
            "https://example.com/data/matrix.h5",
            body=b"fake file content",
            status=200,
            headers={"Content-Length": "17"},
        )
        output_path = str(tmp_path / "matrix.h5")
        download_file("https://example.com/data/matrix.h5", output_path)
        assert os.path.exists(output_path)
        with open(output_path, "rb") as f:
            assert f.read() == b"fake file content"


class TestLogger:
    def test_setup_logger(self, tmp_path):
        log_file = str(tmp_path / "test.log")
        logger = setup_logger("test_logger", log_file)
        logger.info("hello")
        with open(log_file) as f:
            content = f.read()
        assert "hello" in content

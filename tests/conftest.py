import json
import os
import tempfile

import pytest


@pytest.fixture
def tmp_output_dir(tmp_path):
    """Temporary directory for download outputs."""
    return str(tmp_path / "output")


@pytest.fixture
def tmp_cache_dir(tmp_path):
    """Temporary directory for search cache."""
    cache_dir = tmp_path / "search_cache"
    cache_dir.mkdir()
    return str(cache_dir)


@pytest.fixture
def sample_search_result():
    """A single SearchResult dict matching the spec schema."""
    return {
        "accession": "GSE123456",
        "source": "GEO",
        "title": "Single-cell RNA-seq of glioblastoma tumors",
        "organism": "Homo sapiens",
        "omic_type": "scRNAseq",
        "platform": "10x Chromium 3' v3",
        "disease": "glioblastoma",
        "tissue": "brain",
        "sample_count": 12,
        "condition_groups": ["tumor", "normal"],
        "data_files": [
            {"type": "count_matrix", "format": "h5", "size_mb": 450}
        ],
        "paper": {
            "title": "Glioblastoma single-cell landscape",
            "authors": "Smith J, Doe A",
            "journal": "Nature",
            "doi": "10.1038/example",
            "date": "2025-06-15",
            "abstract": "We performed scRNAseq on glioblastoma tumors."
        },
        "metadata_quality": "good",
        "date_submitted": "2025-05-01"
    }


SEARCH_RESULT_REQUIRED_KEYS = {
    "accession", "source", "title", "organism", "omic_type",
    "platform", "disease", "tissue", "sample_count",
    "condition_groups", "data_files", "paper",
    "metadata_quality", "date_submitted"
}

PAPER_REQUIRED_KEYS = {
    "title", "authors", "journal", "doi", "date", "abstract"
}

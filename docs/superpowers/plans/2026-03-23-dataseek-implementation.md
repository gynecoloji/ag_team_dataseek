# Dataseek Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a multi-agent system that searches 5 public omics repositories (GEO, CCLE, Xena UCSC, DepMap, Single Cell Portal), presents tiered results, and downloads selected datasets with metadata and pipeline-ready starter configs.

**Architecture:** Hybrid skill + Python script approach. Two slash commands (`/dataseek-search`, `/dataseek-download`) orchestrate 7 agents (1 coordinator, 5 source searchers, 1 downloader). Skills dispatch subagents in parallel; Python scripts handle API queries, response parsing, and file downloads. Follows the same patterns as the sibling literature_search and scRNAseq projects.

**Tech Stack:** Python 3.10+, biopython, requests, pandas, xenaPython. Conda environment `dataseek`. Claude Code skills + CLAUDE.md agent orchestration.

**Spec:** `docs/superpowers/specs/2026-03-23-dataseek-design.md`

---

## File Structure

```
dataseek/
├── CLAUDE.md                              # Agent orchestration spec (authoritative)
├── environment.yml                        # Conda environment definition
├── .claude/
│   └── skills/
│       ├── dataseek-search/
│       │   └── SKILL.md                   # Search skill definition
│       └── dataseek-download/
│           └── SKILL.md                   # Download skill definition
├── scripts/
│   ├── utils.py                           # Shared: HTTP client, retry, logging, DOI resolver, cache
│   ├── search_geo.py                      # GEO Entrez E-utilities search
│   ├── search_ccle.py                     # CCLE DepMap Portal search
│   ├── search_xena.py                     # UCSC Xena search
│   ├── search_depmap.py                   # DepMap search
│   ├── search_scp.py                      # Single Cell Portal search
│   ├── download_geo.py                    # GEO dataset downloader
│   ├── download_ccle.py                   # CCLE dataset downloader
│   ├── download_xena.py                   # Xena dataset downloader
│   ├── download_depmap.py                 # DepMap dataset downloader
│   └── download_scp.py                    # SCP dataset downloader
├── tests/
│   ├── conftest.py                        # Shared fixtures (mock responses, tmp dirs)
│   ├── test_utils.py                      # Tests for shared utilities
│   ├── test_search_geo.py                 # Tests for GEO search
│   ├── test_search_ccle.py                # Tests for CCLE search
│   ├── test_search_xena.py                # Tests for Xena search
│   ├── test_search_depmap.py              # Tests for DepMap search
│   ├── test_search_scp.py                 # Tests for SCP search
│   ├── test_download_geo.py               # Tests for GEO download
│   ├── test_download_ccle.py              # Tests for CCLE download
│   ├── test_download_xena.py              # Tests for Xena download
│   ├── test_download_depmap.py            # Tests for DepMap download
│   └── test_download_scp.py              # Tests for SCP download
├── results/
│   ├── search_cache/                      # Cached search results (JSON)
│   └── reports/                           # Per-dataset summary reports
└── downloads/                             # Downloaded datasets
```

---

## Task 1: Project Scaffolding

**Files:**
- Create: `environment.yml`
- Create: `scripts/__init__.py` (empty)
- Create: `tests/__init__.py` (empty)
- Create: `tests/conftest.py`
- Create: directory structure for `results/search_cache/`, `results/reports/`, `downloads/`

- [ ] **Step 1: Create environment.yml**

```yaml
name: dataseek
channels:
  - conda-forge
  - defaults
dependencies:
  - python>=3.10
  - biopython
  - requests
  - pandas
  - pip
  - pip:
    - xenaPython
  - pytest
  - pytest-mock
  - responses
```

- [ ] **Step 2: Create directory structure**

```bash
mkdir -p scripts tests results/search_cache results/reports downloads
touch scripts/__init__.py tests/__init__.py
```

- [ ] **Step 3: Create tests/conftest.py with shared fixtures**

```python
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
```

- [ ] **Step 4: Create the conda environment**

```bash
conda env create -f environment.yml
```

- [ ] **Step 5: Verify environment works**

```bash
conda run -n dataseek python -c "import Bio; import requests; import pandas; print('OK')"
```

- [ ] **Step 6: Commit**

```bash
git add environment.yml scripts/__init__.py tests/__init__.py tests/conftest.py
git commit -m "feat: project scaffolding with conda env and test fixtures"
```

---

## Task 2: Shared Utilities (`scripts/utils.py`)

**Files:**
- Create: `scripts/utils.py`
- Create: `tests/test_utils.py`

- [ ] **Step 1: Write failing tests for utils**

```python
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
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
conda run -n dataseek pytest tests/test_utils.py -v
```

Expected: FAIL with import errors (scripts.utils does not exist yet)

- [ ] **Step 3: Implement utils.py**

```python
# scripts/utils.py
"""Shared utilities for dataseek scripts: HTTP, retry, caching, DOI resolution, logging."""

import hashlib
import json
import logging
import os
import time
from datetime import date

import requests


def setup_logger(name: str, log_file: str | None = None, level: int = logging.INFO) -> logging.Logger:
    """Create a logger with timestamp formatting, optionally writing to a file."""
    logger = logging.getLogger(name)
    logger.setLevel(level)
    formatter = logging.Formatter("%(asctime)s [%(name)s] %(levelname)s: %(message)s")

    if not logger.handlers:
        console = logging.StreamHandler()
        console.setFormatter(formatter)
        logger.addHandler(console)

        if log_file:
            os.makedirs(os.path.dirname(log_file), exist_ok=True)
            fh = logging.FileHandler(log_file)
            fh.setFormatter(formatter)
            logger.addHandler(fh)

    return logger


def fetch_with_retry(
    url: str,
    params: dict | None = None,
    headers: dict | None = None,
    max_retries: int = 3,
    base_delay: float = 1.0,
    timeout: int = 30,
) -> requests.Response:
    """GET request with exponential backoff retry on 5xx or connection errors."""
    for attempt in range(max_retries):
        try:
            resp = requests.get(url, params=params, headers=headers, timeout=timeout)
            if resp.status_code < 500:
                return resp
        except requests.RequestException:
            if attempt == max_retries - 1:
                raise
        if attempt < max_retries - 1:
            time.sleep(base_delay * (2 ** attempt))
    raise Exception(f"Request to {url} failed after {max_retries} retries (last status: {resp.status_code})")


def resolve_doi(doi: str) -> dict:
    """Resolve a DOI to paper metadata via CrossRef API. Returns empty fields on failure."""
    empty = {"title": "", "authors": "", "journal": "", "doi": doi, "date": "", "abstract": ""}
    try:
        resp = fetch_with_retry(
            f"https://api.crossref.org/works/{doi}",
            headers={"Accept": "application/json"},
            max_retries=2,
            base_delay=0.5,
        )
        if resp.status_code != 200:
            return empty
        msg = resp.json().get("message", {})
        titles = msg.get("title", [""])
        authors_list = msg.get("author", [])
        authors_str = ", ".join(
            f"{a.get('given', '')} {a.get('family', '')}".strip() for a in authors_list
        )
        journal = msg.get("container-title", [""])[0] if msg.get("container-title") else ""
        date_parts = msg.get("published-print", msg.get("published-online", {})).get("date-parts", [[]])
        if date_parts and date_parts[0]:
            parts = date_parts[0]
            date_str = "-".join(str(p).zfill(2) for p in parts)
        else:
            date_str = ""
        abstract = msg.get("abstract", "")
        return {
            "title": titles[0] if titles else "",
            "authors": authors_str,
            "journal": journal,
            "doi": msg.get("DOI", doi),
            "date": date_str,
            "abstract": abstract,
        }
    except Exception:
        return empty


def build_cache_filename(omic: str, params: dict) -> str:
    """Build a deterministic cache filename from omic type and query params."""
    sorted_params = json.dumps(params, sort_keys=True)
    param_hash = hashlib.sha256(sorted_params.encode()).hexdigest()[:12]
    today = str(date.today())
    return f"{today}_{omic}_{param_hash}.json"


def save_search_cache(results: list[dict], omic: str, params: dict, cache_dir: str) -> str:
    """Save search results to cache directory. Returns the cache file path."""
    os.makedirs(cache_dir, exist_ok=True)
    filename = build_cache_filename(omic, params)
    filepath = os.path.join(cache_dir, filename)
    with open(filepath, "w") as f:
        json.dump({"omic": omic, "params": params, "date": str(date.today()), "results": results}, f, indent=2)
    return filepath


def load_search_cache(omic: str, params: dict, cache_dir: str) -> list[dict] | None:
    """Load cached search results. Returns None if cache miss."""
    filename = build_cache_filename(omic, params)
    filepath = os.path.join(cache_dir, filename)
    if not os.path.exists(filepath):
        return None
    with open(filepath) as f:
        data = json.load(f)
    return data.get("results")


def download_file(url: str, output_path: str, chunk_size: int = 8192) -> None:
    """Download a file with streaming and resume support."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    headers = {}
    existing_size = 0
    if os.path.exists(output_path):
        existing_size = os.path.getsize(output_path)
        headers["Range"] = f"bytes={existing_size}-"

    resp = requests.get(url, headers=headers, stream=True, timeout=60)
    resp.raise_for_status()

    mode = "ab" if existing_size and resp.status_code == 206 else "wb"
    with open(output_path, mode) as f:
        for chunk in resp.iter_content(chunk_size=chunk_size):
            f.write(chunk)
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
conda run -n dataseek pytest tests/test_utils.py -v
```

Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add scripts/utils.py tests/test_utils.py
git commit -m "feat: shared utilities with HTTP retry, DOI resolver, caching, file download"
```

---

## Task 3: GEO Search Script

**Files:**
- Create: `scripts/search_geo.py`
- Create: `tests/test_search_geo.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_search_geo.py
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


# Sample GEO XML response for testing
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
        # Mock esearch
        responses.add(
            responses.GET,
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            body=SAMPLE_GEO_ESEARCH_RESPONSE,
            status=200,
        )
        # Mock esummary
        responses.add(
            responses.GET,
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
            body=SAMPLE_GEO_ESUMMARY_RESPONSE,
            status=200,
        )
        results = search_geo(omic="scRNAseq", organism="human", disease="glioblastoma", max_results=10)
        assert isinstance(results, list)
        assert len(results) >= 1
        assert results[0]["accession"] == "GSE123456"

    @responses.activate
    def test_returns_empty_on_no_results(self):
        empty_response = """<?xml version="1.0"?>
        <eSearchResult><Count>0</Count><RetMax>0</RetMax><IdList></IdList></eSearchResult>"""
        responses.add(
            responses.GET,
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            body=empty_response,
            status=200,
        )
        results = search_geo(omic="scRNAseq", organism="human", disease="nonexistent_disease_xyz")
        assert results == []
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
conda run -n dataseek pytest tests/test_search_geo.py -v
```

Expected: FAIL with import errors

- [ ] **Step 3: Implement search_geo.py**

```python
#!/usr/bin/env python3
"""Search NCBI GEO for omics datasets via Entrez E-utilities."""

import argparse
import json
import os
import sys
import xml.etree.ElementTree as ET

from scripts.utils import fetch_with_retry, resolve_doi, setup_logger

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

OMIC_TO_GEO_TERMS = {
    "scRNAseq": '("single cell" OR "scRNA" OR "10x" OR "Smart-seq") AND "Expression profiling by high throughput sequencing"[Filter]',
    "bulkRNAseq": '"Expression profiling by high throughput sequencing"[Filter] NOT "single cell" NOT "scRNA"',
    "scMultiomicseq": '("single cell" OR "multiome") AND ("ATAC" OR "multi-omic")',
    "scGenomicseq": '("single cell" OR "single-cell") AND ("DNA sequencing" OR "whole genome")',
    "bulkGenomicseq": '("whole genome sequencing" OR "WGS" OR "exome") NOT "single cell"',
    "ATACseq": '("ATAC-seq" OR "ATACseq" OR "chromatin accessibility")',
    "ChIPseq": '("ChIP-seq" OR "ChIPseq" OR "chromatin immunoprecipitation sequencing")',
    "CRISPR": '("CRISPR" OR "genome-wide screen" OR "loss of function screen")',
}

logger = setup_logger("search_geo")


def build_esearch_query(
    omic: str,
    organism: str = "human",
    disease: str | None = None,
    tissue: str | None = None,
    date_from: str | None = None,
    date_to: str | None = None,
) -> str:
    """Build an Entrez search query string for GEO DataSets (gds) database."""
    parts = []
    omic_term = OMIC_TO_GEO_TERMS.get(omic, f'"{omic}"')
    parts.append(omic_term)

    organism_map = {"human": "Homo sapiens", "mouse": "Mus musculus"}
    org = organism_map.get(organism.lower(), organism)
    parts.append(f'"{org}"[Organism]')

    if disease:
        parts.append(f'"{disease}"')
    if tissue:
        parts.append(f'"{tissue}"')
    if date_from or date_to:
        df = date_from.replace("-", "/") if date_from else "1900/01/01"
        dt = date_to.replace("-", "/") if date_to else "3000/12/31"
        parts.append(f'"{df}"[PDAT] : "{dt}"[PDAT]')

    return " AND ".join(parts)


def _get_item_text(doc: ET.Element, name: str) -> str:
    """Extract text from a GEO eSummary Item element by Name attribute."""
    for item in doc.findall("Item"):
        if item.get("Name") == name:
            return (item.text or "").strip()
    return ""


def _get_item_int(doc: ET.Element, name: str) -> int:
    """Extract integer from a GEO eSummary Item element."""
    text = _get_item_text(doc, name)
    try:
        return int(text)
    except (ValueError, TypeError):
        return 0


def _get_pubmed_ids(doc: ET.Element) -> list[str]:
    """Extract PubMed IDs from eSummary."""
    for item in doc.findall("Item"):
        if item.get("Name") == "PubMedIds":
            return [sub.text.strip() for sub in item.findall("Item") if sub.text]
    return []


def _infer_platform(gpl: str, gds_type: str) -> str:
    """Infer sequencing platform from GPL ID and GDS type."""
    if not gpl:
        return "unknown"
    return gpl


def _infer_data_files(supp_file: str) -> list[dict]:
    """Infer available data files from supplementary file info."""
    files = []
    if not supp_file:
        return [{"type": "supplementary", "format": "unknown", "size_mb": 0}]
    for fmt in supp_file.split(";"):
        fmt = fmt.strip().upper()
        if fmt:
            files.append({"type": "count_matrix" if fmt in ("H5", "MTX", "CSV", "TSV") else "supplementary", "format": fmt.lower(), "size_mb": 0})
    return files if files else [{"type": "supplementary", "format": "unknown", "size_mb": 0}]


def _assess_metadata_quality(doc: ET.Element) -> str:
    """Assess metadata completeness."""
    has_title = bool(_get_item_text(doc, "title"))
    has_summary = bool(_get_item_text(doc, "summary"))
    has_samples = _get_item_int(doc, "n_samples") > 0
    has_pubmed = bool(_get_pubmed_ids(doc))
    score = sum([has_title, has_summary, has_samples, has_pubmed])
    if score >= 4:
        return "good"
    elif score >= 2:
        return "partial"
    return "minimal"


def parse_geo_record(doc: ET.Element, omic: str) -> dict:
    """Parse an eSummary DocSum element into a SearchResult dict."""
    accession = _get_item_text(doc, "Accession")
    if not accession:
        accession = f"GSE{_get_item_text(doc, 'GSE')}" if _get_item_text(doc, "GSE") else ""

    pdat = _get_item_text(doc, "PDAT").replace("/", "-")
    summary = _get_item_text(doc, "summary")

    # Paper info from PubMed IDs — will be resolved by the agent if needed
    pubmed_ids = _get_pubmed_ids(doc)
    paper = {
        "title": "",
        "authors": "",
        "journal": "",
        "doi": "",
        "date": "",
        "abstract": summary,
        "pubmed_id": pubmed_ids[0] if pubmed_ids else "",
    }

    return {
        "accession": accession,
        "source": "GEO",
        "title": _get_item_text(doc, "title"),
        "organism": _get_item_text(doc, "taxon"),
        "omic_type": omic,
        "platform": _infer_platform(_get_item_text(doc, "GPL"), _get_item_text(doc, "gdsType")),
        "disease": "",
        "tissue": "",
        "sample_count": _get_item_int(doc, "n_samples"),
        "condition_groups": [],
        "data_files": _infer_data_files(_get_item_text(doc, "suppFile")),
        "paper": paper,
        "metadata_quality": _assess_metadata_quality(doc),
        "date_submitted": pdat,
    }


def search_geo(
    omic: str,
    organism: str = "human",
    disease: str | None = None,
    tissue: str | None = None,
    max_results: int = 50,
    date_from: str | None = None,
    date_to: str | None = None,
) -> list[dict]:
    """Search GEO and return a list of SearchResult dicts."""
    query = build_esearch_query(omic, organism, disease, tissue, date_from, date_to)
    logger.info(f"GEO query: {query}")

    api_key = os.environ.get("NCBI_API_KEY", "")
    base_params = {"api_key": api_key} if api_key else {}

    # Step 1: esearch to get GEO IDs
    search_params = {
        **base_params,
        "db": "gds",
        "term": query,
        "retmax": max_results,
        "usehistory": "n",
    }
    resp = fetch_with_retry(ESEARCH_URL, params=search_params)
    root = ET.fromstring(resp.text)

    id_list = root.find("IdList")
    if id_list is None or not id_list.findall("Id"):
        logger.info("No GEO results found.")
        return []

    ids = [id_elem.text for id_elem in id_list.findall("Id") if id_elem.text]
    logger.info(f"Found {len(ids)} GEO IDs")

    # Step 2: esummary to get metadata
    summary_params = {
        **base_params,
        "db": "gds",
        "id": ",".join(ids),
    }
    resp = fetch_with_retry(ESUMMARY_URL, params=summary_params)
    root = ET.fromstring(resp.text)

    results = []
    for doc in root.findall("DocSum"):
        # Only include GSE series entries
        accession = _get_item_text(doc, "Accession")
        if accession and accession.startswith("GSE"):
            result = parse_geo_record(doc, omic)
            results.append(result)

    logger.info(f"Parsed {len(results)} GEO series records")
    return results


def main():
    parser = argparse.ArgumentParser(description="Search NCBI GEO for omics datasets")
    parser.add_argument("--omic", required=True, help="Omic type (scRNAseq, bulkRNAseq, etc.)")
    parser.add_argument("--organism", default="human", help="Organism (human, mouse)")
    parser.add_argument("--disease", default=None, help="Disease filter")
    parser.add_argument("--tissue", default=None, help="Tissue filter")
    parser.add_argument("--max-results", type=int, default=50, help="Max results")
    parser.add_argument("--date-from", default=None, help="Start date (YYYY-MM-DD)")
    parser.add_argument("--date-to", default=None, help="End date (YYYY-MM-DD)")
    args = parser.parse_args()

    results = search_geo(
        omic=args.omic, organism=args.organism, disease=args.disease,
        tissue=args.tissue, max_results=args.max_results,
        date_from=args.date_from, date_to=args.date_to,
    )
    json.dump(results, sys.stdout, indent=2)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
conda run -n dataseek pytest tests/test_search_geo.py -v
```

Expected: All PASS

- [ ] **Step 5: Commit**

```bash
git add scripts/search_geo.py tests/test_search_geo.py
git commit -m "feat: GEO search via Entrez E-utilities with query building and XML parsing"
```

---

## Task 4: CCLE Search Script

**Files:**
- Create: `scripts/search_ccle.py`
- Create: `tests/test_search_ccle.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_search_ccle.py
import json

import pytest
import responses

from scripts.search_ccle import search_ccle, parse_ccle_dataset
from tests.conftest import SEARCH_RESULT_REQUIRED_KEYS


SAMPLE_DEPMAP_DOWNLOAD_LIST = [
    {
        "fileName": "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
        "fileDescription": "Gene expression TPM values (log2(TPM+1)) for protein coding genes",
        "releaseName": "DepMap Public 24Q4",
        "tagsUrl": "expression",
        "size": 450000000,
    },
    {
        "fileName": "OmicsCNGene.csv",
        "fileDescription": "Gene-level copy number data",
        "releaseName": "DepMap Public 24Q4",
        "tagsUrl": "copy_number",
        "size": 120000000,
    },
    {
        "fileName": "OmicsSomaticMutations.csv",
        "fileDescription": "Somatic mutation calls",
        "releaseName": "DepMap Public 24Q4",
        "tagsUrl": "mutation",
        "size": 80000000,
    },
]


class TestParseCcleDataset:
    def test_parses_download_entry(self):
        entry = SAMPLE_DEPMAP_DOWNLOAD_LIST[0]
        result = parse_ccle_dataset(entry, omic="bulkRNAseq")
        assert result["source"] == "CCLE"
        assert result["organism"] == "Homo sapiens"
        for key in SEARCH_RESULT_REQUIRED_KEYS:
            assert key in result


class TestSearchCcle:
    @responses.activate
    def test_returns_matching_datasets(self):
        responses.add(
            responses.GET,
            "https://depmap.org/portal/api/download/all",
            json=SAMPLE_DEPMAP_DOWNLOAD_LIST,
            status=200,
        )
        results = search_ccle(omic="bulkRNAseq")
        assert isinstance(results, list)
        assert len(results) >= 1
        assert all(r["source"] == "CCLE" for r in results)

    @responses.activate
    def test_filters_by_omic_type(self):
        responses.add(
            responses.GET,
            "https://depmap.org/portal/api/download/all",
            json=SAMPLE_DEPMAP_DOWNLOAD_LIST,
            status=200,
        )
        results = search_ccle(omic="bulkRNAseq")
        # Should match expression data
        assert any("expression" in r["title"].lower() or "expression" in r.get("accession", "").lower() for r in results)
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
conda run -n dataseek pytest tests/test_search_ccle.py -v
```

- [ ] **Step 3: Implement search_ccle.py**

```python
#!/usr/bin/env python3
"""Search Cancer Cell Line Encyclopedia (CCLE) via DepMap Portal API."""

import argparse
import json
import sys

from scripts.utils import fetch_with_retry, setup_logger

DEPMAP_DOWNLOAD_API = "https://depmap.org/portal/api/download/all"

OMIC_TO_CCLE_TAGS = {
    "bulkRNAseq": ["expression"],
    "bulkGenomicseq": ["copy_number", "mutation", "wgs", "wes"],
}

logger = setup_logger("search_ccle")


def parse_ccle_dataset(entry: dict, omic: str) -> dict:
    """Parse a DepMap download entry into a SearchResult dict."""
    size_mb = round(entry.get("size", 0) / 1_000_000, 1)
    filename = entry.get("fileName", "")
    release = entry.get("releaseName", "")

    return {
        "accession": f"CCLE_{filename.replace('.csv', '').replace('.', '_')}",
        "source": "CCLE",
        "title": entry.get("fileDescription", filename),
        "organism": "Homo sapiens",
        "omic_type": omic,
        "platform": "multiple (cell line panel)",
        "disease": "cancer cell lines",
        "tissue": "multiple",
        "sample_count": 0,
        "condition_groups": [],
        "data_files": [{"type": "processed_matrix", "format": "csv", "size_mb": size_mb}],
        "paper": {
            "title": "Cancer Cell Line Encyclopedia",
            "authors": "Broad Institute",
            "journal": "DepMap Portal",
            "doi": "",
            "date": "",
            "abstract": f"Release: {release}. {entry.get('fileDescription', '')}",
        },
        "metadata_quality": "good",
        "date_submitted": "",
    }


def search_ccle(
    omic: str,
    disease: str | None = None,
    max_results: int = 50,
) -> list[dict]:
    """Search CCLE datasets on DepMap Portal. Returns list of SearchResult dicts."""
    tags = OMIC_TO_CCLE_TAGS.get(omic, [])
    if not tags:
        logger.info(f"No CCLE mapping for omic type: {omic}")
        return []

    logger.info(f"Querying DepMap download API for tags: {tags}")
    resp = fetch_with_retry(DEPMAP_DOWNLOAD_API)
    if resp.status_code != 200:
        logger.warning(f"DepMap API returned status {resp.status_code}")
        return []

    all_entries = resp.json()
    results = []
    for entry in all_entries:
        entry_tags = entry.get("tagsUrl", "")
        if any(tag in entry_tags for tag in tags):
            result = parse_ccle_dataset(entry, omic)
            if disease and disease.lower() not in result["title"].lower():
                continue
            results.append(result)

    logger.info(f"Found {len(results)} CCLE datasets")
    return results[:max_results]


def main():
    parser = argparse.ArgumentParser(description="Search CCLE via DepMap Portal")
    parser.add_argument("--omic", required=True)
    parser.add_argument("--disease", default=None)
    parser.add_argument("--max-results", type=int, default=50)
    args = parser.parse_args()

    results = search_ccle(omic=args.omic, disease=args.disease, max_results=args.max_results)
    json.dump(results, sys.stdout, indent=2)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
conda run -n dataseek pytest tests/test_search_ccle.py -v
```

- [ ] **Step 5: Commit**

```bash
git add scripts/search_ccle.py tests/test_search_ccle.py
git commit -m "feat: CCLE search via DepMap Portal download API"
```

---

## Task 5: Xena UCSC Search Script

**Files:**
- Create: `scripts/search_xena.py`
- Create: `tests/test_search_xena.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_search_xena.py
import json

import pytest

from scripts.search_xena import search_xena, parse_xena_dataset, filter_datasets_by_omic
from tests.conftest import SEARCH_RESULT_REQUIRED_KEYS


SAMPLE_XENA_COHORT = {
    "name": "TCGA-GBM",
    "label": "TCGA Glioblastoma (GBM)",
    "host": "https://tcga.xenahubs.net",
}

SAMPLE_XENA_DATASET = {
    "name": "TCGA-GBM.htseq_fpkm.tsv",
    "label": "Gene Expression RNAseq - HTSeq - FPKM",
    "cohort": "TCGA Glioblastoma (GBM)",
    "type": "genomicMatrix",
    "dataSubType": "gene expression RNAseq",
    "unit": "fpkm",
    "host": "https://tcga.xenahubs.net",
    "sampleCount": 174,
}


class TestParsXenaDataset:
    def test_parses_to_search_result(self):
        result = parse_xena_dataset(SAMPLE_XENA_DATASET, omic="bulkRNAseq")
        assert result["source"] == "Xena"
        assert result["sample_count"] == 174
        assert "TCGA" in result["title"] or "TCGA" in result["accession"]
        for key in SEARCH_RESULT_REQUIRED_KEYS:
            assert key in result


class TestFilterDatasetsByOmic:
    def test_filters_expression_for_bulkRNAseq(self):
        datasets = [
            {"dataSubType": "gene expression RNAseq", "label": "RNA"},
            {"dataSubType": "copy number", "label": "CN"},
            {"dataSubType": "somatic mutation", "label": "Mut"},
        ]
        filtered = filter_datasets_by_omic(datasets, "bulkRNAseq")
        assert len(filtered) == 1
        assert filtered[0]["label"] == "RNA"

    def test_filters_cnv_for_bulkGenomicseq(self):
        datasets = [
            {"dataSubType": "gene expression RNAseq", "label": "RNA"},
            {"dataSubType": "copy number", "label": "CN"},
        ]
        filtered = filter_datasets_by_omic(datasets, "bulkGenomicseq")
        assert len(filtered) == 1
        assert filtered[0]["label"] == "CN"
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
conda run -n dataseek pytest tests/test_search_xena.py -v
```

- [ ] **Step 3: Implement search_xena.py**

```python
#!/usr/bin/env python3
"""Search UCSC Xena for omics datasets via xenaPython."""

import argparse
import json
import sys

try:
    import xenaPython as xena
except ImportError:
    xena = None

from scripts.utils import fetch_with_retry, setup_logger

XENA_HUB_URL = "https://ucscpublic.xenahubs.net"
TCGA_HUB_URL = "https://tcga.xenahubs.net"
TOIL_HUB_URL = "https://toil.xenahubs.net"

ALL_HUBS = [TCGA_HUB_URL, XENA_HUB_URL, TOIL_HUB_URL]

OMIC_TO_XENA_SUBTYPES = {
    "bulkRNAseq": ["gene expression rnaseq", "gene expression", "rnaseq"],
    "bulkGenomicseq": ["copy number", "somatic mutation", "mutation", "cnv"],
}

logger = setup_logger("search_xena")


def filter_datasets_by_omic(datasets: list[dict], omic: str) -> list[dict]:
    """Filter Xena datasets by omic-relevant data subtypes."""
    subtypes = OMIC_TO_XENA_SUBTYPES.get(omic, [])
    if not subtypes:
        return []
    return [
        ds for ds in datasets
        if any(st in ds.get("dataSubType", "").lower() for st in subtypes)
    ]


def parse_xena_dataset(ds: dict, omic: str) -> dict:
    """Parse a Xena dataset metadata dict into a SearchResult."""
    cohort = ds.get("cohort", "")
    name = ds.get("name", "")
    label = ds.get("label", "")

    return {
        "accession": f"XENA_{name.replace('/', '_').replace('.', '_')}",
        "source": "Xena",
        "title": f"{cohort} - {label}",
        "organism": "Homo sapiens",
        "omic_type": omic,
        "platform": ds.get("dataSubType", "unknown"),
        "disease": cohort,
        "tissue": "",
        "sample_count": ds.get("sampleCount", 0),
        "condition_groups": [],
        "data_files": [{"type": "genomic_matrix", "format": "tsv", "size_mb": 0}],
        "paper": {
            "title": cohort,
            "authors": "",
            "journal": "",
            "doi": "",
            "date": "",
            "abstract": f"Xena dataset: {label}. Cohort: {cohort}.",
        },
        "metadata_quality": "good" if ds.get("sampleCount", 0) > 0 else "partial",
        "date_submitted": "",
    }


def search_xena(
    omic: str,
    organism: str = "human",
    disease: str | None = None,
    tissue: str | None = None,
    max_results: int = 50,
) -> list[dict]:
    """Search UCSC Xena hubs for matching datasets. Returns list of SearchResult dicts."""
    if xena is None:
        logger.warning("xenaPython not installed, falling back to API")
        return _search_xena_api(omic, disease, max_results)

    results = []
    for hub in ALL_HUBS:
        try:
            datasets = xena.all_datasets(hub)
            if not datasets:
                continue
            filtered = filter_datasets_by_omic(datasets, omic)
            for ds in filtered:
                ds["host"] = hub
                if disease and disease.lower() not in ds.get("cohort", "").lower():
                    continue
                results.append(parse_xena_dataset(ds, omic))
        except Exception as e:
            logger.warning(f"Failed to query Xena hub {hub}: {e}")

    logger.info(f"Found {len(results)} Xena datasets")
    return results[:max_results]


def _search_xena_api(omic: str, disease: str | None, max_results: int) -> list[dict]:
    """Fallback: query Xena via REST API if xenaPython unavailable."""
    results = []
    for hub in ALL_HUBS:
        try:
            resp = fetch_with_retry(f"{hub}/data/", params={"format": "json"})
            if resp.status_code != 200:
                continue
            datasets = resp.json()
            filtered = filter_datasets_by_omic(datasets, omic)
            for ds in filtered:
                if disease and disease.lower() not in ds.get("cohort", "").lower():
                    continue
                results.append(parse_xena_dataset(ds, omic))
        except Exception as e:
            logger.warning(f"Xena API fallback failed for {hub}: {e}")
    return results[:max_results]


def main():
    parser = argparse.ArgumentParser(description="Search UCSC Xena for omics datasets")
    parser.add_argument("--omic", required=True)
    parser.add_argument("--organism", default="human")
    parser.add_argument("--disease", default=None)
    parser.add_argument("--tissue", default=None)
    parser.add_argument("--max-results", type=int, default=50)
    args = parser.parse_args()

    results = search_xena(omic=args.omic, organism=args.organism, disease=args.disease,
                          tissue=args.tissue, max_results=args.max_results)
    json.dump(results, sys.stdout, indent=2)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
conda run -n dataseek pytest tests/test_search_xena.py -v
```

- [ ] **Step 5: Commit**

```bash
git add scripts/search_xena.py tests/test_search_xena.py
git commit -m "feat: Xena UCSC search with xenaPython and REST API fallback"
```

---

## Task 6: DepMap Search Script

**Files:**
- Create: `scripts/search_depmap.py`
- Create: `tests/test_search_depmap.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_search_depmap.py
import json

import pytest
import responses

from scripts.search_depmap import search_depmap, parse_depmap_dataset
from tests.conftest import SEARCH_RESULT_REQUIRED_KEYS


SAMPLE_DEPMAP_DATASETS = [
    {
        "fileName": "CRISPRGeneEffect.csv",
        "fileDescription": "CRISPR knockout gene effect scores (Chronos)",
        "releaseName": "DepMap Public 24Q4",
        "tagsUrl": "crispr",
        "size": 300000000,
    },
    {
        "fileName": "CRISPRGeneDependency.csv",
        "fileDescription": "CRISPR gene dependency probability scores",
        "releaseName": "DepMap Public 24Q4",
        "tagsUrl": "crispr",
        "size": 250000000,
    },
    {
        "fileName": "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
        "fileDescription": "Gene expression TPM values",
        "releaseName": "DepMap Public 24Q4",
        "tagsUrl": "expression",
        "size": 450000000,
    },
]


class TestParseDepmapDataset:
    def test_parses_to_search_result(self):
        result = parse_depmap_dataset(SAMPLE_DEPMAP_DATASETS[0], omic="CRISPR")
        assert result["source"] == "DepMap"
        for key in SEARCH_RESULT_REQUIRED_KEYS:
            assert key in result


class TestSearchDepmap:
    @responses.activate
    def test_filters_crispr_datasets(self):
        responses.add(
            responses.GET,
            "https://depmap.org/portal/api/download/all",
            json=SAMPLE_DEPMAP_DATASETS,
            status=200,
        )
        results = search_depmap(omic="CRISPR")
        assert len(results) >= 1
        assert all(r["source"] == "DepMap" for r in results)

    @responses.activate
    def test_returns_empty_for_unmapped_omic(self):
        responses.add(
            responses.GET,
            "https://depmap.org/portal/api/download/all",
            json=SAMPLE_DEPMAP_DATASETS,
            status=200,
        )
        results = search_depmap(omic="ChIPseq")
        assert results == []
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
conda run -n dataseek pytest tests/test_search_depmap.py -v
```

- [ ] **Step 3: Implement search_depmap.py**

```python
#!/usr/bin/env python3
"""Search DepMap portal for functional genomics datasets."""

import argparse
import json
import sys

from scripts.utils import fetch_with_retry, setup_logger

DEPMAP_DOWNLOAD_API = "https://depmap.org/portal/api/download/all"

OMIC_TO_DEPMAP_TAGS = {
    "CRISPR": ["crispr"],
    "bulkGenomicseq": ["copy_number", "mutation", "wgs", "wes", "fusion"],
}

logger = setup_logger("search_depmap")


def parse_depmap_dataset(entry: dict, omic: str) -> dict:
    """Parse a DepMap download entry into a SearchResult dict."""
    size_mb = round(entry.get("size", 0) / 1_000_000, 1)
    filename = entry.get("fileName", "")
    release = entry.get("releaseName", "")

    return {
        "accession": f"DepMap_{filename.replace('.csv', '').replace('.', '_')}",
        "source": "DepMap",
        "title": entry.get("fileDescription", filename),
        "organism": "Homo sapiens",
        "omic_type": omic,
        "platform": "DepMap (Broad Institute)",
        "disease": "cancer cell lines",
        "tissue": "multiple",
        "sample_count": 0,
        "condition_groups": [],
        "data_files": [{"type": "processed_matrix", "format": "csv", "size_mb": size_mb}],
        "paper": {
            "title": "Cancer Dependency Map",
            "authors": "Broad Institute",
            "journal": "DepMap Portal",
            "doi": "",
            "date": "",
            "abstract": f"Release: {release}. {entry.get('fileDescription', '')}",
        },
        "metadata_quality": "good",
        "date_submitted": "",
    }


def search_depmap(
    omic: str,
    disease: str | None = None,
    max_results: int = 50,
) -> list[dict]:
    """Search DepMap for matching datasets. Returns list of SearchResult dicts."""
    tags = OMIC_TO_DEPMAP_TAGS.get(omic, [])
    if not tags:
        logger.info(f"No DepMap mapping for omic type: {omic}")
        return []

    logger.info(f"Querying DepMap API for tags: {tags}")
    resp = fetch_with_retry(DEPMAP_DOWNLOAD_API)
    if resp.status_code != 200:
        logger.warning(f"DepMap API returned status {resp.status_code}")
        return []

    all_entries = resp.json()
    results = []
    for entry in all_entries:
        entry_tags = entry.get("tagsUrl", "")
        if any(tag in entry_tags for tag in tags):
            result = parse_depmap_dataset(entry, omic)
            results.append(result)

    logger.info(f"Found {len(results)} DepMap datasets")
    return results[:max_results]


def main():
    parser = argparse.ArgumentParser(description="Search DepMap for omics datasets")
    parser.add_argument("--omic", required=True)
    parser.add_argument("--disease", default=None)
    parser.add_argument("--max-results", type=int, default=50)
    args = parser.parse_args()

    results = search_depmap(omic=args.omic, disease=args.disease, max_results=args.max_results)
    json.dump(results, sys.stdout, indent=2)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
conda run -n dataseek pytest tests/test_search_depmap.py -v
```

- [ ] **Step 5: Commit**

```bash
git add scripts/search_depmap.py tests/test_search_depmap.py
git commit -m "feat: DepMap search for CRISPR screens and genomic datasets"
```

---

## Task 7: Single Cell Portal Search Script

**Files:**
- Create: `scripts/search_scp.py`
- Create: `tests/test_search_scp.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_search_scp.py
import json

import pytest
import responses

from scripts.search_scp import search_scp, parse_scp_study
from tests.conftest import SEARCH_RESULT_REQUIRED_KEYS


SAMPLE_SCP_RESPONSE = {
    "studies": [
        {
            "accession": "SCP1234",
            "name": "Single-cell atlas of glioblastoma",
            "description": "Comprehensive single-cell transcriptomic atlas of glioblastoma tumors",
            "cell_count": 50000,
            "gene_count": 20000,
            "disease": ["glioblastoma"],
            "organ": ["brain"],
            "species": ["Homo sapiens"],
            "library_preparation_protocol": ["10x 3' v3"],
            "study_url": "https://singlecell.broadinstitute.org/single_cell/study/SCP1234",
        }
    ],
    "total_studies": 1,
}


class TestParsScpStudy:
    def test_parses_to_search_result(self):
        study = SAMPLE_SCP_RESPONSE["studies"][0]
        result = parse_scp_study(study, omic="scRNAseq")
        assert result["accession"] == "SCP1234"
        assert result["source"] == "SCP"
        assert result["organism"] == "Homo sapiens"
        for key in SEARCH_RESULT_REQUIRED_KEYS:
            assert key in result


class TestSearchScp:
    @responses.activate
    def test_returns_matching_studies(self):
        responses.add(
            responses.GET,
            "https://singlecell.broadinstitute.org/single_cell/api/v1/search",
            json=SAMPLE_SCP_RESPONSE,
            status=200,
        )
        results = search_scp(omic="scRNAseq", disease="glioblastoma", organism="human")
        assert isinstance(results, list)
        assert len(results) >= 1
        assert results[0]["source"] == "SCP"

    @responses.activate
    def test_returns_empty_on_no_results(self):
        responses.add(
            responses.GET,
            "https://singlecell.broadinstitute.org/single_cell/api/v1/search",
            json={"studies": [], "total_studies": 0},
            status=200,
        )
        results = search_scp(omic="scRNAseq", disease="nonexistent")
        assert results == []
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
conda run -n dataseek pytest tests/test_search_scp.py -v
```

- [ ] **Step 3: Implement search_scp.py**

```python
#!/usr/bin/env python3
"""Search Single Cell Portal (Broad Institute) for single-cell datasets."""

import argparse
import json
import os
import sys

from scripts.utils import fetch_with_retry, setup_logger

SCP_API_BASE = "https://singlecell.broadinstitute.org/single_cell/api/v1"

logger = setup_logger("search_scp")


def parse_scp_study(study: dict, omic: str) -> dict:
    """Parse an SCP study dict into a SearchResult."""
    species = study.get("species", [])
    organism = species[0] if species else ""
    organs = study.get("organ", [])
    diseases = study.get("disease", [])
    protocols = study.get("library_preparation_protocol", [])

    return {
        "accession": study.get("accession", ""),
        "source": "SCP",
        "title": study.get("name", ""),
        "organism": organism,
        "omic_type": omic,
        "platform": ", ".join(protocols) if protocols else "unknown",
        "disease": ", ".join(diseases) if diseases else "",
        "tissue": ", ".join(organs) if organs else "",
        "sample_count": 0,
        "condition_groups": [],
        "data_files": [{"type": "study_files", "format": "mixed", "size_mb": 0}],
        "paper": {
            "title": "",
            "authors": "",
            "journal": "",
            "doi": "",
            "date": "",
            "abstract": study.get("description", ""),
        },
        "metadata_quality": "good" if study.get("description") else "partial",
        "date_submitted": "",
    }


def search_scp(
    omic: str,
    organism: str = "human",
    disease: str | None = None,
    tissue: str | None = None,
    max_results: int = 50,
) -> list[dict]:
    """Search Single Cell Portal. Returns list of SearchResult dicts."""
    token = os.environ.get("SCP_TOKEN", "")
    headers = {}
    if token:
        headers["Authorization"] = f"Bearer {token}"

    # Build search terms
    terms = []
    organism_map = {"human": "Homo sapiens", "mouse": "Mus musculus"}
    org = organism_map.get(organism.lower(), organism)
    terms.append(org)
    if disease:
        terms.append(disease)
    if tissue:
        terms.append(tissue)

    params = {
        "terms": " ".join(terms),
        "limit": max_results,
    }

    logger.info(f"Searching SCP: {params}")
    resp = fetch_with_retry(f"{SCP_API_BASE}/search", params=params, headers=headers)
    if resp.status_code != 200:
        logger.warning(f"SCP API returned status {resp.status_code}")
        return []

    data = resp.json()
    studies = data.get("studies", [])

    results = []
    for study in studies:
        result = parse_scp_study(study, omic)
        results.append(result)

    logger.info(f"Found {len(results)} SCP studies")
    return results


def main():
    parser = argparse.ArgumentParser(description="Search Single Cell Portal")
    parser.add_argument("--omic", required=True)
    parser.add_argument("--organism", default="human")
    parser.add_argument("--disease", default=None)
    parser.add_argument("--tissue", default=None)
    parser.add_argument("--max-results", type=int, default=50)
    args = parser.parse_args()

    results = search_scp(omic=args.omic, organism=args.organism, disease=args.disease,
                         tissue=args.tissue, max_results=args.max_results)
    json.dump(results, sys.stdout, indent=2)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
conda run -n dataseek pytest tests/test_search_scp.py -v
```

- [ ] **Step 5: Commit**

```bash
git add scripts/search_scp.py tests/test_search_scp.py
git commit -m "feat: Single Cell Portal search via REST API"
```

---

## Task 8: GEO Download Script

**Files:**
- Create: `scripts/download_geo.py`
- Create: `tests/test_download_geo.py`

- [ ] **Step 1: Write failing tests**

```python
# tests/test_download_geo.py
import json
import os

import pytest
import responses

from scripts.download_geo import download_geo, build_geo_ftp_urls, generate_summary_report


SAMPLE_GEO_MINIML = """<?xml version="1.0" encoding="UTF-8"?>
<MINiML xmlns="http://www.ncbi.nlm.nih.gov/geo/info/MINiML" version="0.5.0">
    <Series iid="GSE123456">
        <Title>Single-cell RNA-seq of glioblastoma</Title>
        <Summary>We performed scRNAseq on glioblastoma tumors.</Summary>
        <Overall-Design>12 tumor samples, 6 normal controls</Overall-Design>
        <Pubmed-ID>38000001</Pubmed-ID>
        <Supplementary-Data type="unknown">
            ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123456/suppl/GSE123456_RAW.tar
        </Supplementary-Data>
    </Series>
    <Sample iid="GSM100001">
        <Title>Tumor_sample_1</Title>
        <Channel position="1">
            <Source>brain tumor</Source>
            <Organism>Homo sapiens</Organism>
        </Channel>
    </Sample>
    <Sample iid="GSM100002">
        <Title>Normal_sample_1</Title>
        <Channel position="1">
            <Source>normal brain</Source>
            <Organism>Homo sapiens</Organism>
        </Channel>
    </Sample>
</MINiML>"""


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
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
conda run -n dataseek pytest tests/test_download_geo.py -v
```

- [ ] **Step 3: Implement download_geo.py**

```python
#!/usr/bin/env python3
"""Download GEO datasets with metadata and supplementary files."""

import argparse
import json
import os
import sys
import xml.etree.ElementTree as ET

from scripts.utils import fetch_with_retry, download_file, resolve_doi, setup_logger

GEO_MINIML_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/{range}/{accession}/miniml/{accession}_family.xml.tgz"
GEO_SOFT_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
GEO_SUPP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/{range}/{accession}/suppl/"

logger = setup_logger("download_geo")

MINIML_NS = {"geo": "http://www.ncbi.nlm.nih.gov/geo/info/MINiML"}


def _accession_range(accession: str) -> str:
    """Convert GSE123456 -> GSE123nnn for FTP path."""
    numeric = accession.replace("GSE", "")
    return f"GSE{numeric[:-3]}nnn"


def build_geo_ftp_urls(accession: str) -> list[str]:
    """Build FTP/HTTPS URLs for GEO supplementary data."""
    acc_range = _accession_range(accession)
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{acc_range}/{accession}/suppl/"
    return [base, f"{base}{accession}_RAW.tar"]


def _fetch_geo_metadata(accession: str) -> dict:
    """Fetch GEO series metadata via SOFT format."""
    params = {"acc": accession, "targ": "self", "form": "text", "view": "brief"}
    resp = fetch_with_retry(GEO_SOFT_URL, params=params)
    if resp.status_code != 200:
        return {}

    metadata = {"title": "", "summary": "", "overall_design": "", "pubmed_ids": [], "samples": [], "platform": "", "organism": ""}
    for line in resp.text.split("\n"):
        if line.startswith("!Series_title"):
            metadata["title"] = line.split("=", 1)[1].strip() if "=" in line else ""
        elif line.startswith("!Series_summary"):
            metadata["summary"] += line.split("=", 1)[1].strip() + " " if "=" in line else ""
        elif line.startswith("!Series_overall_design"):
            metadata["overall_design"] = line.split("=", 1)[1].strip() if "=" in line else ""
        elif line.startswith("!Series_pubmed_id"):
            val = line.split("=", 1)[1].strip() if "=" in line else ""
            if val:
                metadata["pubmed_ids"].append(val)
        elif line.startswith("!Series_sample_id"):
            val = line.split("=", 1)[1].strip() if "=" in line else ""
            if val:
                metadata["samples"].append(val)
        elif line.startswith("!Series_platform_id"):
            metadata["platform"] = line.split("=", 1)[1].strip() if "=" in line else ""
        elif line.startswith("!Platform_organism"):
            metadata["organism"] = line.split("=", 1)[1].strip() if "=" in line else ""

    return metadata


def _fetch_supplementary_file_list(accession: str) -> list[dict]:
    """List supplementary files available for a GEO accession."""
    acc_range = _accession_range(accession)
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{acc_range}/{accession}/suppl/"
    try:
        resp = fetch_with_retry(url)
        if resp.status_code != 200:
            return []
        # Parse the FTP directory listing (HTML)
        files = []
        for line in resp.text.split("\n"):
            if 'href="' in line and accession in line:
                href_start = line.index('href="') + 6
                href_end = line.index('"', href_start)
                filename = line[href_start:href_end]
                if filename and not filename.endswith("/"):
                    files.append({"name": filename, "url": url + filename, "size_mb": 0})
        return files
    except Exception as e:
        logger.warning(f"Failed to list supplementary files: {e}")
        return []


def generate_summary_report(
    accession: str, title: str, organism: str, omic_type: str,
    platform: str, disease: str, tissue: str, sample_count: int,
    condition_groups: list, data_files: list, paper: dict,
    metadata_quality: str,
) -> str:
    """Generate a markdown summary report for a downloaded dataset."""
    conditions = ", ".join(condition_groups) if condition_groups else "unknown"
    files_list = "\n".join(f"  - {f.get('name', 'unknown')} ({f.get('size_mb', 0)} MB)" for f in data_files)

    return f"""# {accession}

## Dataset Metadata
- **Accession:** {accession}
- **Title:** {title}
- **Organism:** {organism}
- **Tissue:** {tissue or 'not specified'}
- **Disease/Condition:** {disease or 'not specified'}
- **Sample Count:** {sample_count}

## Experimental Design
- **Omic Type:** {omic_type}
- **Platform:** {platform}
- **Conditions/Groups:** {conditions}

## Data Availability
{files_list if files_list.strip() else '  - No supplementary files listed'}

## Original Paper
- **Title:** {paper.get('title', 'N/A')}
- **Authors:** {paper.get('authors', 'N/A')}
- **Journal:** {paper.get('journal', 'N/A')}
- **DOI:** {paper.get('doi', 'N/A')}
- **Date:** {paper.get('date', 'N/A')}
- **Abstract:** {paper.get('abstract', 'N/A')}

## Quality Indicators
- **Metadata Quality:** {metadata_quality}
- **Samples per Group:** {conditions}
"""


def download_geo(accession: str, output_dir: str, omic_type: str = "unknown") -> dict:
    """Download GEO dataset files and metadata. Returns download summary."""
    os.makedirs(output_dir, exist_ok=True)
    data_dir = os.path.join(output_dir, "data")
    meta_dir = os.path.join(output_dir, "metadata")
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(meta_dir, exist_ok=True)

    logger.info(f"Downloading GEO dataset: {accession}")

    # 1. Fetch metadata
    metadata = _fetch_geo_metadata(accession)
    with open(os.path.join(meta_dir, "study_metadata.json"), "w") as f:
        json.dump(metadata, f, indent=2)

    # 2. Resolve paper info
    paper = {"title": "", "authors": "", "journal": "", "doi": "", "date": "", "abstract": metadata.get("summary", "")}
    if metadata.get("pubmed_ids"):
        # Try to get DOI from PubMed
        pmid = metadata["pubmed_ids"][0]
        try:
            resp = fetch_with_retry(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                params={"db": "pubmed", "id": pmid, "retmode": "json"},
            )
            if resp.status_code == 200:
                pub_data = resp.json().get("result", {}).get(pmid, {})
                doi = ""
                for aid in pub_data.get("articleids", []):
                    if aid.get("idtype") == "doi":
                        doi = aid.get("value", "")
                        break
                if doi:
                    paper = resolve_doi(doi)
                else:
                    paper["title"] = pub_data.get("title", "")
                    paper["authors"] = ", ".join(a.get("name", "") for a in pub_data.get("authors", []))
                    paper["journal"] = pub_data.get("fulljournalname", "")
                    paper["date"] = pub_data.get("pubdate", "")
        except Exception as e:
            logger.warning(f"Failed to resolve paper for PMID {pmid}: {e}")

    with open(os.path.join(meta_dir, "paper.json"), "w") as f:
        json.dump(paper, f, indent=2)

    # 3. List and download supplementary files
    supp_files = _fetch_supplementary_file_list(accession)
    downloaded_files = []
    for sf in supp_files:
        dest = os.path.join(data_dir, sf["name"])
        try:
            logger.info(f"Downloading {sf['name']}...")
            download_file(sf["url"], dest)
            size_mb = round(os.path.getsize(dest) / 1_000_000, 1)
            downloaded_files.append({"name": sf["name"], "size_mb": size_mb})
        except Exception as e:
            logger.warning(f"Failed to download {sf['name']}: {e}")

    # 4. Generate sample metadata CSV
    sample_metadata_path = os.path.join(meta_dir, "sample_metadata.csv")
    with open(sample_metadata_path, "w") as f:
        f.write("sample_id,condition,batch\n")
        for sample_id in metadata.get("samples", []):
            f.write(f"{sample_id},TODO,TODO\n")

    return {
        "accession": accession,
        "output_dir": output_dir,
        "files_downloaded": len(downloaded_files),
        "data_files": downloaded_files,
        "paper": paper,
        "metadata": metadata,
    }


def main():
    parser = argparse.ArgumentParser(description="Download GEO dataset")
    parser.add_argument("--accession", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--omic-type", default="unknown")
    args = parser.parse_args()

    result = download_geo(args.accession, args.output_dir, args.omic_type)
    json.dump(result, sys.stdout, indent=2)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
conda run -n dataseek pytest tests/test_download_geo.py -v
```

- [ ] **Step 5: Commit**

```bash
git add scripts/download_geo.py tests/test_download_geo.py
git commit -m "feat: GEO dataset downloader with metadata extraction and summary reports"
```

---

## Task 9: CCLE, Xena, DepMap, SCP Download Scripts

**Files:**
- Create: `scripts/download_ccle.py`, `scripts/download_xena.py`, `scripts/download_depmap.py`, `scripts/download_scp.py`
- Create: `tests/test_download_ccle.py`, `tests/test_download_xena.py`, `tests/test_download_depmap.py`, `tests/test_download_scp.py`

All four download scripts share the same pattern: fetch file by URL, save metadata, generate summary. They differ only in URL construction and metadata parsing.

- [ ] **Step 1: Write failing tests for all four downloaders**

```python
# tests/test_download_ccle.py
import json
import os
import pytest
import responses
from scripts.download_ccle import download_ccle

class TestDownloadCcle:
    @responses.activate
    def test_downloads_to_output_dir(self, tmp_path):
        output_dir = str(tmp_path / "CCLE_test")
        responses.add(
            responses.GET,
            "https://depmap.org/portal/api/download/all",
            json=[{"fileName": "OmicsExpression.csv", "downloadUrl": "https://depmap.org/portal/api/download/file/OmicsExpression.csv", "size": 100}],
            status=200,
        )
        responses.add(
            responses.GET,
            "https://depmap.org/portal/api/download/file/OmicsExpression.csv",
            body=b"gene,sample1\nTP53,5.2",
            status=200,
        )
        result = download_ccle("CCLE_OmicsExpression", output_dir)
        assert result["accession"] == "CCLE_OmicsExpression"
        assert os.path.isdir(output_dir)
```

```python
# tests/test_download_xena.py
import json
import os
import pytest
import responses
from scripts.download_xena import download_xena

class TestDownloadXena:
    def test_builds_download_url(self):
        from scripts.download_xena import build_xena_download_url
        url = build_xena_download_url("https://tcga.xenahubs.net", "TCGA-GBM.htseq_fpkm.tsv")
        assert "tcga.xenahubs.net" in url
        assert "TCGA-GBM" in url
```

```python
# tests/test_download_depmap.py
import json
import os
import pytest
import responses
from scripts.download_depmap import download_depmap

class TestDownloadDepmap:
    @responses.activate
    def test_downloads_to_output_dir(self, tmp_path):
        output_dir = str(tmp_path / "DepMap_test")
        responses.add(
            responses.GET,
            "https://depmap.org/portal/api/download/all",
            json=[{"fileName": "CRISPRGeneEffect.csv", "downloadUrl": "https://depmap.org/portal/api/download/file/CRISPRGeneEffect.csv", "size": 100}],
            status=200,
        )
        responses.add(
            responses.GET,
            "https://depmap.org/portal/api/download/file/CRISPRGeneEffect.csv",
            body=b"gene,cell1\nTP53,-0.5",
            status=200,
        )
        result = download_depmap("DepMap_CRISPRGeneEffect", output_dir)
        assert result["accession"] == "DepMap_CRISPRGeneEffect"
```

```python
# tests/test_download_scp.py
import json
import os
import pytest
import responses
from scripts.download_scp import download_scp

class TestDownloadScp:
    @responses.activate
    def test_downloads_study_files(self, tmp_path):
        output_dir = str(tmp_path / "SCP_test")
        study_info = {
            "accession": "SCP1234",
            "name": "Test Study",
            "description": "A test study",
            "study_files": [
                {"name": "expression_matrix.tsv.gz", "download_url": "https://singlecell.broadinstitute.org/data/SCP1234/expression_matrix.tsv.gz"}
            ],
        }
        responses.add(
            responses.GET,
            "https://singlecell.broadinstitute.org/single_cell/api/v1/studies/SCP1234",
            json=study_info,
            status=200,
        )
        responses.add(
            responses.GET,
            "https://singlecell.broadinstitute.org/data/SCP1234/expression_matrix.tsv.gz",
            body=b"fake compressed data",
            status=200,
        )
        result = download_scp("SCP1234", output_dir)
        assert result["accession"] == "SCP1234"
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
conda run -n dataseek pytest tests/test_download_ccle.py tests/test_download_xena.py tests/test_download_depmap.py tests/test_download_scp.py -v
```

- [ ] **Step 3: Implement download_ccle.py**

```python
#!/usr/bin/env python3
"""Download CCLE datasets from DepMap Portal."""

import argparse
import json
import os
import sys

from scripts.utils import fetch_with_retry, download_file, setup_logger
from scripts.download_geo import generate_summary_report

DEPMAP_DOWNLOAD_API = "https://depmap.org/portal/api/download/all"

logger = setup_logger("download_ccle")


def download_ccle(accession: str, output_dir: str, omic_type: str = "unknown") -> dict:
    """Download a CCLE dataset. Accession format: CCLE_{filename_stem}."""
    os.makedirs(output_dir, exist_ok=True)
    data_dir = os.path.join(output_dir, "data")
    meta_dir = os.path.join(output_dir, "metadata")
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(meta_dir, exist_ok=True)

    # Derive filename from accession
    filename_stem = accession.replace("CCLE_", "")
    logger.info(f"Downloading CCLE dataset: {accession}")

    # Get download URL from DepMap API
    resp = fetch_with_retry(DEPMAP_DOWNLOAD_API)
    entries = resp.json() if resp.status_code == 200 else []

    downloaded_files = []
    for entry in entries:
        entry_stem = entry.get("fileName", "").replace(".csv", "").replace(".", "_")
        if entry_stem == filename_stem or entry.get("fileName", "").replace(".csv", "") == filename_stem.replace("_", "."):
            url = entry.get("downloadUrl", "")
            if url:
                dest = os.path.join(data_dir, entry["fileName"])
                try:
                    download_file(url, dest)
                    size_mb = round(os.path.getsize(dest) / 1_000_000, 1)
                    downloaded_files.append({"name": entry["fileName"], "size_mb": size_mb})
                except Exception as e:
                    logger.warning(f"Download failed: {e}")

    # Save metadata
    meta = {"accession": accession, "source": "CCLE", "organism": "Homo sapiens", "disease": "cancer cell lines"}
    with open(os.path.join(meta_dir, "study_metadata.json"), "w") as f:
        json.dump(meta, f, indent=2)
    with open(os.path.join(meta_dir, "paper.json"), "w") as f:
        json.dump({"title": "CCLE", "authors": "Broad Institute", "journal": "DepMap", "doi": "", "date": "", "abstract": ""}, f, indent=2)

    return {"accession": accession, "output_dir": output_dir, "files_downloaded": len(downloaded_files), "data_files": downloaded_files}


def main():
    parser = argparse.ArgumentParser(description="Download CCLE dataset")
    parser.add_argument("--accession", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--omic-type", default="unknown")
    args = parser.parse_args()
    result = download_ccle(args.accession, args.output_dir, args.omic_type)
    json.dump(result, sys.stdout, indent=2)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Implement download_xena.py**

```python
#!/usr/bin/env python3
"""Download UCSC Xena datasets."""

import argparse
import json
import os
import sys

from scripts.utils import fetch_with_retry, download_file, setup_logger

logger = setup_logger("download_xena")


def build_xena_download_url(host: str, dataset_name: str) -> str:
    """Build the download URL for a Xena dataset."""
    return f"{host}/download/{dataset_name}"


def download_xena(accession: str, output_dir: str, omic_type: str = "unknown") -> dict:
    """Download a Xena dataset. Accession format: XENA_{dataset_name}."""
    os.makedirs(output_dir, exist_ok=True)
    data_dir = os.path.join(output_dir, "data")
    meta_dir = os.path.join(output_dir, "metadata")
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(meta_dir, exist_ok=True)

    dataset_name = accession.replace("XENA_", "").replace("_", ".")
    logger.info(f"Downloading Xena dataset: {dataset_name}")

    # Try each known hub
    hubs = ["https://tcga.xenahubs.net", "https://ucscpublic.xenahubs.net", "https://toil.xenahubs.net"]
    downloaded_files = []
    for hub in hubs:
        url = build_xena_download_url(hub, dataset_name)
        dest = os.path.join(data_dir, os.path.basename(dataset_name))
        try:
            download_file(url, dest)
            if os.path.exists(dest) and os.path.getsize(dest) > 0:
                size_mb = round(os.path.getsize(dest) / 1_000_000, 1)
                downloaded_files.append({"name": os.path.basename(dataset_name), "size_mb": size_mb})
                break
        except Exception:
            continue

    meta = {"accession": accession, "source": "Xena", "organism": "Homo sapiens"}
    with open(os.path.join(meta_dir, "study_metadata.json"), "w") as f:
        json.dump(meta, f, indent=2)
    with open(os.path.join(meta_dir, "paper.json"), "w") as f:
        json.dump({"title": "", "authors": "", "journal": "", "doi": "", "date": "", "abstract": ""}, f, indent=2)

    return {"accession": accession, "output_dir": output_dir, "files_downloaded": len(downloaded_files), "data_files": downloaded_files}


def main():
    parser = argparse.ArgumentParser(description="Download Xena dataset")
    parser.add_argument("--accession", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--omic-type", default="unknown")
    args = parser.parse_args()
    result = download_xena(args.accession, args.output_dir, args.omic_type)
    json.dump(result, sys.stdout, indent=2)


if __name__ == "__main__":
    main()
```

- [ ] **Step 5: Implement download_depmap.py**

```python
#!/usr/bin/env python3
"""Download DepMap datasets."""

import argparse
import json
import os
import sys

from scripts.utils import fetch_with_retry, download_file, setup_logger

DEPMAP_DOWNLOAD_API = "https://depmap.org/portal/api/download/all"

logger = setup_logger("download_depmap")


def download_depmap(accession: str, output_dir: str, omic_type: str = "unknown") -> dict:
    """Download a DepMap dataset. Accession format: DepMap_{filename_stem}."""
    os.makedirs(output_dir, exist_ok=True)
    data_dir = os.path.join(output_dir, "data")
    meta_dir = os.path.join(output_dir, "metadata")
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(meta_dir, exist_ok=True)

    filename_stem = accession.replace("DepMap_", "")
    logger.info(f"Downloading DepMap dataset: {accession}")

    resp = fetch_with_retry(DEPMAP_DOWNLOAD_API)
    entries = resp.json() if resp.status_code == 200 else []

    downloaded_files = []
    for entry in entries:
        entry_stem = entry.get("fileName", "").replace(".csv", "").replace(".", "_")
        if entry_stem == filename_stem:
            url = entry.get("downloadUrl", "")
            if url:
                dest = os.path.join(data_dir, entry["fileName"])
                try:
                    download_file(url, dest)
                    size_mb = round(os.path.getsize(dest) / 1_000_000, 1)
                    downloaded_files.append({"name": entry["fileName"], "size_mb": size_mb})
                except Exception as e:
                    logger.warning(f"Download failed: {e}")

    meta = {"accession": accession, "source": "DepMap", "organism": "Homo sapiens", "disease": "cancer cell lines"}
    with open(os.path.join(meta_dir, "study_metadata.json"), "w") as f:
        json.dump(meta, f, indent=2)
    with open(os.path.join(meta_dir, "paper.json"), "w") as f:
        json.dump({"title": "DepMap", "authors": "Broad Institute", "journal": "DepMap Portal", "doi": "", "date": "", "abstract": ""}, f, indent=2)

    return {"accession": accession, "output_dir": output_dir, "files_downloaded": len(downloaded_files), "data_files": downloaded_files}


def main():
    parser = argparse.ArgumentParser(description="Download DepMap dataset")
    parser.add_argument("--accession", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--omic-type", default="unknown")
    args = parser.parse_args()
    result = download_depmap(args.accession, args.output_dir, args.omic_type)
    json.dump(result, sys.stdout, indent=2)


if __name__ == "__main__":
    main()
```

- [ ] **Step 6: Implement download_scp.py**

```python
#!/usr/bin/env python3
"""Download Single Cell Portal study files."""

import argparse
import json
import os
import sys

from scripts.utils import fetch_with_retry, download_file, setup_logger

SCP_API_BASE = "https://singlecell.broadinstitute.org/single_cell/api/v1"

logger = setup_logger("download_scp")


def download_scp(accession: str, output_dir: str, omic_type: str = "unknown") -> dict:
    """Download an SCP study. Accession format: SCP####."""
    os.makedirs(output_dir, exist_ok=True)
    data_dir = os.path.join(output_dir, "data")
    meta_dir = os.path.join(output_dir, "metadata")
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(meta_dir, exist_ok=True)

    logger.info(f"Downloading SCP study: {accession}")

    token = os.environ.get("SCP_TOKEN", "")
    headers = {"Authorization": f"Bearer {token}"} if token else {}

    # Get study info
    resp = fetch_with_retry(f"{SCP_API_BASE}/studies/{accession}", headers=headers)
    if resp.status_code != 200:
        logger.warning(f"SCP API returned {resp.status_code} for {accession}")
        return {"accession": accession, "output_dir": output_dir, "files_downloaded": 0, "data_files": []}

    study = resp.json()

    # Save metadata
    meta = {"accession": accession, "source": "SCP", "name": study.get("name", ""), "description": study.get("description", "")}
    with open(os.path.join(meta_dir, "study_metadata.json"), "w") as f:
        json.dump(meta, f, indent=2)
    with open(os.path.join(meta_dir, "paper.json"), "w") as f:
        json.dump({"title": "", "authors": "", "journal": "", "doi": "", "date": "", "abstract": study.get("description", "")}, f, indent=2)

    # Download study files
    downloaded_files = []
    for sf in study.get("study_files", []):
        url = sf.get("download_url", "")
        name = sf.get("name", "unknown_file")
        if url:
            dest = os.path.join(data_dir, name)
            try:
                download_file(url, dest)
                size_mb = round(os.path.getsize(dest) / 1_000_000, 1)
                downloaded_files.append({"name": name, "size_mb": size_mb})
            except Exception as e:
                logger.warning(f"Failed to download {name}: {e}")

    # Generate sample metadata stub
    with open(os.path.join(meta_dir, "sample_metadata.csv"), "w") as f:
        f.write("sample_id,condition,batch\n")

    return {"accession": accession, "output_dir": output_dir, "files_downloaded": len(downloaded_files), "data_files": downloaded_files}


def main():
    parser = argparse.ArgumentParser(description="Download SCP study")
    parser.add_argument("--accession", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--omic-type", default="unknown")
    args = parser.parse_args()
    result = download_scp(args.accession, args.output_dir, args.omic_type)
    json.dump(result, sys.stdout, indent=2)


if __name__ == "__main__":
    main()
```

- [ ] **Step 7: Run all download tests**

```bash
conda run -n dataseek pytest tests/test_download_ccle.py tests/test_download_xena.py tests/test_download_depmap.py tests/test_download_scp.py -v
```

- [ ] **Step 8: Commit**

```bash
git add scripts/download_ccle.py scripts/download_xena.py scripts/download_depmap.py scripts/download_scp.py tests/test_download_ccle.py tests/test_download_xena.py tests/test_download_depmap.py tests/test_download_scp.py
git commit -m "feat: download scripts for CCLE, Xena, DepMap, and SCP"
```

---

## Task 10: Dataseek Search Skill

**Files:**
- Create: `.claude/skills/dataseek-search/SKILL.md`

- [ ] **Step 1: Create skill directory**

```bash
mkdir -p .claude/skills/dataseek-search
```

- [ ] **Step 2: Write the search skill**

Create `.claude/skills/dataseek-search/SKILL.md` with the full skill definition. This skill:

1. Parses `--omic` (required) and optional filters from user args
2. Applies source routing table to determine which sources to search
3. Dispatches source agents in parallel using the Agent tool (one message, multiple Agent calls)
4. Each source agent runs: `conda run -n dataseek python scripts/search_{source}.py --omic {omic} ...` and returns the JSON output
5. Coordinator merges all results, deduplicates by accession/title, ranks by sample_count + metadata_quality + recency
6. Saves results to `results/search_cache/{date}_{omic}_{hash}.json` via `scripts/utils.py`
7. Presents tiered output: top 10 with full cards, remaining as one-line table

The skill file should contain:
- The source routing table (copied from spec)
- The SearchResult schema
- Exact agent dispatch prompts for each source
- The merge/deduplicate/rank algorithm
- The output formatting template for both tiers
- `--min-samples` filtering: applied as post-filter by the coordinator after merging results (filter out results where `sample_count < min_samples`)
- `--date-from`/`--date-to`: only passed to GEO search script (other sources lack date-based API query support); coordinator notes this in output if date filters were specified

- [ ] **Step 3: Commit**

```bash
git add .claude/skills/dataseek-search/SKILL.md
git commit -m "feat: dataseek-search skill for parallel multi-source omics search"
```

---

## Task 11: Dataseek Download Skill

**Files:**
- Create: `.claude/skills/dataseek-download/SKILL.md`

- [ ] **Step 1: Create skill directory**

```bash
mkdir -p .claude/skills/dataseek-download
```

- [ ] **Step 2: Write the download skill**

Create `.claude/skills/dataseek-download/SKILL.md` with the full skill definition. This skill:

1. Takes one or more accession IDs as args
2. For each accession, looks up source from latest search cache files in `results/search_cache/`
3. If cache miss, determines source from accession prefix (GSE→GEO, CCLE_→CCLE, XENA_→Xena, DepMap_→DepMap, SCP→SCP)
4. For each dataset, runs: `conda run -n dataseek python scripts/download_{source}.py --accession {accession} --output-dir downloads/{accession} --omic-type {omic}`
5. After download, generates summary report to `results/reports/{accession}_summary.md`
6. Generates pipeline config in `downloads/{accession}/pipeline_config/` based on omic type
7. Reports download status: files downloaded, total size, report location

The skill should contain:
- Source inference from accession prefix
- Pipeline config generation rules per omic type
- The summary report template
- Error handling (retry, skip, report)

- [ ] **Step 3: Commit**

```bash
git add .claude/skills/dataseek-download/SKILL.md
git commit -m "feat: dataseek-download skill for dataset retrieval and config generation"
```

---

## Task 12: CLAUDE.md Agent Orchestration Spec

**Files:**
- Create: `CLAUDE.md`

- [ ] **Step 1: Write CLAUDE.md**

The CLAUDE.md should define:

1. **Project overview** — what dataseek does
2. **Directory structure** — reference to spec
3. **Conda environment** — `dataseek`, how to activate, key packages
4. **Agent definitions table** — all 7 agents with:
   - Agent name
   - Role/trigger
   - Input
   - Output
   - Required tools/scripts
5. **Source routing table** (from spec)
6. **SearchResult schema** (from spec)
7. **Slash commands** — `/dataseek-search` and `/dataseek-download` with full parameter docs
8. **Script invocation pattern** — `conda run -n dataseek python scripts/{script}.py --args`
9. **Error handling rules** (from spec)
10. **Data flow diagram** showing search → cache → download → pipeline_config

- [ ] **Step 2: Commit**

```bash
git add CLAUDE.md
git commit -m "feat: CLAUDE.md orchestration spec for dataseek agent team"
```

---

## Task 13: Integration Test

**Files:**
- Create: `tests/test_integration.py`

- [ ] **Step 1: Write integration test**

```python
# tests/test_integration.py
"""Integration tests: verify scripts can be invoked via CLI and return valid JSON."""

import json
import subprocess
import sys

import pytest


SCRIPTS = [
    ("scripts/search_geo.py", ["--omic", "scRNAseq", "--organism", "human", "--max-results", "1"]),
    ("scripts/search_ccle.py", ["--omic", "bulkRNAseq", "--max-results", "1"]),
    ("scripts/search_xena.py", ["--omic", "bulkRNAseq", "--max-results", "1"]),
    ("scripts/search_depmap.py", ["--omic", "CRISPR", "--max-results", "1"]),
    ("scripts/search_scp.py", ["--omic", "scRNAseq", "--max-results", "1"]),
]


@pytest.mark.integration
class TestScriptCLI:
    @pytest.mark.parametrize("script,args", SCRIPTS)
    def test_script_returns_valid_json(self, script, args):
        """Each search script should return valid JSON when invoked."""
        result = subprocess.run(
            ["conda", "run", "-n", "dataseek", "python", script] + args,
            capture_output=True, text=True, timeout=60,
        )
        assert result.returncode == 0, f"Script failed: {result.stderr}"
        data = json.loads(result.stdout)
        assert isinstance(data, list)
```

- [ ] **Step 2: Run integration tests** (requires network access)

```bash
conda run -n dataseek pytest tests/test_integration.py -v -m integration
```

- [ ] **Step 3: Commit**

```bash
git add tests/test_integration.py
git commit -m "test: integration tests for all search script CLIs"
```

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
    verify: bool | str = True,
) -> requests.Response:
    """GET request with exponential backoff retry on 5xx or connection errors."""
    import urllib3
    import warnings
    for attempt in range(max_retries):
        try:
            resp = requests.get(url, params=params, headers=headers, timeout=timeout, verify=verify)
            if resp.status_code < 500:
                return resp
        except requests.exceptions.SSLError:
            # Retry once with SSL verification disabled on SSL errors
            if attempt == 0:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", urllib3.exceptions.InsecureRequestWarning)
                    try:
                        resp = requests.get(url, params=params, headers=headers, timeout=timeout, verify=False)
                        if resp.status_code < 500:
                            return resp
                    except requests.RequestException:
                        pass
            if attempt == max_retries - 1:
                raise
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

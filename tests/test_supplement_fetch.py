"""Tests for supplement_fetch.py — PMC supplementary table extraction."""

import csv
import os
import tempfile

import pandas as pd
import pytest
import responses

from scripts.supplement_fetch import (
    IDCONV_URL,
    PMC_OA_URL,
    fetch_supplementary_tables,
    match_supplement_to_samples,
    parse_tabular_file,
    resolve_pmcid,
)


# ---------------------------------------------------------------------------
# TestResolvePmcid
# ---------------------------------------------------------------------------

class TestResolvePmcid:
    @responses.activate
    def test_resolves_pmid_to_pmcid(self):
        """Mock IDCONV API returning PMC10500001 for a PMID lookup."""
        xml_body = (
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<pmcids status="ok">'
            '  <record requested-id="38000001" pmcid="PMC10500001" pmid="38000001" doi="10.1038/example"/>'
            '</pmcids>'
        )
        responses.add(responses.GET, IDCONV_URL, body=xml_body, status=200)
        result = resolve_pmcid(pmid="38000001")
        assert result == "PMC10500001"

    @responses.activate
    def test_returns_none_on_no_match(self):
        """Mock IDCONV API returning a record with no pmcid attribute."""
        xml_body = (
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<pmcids status="ok">'
            '  <record requested-id="99999999" pmid="99999999"/>'
            '</pmcids>'
        )
        responses.add(responses.GET, IDCONV_URL, body=xml_body, status=200)
        result = resolve_pmcid(pmid="99999999")
        assert result is None


# ---------------------------------------------------------------------------
# TestParseTabularFile
# ---------------------------------------------------------------------------

class TestParseTabularFile:
    def test_parses_csv(self, tmp_path):
        """Write a temp CSV and verify parse_tabular_file returns a DataFrame."""
        csv_file = tmp_path / "samples.csv"
        csv_file.write_text("sample_id,tissue,age\nS1,lung,45\nS2,brain,60\n")
        df = parse_tabular_file(str(csv_file))
        assert df is not None
        assert isinstance(df, pd.DataFrame)
        assert list(df.columns) == ["sample_id", "tissue", "age"]
        assert len(df) == 2

    def test_parses_tsv(self, tmp_path):
        """Write a temp TSV and verify parse_tabular_file returns a DataFrame."""
        tsv_file = tmp_path / "samples.tsv"
        tsv_file.write_text("patient_id\tsex\tstage\nP1\tM\tII\nP2\tF\tIII\n")
        df = parse_tabular_file(str(tsv_file))
        assert df is not None
        assert isinstance(df, pd.DataFrame)
        assert "patient_id" in df.columns
        assert len(df) == 2

    def test_returns_none_for_invalid_file(self, tmp_path):
        """A file that doesn't exist or has an unsupported extension returns None."""
        result = parse_tabular_file(str(tmp_path / "nonexistent.csv"))
        assert result is None


# ---------------------------------------------------------------------------
# TestMatchSupplementToSamples
# ---------------------------------------------------------------------------

class TestMatchSupplementToSamples:
    def _make_df(self):
        """Build a DataFrame with clinical-looking columns."""
        return pd.DataFrame({
            "Patient ID": ["PT001", "PT002", "PT003"],
            "Tissue Site": ["lung", "brain", "liver"],
            "Treatment": ["chemo", "none", "chemo"],
            "Stage": ["II", "III", "I"],
        })

    def test_matches_sample_ids_and_maps_columns(self):
        """DataFrame with Patient ID, Tissue Site, Treatment, Stage — all should map."""
        df = self._make_df()
        result = match_supplement_to_samples(df, sample_ids=["PT001", "PT002"])
        assert result is not None
        assert isinstance(result, list)
        # Should only return rows matching the provided sample IDs
        returned_ids = {row.get("sample_id") for row in result}
        assert "PT001" in returned_ids
        assert "PT002" in returned_ids
        assert "PT003" not in returned_ids
        # Schema columns should appear in at least one row
        schema_keys_present = set()
        for row in result:
            schema_keys_present.update(row.keys())
        assert any(k in schema_keys_present for k in ("tissue_site", "treatment_status", "stage"))

    def test_returns_none_when_no_sample_id_column(self):
        """DataFrame without any sample-ID-like column returns None."""
        df = pd.DataFrame({
            "alpha": [1, 2],
            "beta": [3, 4],
        })
        result = match_supplement_to_samples(df, sample_ids=["PT001"])
        assert result is None

    def test_returns_none_when_no_clinical_columns(self):
        """DataFrame has a sample ID column but zero clinical columns returns None."""
        df = pd.DataFrame({
            "sample_id": ["S1", "S2"],
            "random_col_xyz": ["foo", "bar"],
        })
        result = match_supplement_to_samples(df, sample_ids=["S1"])
        assert result is None


# ---------------------------------------------------------------------------
# TestFetchSupplementaryTables
# ---------------------------------------------------------------------------

class TestFetchSupplementaryTables:
    def test_returns_none_when_no_doi_or_pmid(self):
        """Calling with neither doi nor pmid returns None immediately."""
        result = fetch_supplementary_tables()
        assert result is None

    @responses.activate
    def test_returns_none_when_pmcid_not_found(self):
        """When IDCONV returns no PMCID, the pipeline returns None."""
        xml_body = (
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<pmcids status="ok">'
            '  <record requested-id="00000000" pmid="00000000"/>'
            '</pmcids>'
        )
        responses.add(responses.GET, IDCONV_URL, body=xml_body, status=200)
        result = fetch_supplementary_tables(pmid="00000000")
        assert result is None

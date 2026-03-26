# tests/test_sample_utils.py
import csv
import os

import pytest

from scripts.sample_utils import (
    STANDARD_COLUMNS,
    SOURCE_EXTRA_COLUMNS,
    SAMPLE_CAP,
    normalize_sample_table,
    write_sample_csv,
    summarize_samples,
    build_sample_table,
    _build_summary_dict,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_rows(n, extra=None):
    """Generate n minimal sample dicts, optionally merging extra keys."""
    rows = []
    for i in range(n):
        row = {
            "sample_id": f"S{i:04d}",
            "tissue_site": "brain",
            "sample_type": "primary",
            "treatment_status": "untreated",
            "disease": "glioblastoma",
        }
        if extra:
            row.update(extra)
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# TestNormalizeSampleTable
# ---------------------------------------------------------------------------

class TestNormalizeSampleTable:

    def test_fills_missing_standard_columns_with_na(self):
        rows = [{"sample_id": "S0001"}]
        result = normalize_sample_table(rows, "GEO")
        assert result[0]["tissue_site"] == "N/A"
        assert result[0]["disease"] == "N/A"
        assert result[0]["cell_type"] == "N/A"

    def test_preserves_existing_values(self):
        rows = [{"sample_id": "S0001", "tissue_site": "brain", "disease": "GBM"}]
        result = normalize_sample_table(rows, "GEO")
        assert result[0]["tissue_site"] == "brain"
        assert result[0]["disease"] == "GBM"

    def test_appends_geo_source_extras(self):
        rows = [{"sample_id": "S0001", "platform_id": "GPL1234"}]
        result = normalize_sample_table(rows, "GEO")
        assert "platform_id" in result[0]
        assert result[0]["platform_id"] == "GPL1234"

    def test_appends_ccle_source_extras(self):
        rows = [{"sample_id": "S0001", "lineage": "breast", "sublineage": "luminal", "cosmic_id": "999"}]
        result = normalize_sample_table(rows, "CCLE")
        assert result[0]["lineage"] == "breast"
        assert result[0]["sublineage"] == "luminal"
        assert result[0]["cosmic_id"] == "999"

    def test_fills_missing_source_extras_with_na(self):
        rows = [{"sample_id": "S0001"}]
        result = normalize_sample_table(rows, "CCLE")
        assert result[0]["lineage"] == "N/A"
        assert result[0]["sublineage"] == "N/A"
        assert result[0]["cosmic_id"] == "N/A"

    def test_fills_missing_source_extras_geo_with_na(self):
        rows = [{"sample_id": "S0001"}]
        result = normalize_sample_table(rows, "GEO")
        assert result[0]["platform_id"] == "N/A"

    def test_caps_at_sample_cap(self):
        rows = _make_rows(SAMPLE_CAP + 50)
        result = normalize_sample_table(rows, "GEO")
        assert len(result) == SAMPLE_CAP

    def test_does_not_cap_when_under_limit(self):
        rows = _make_rows(10)
        result = normalize_sample_table(rows, "GEO")
        assert len(result) == 10

    def test_all_standard_columns_present(self):
        rows = [{"sample_id": "S0001"}]
        result = normalize_sample_table(rows, "GEO")
        for col in STANDARD_COLUMNS:
            assert col in result[0], f"Missing standard column: {col}"

    def test_unknown_source_has_no_extras(self):
        rows = [{"sample_id": "S0001"}]
        result = normalize_sample_table(rows, "UNKNOWN")
        for col in STANDARD_COLUMNS:
            assert col in result[0]
        # No extra columns beyond standard
        assert set(result[0].keys()) == set(STANDARD_COLUMNS)

    def test_xena_source_extras(self):
        rows = [{"sample_id": "S0001", "cohort_name": "TCGA-GBM"}]
        result = normalize_sample_table(rows, "Xena")
        assert result[0]["cohort_name"] == "TCGA-GBM"

    def test_scp_source_extras(self):
        rows = [{"sample_id": "S0001", "cell_count": "5000", "library_prep": "10x"}]
        result = normalize_sample_table(rows, "SCP")
        assert result[0]["cell_count"] == "5000"
        assert result[0]["library_prep"] == "10x"


# ---------------------------------------------------------------------------
# TestWriteSampleCsv
# ---------------------------------------------------------------------------

class TestWriteSampleCsv:

    def _make_table_data(self, n=3, capped=False, extra_rows=0):
        """Build a minimal table_data dict."""
        rows = _make_rows(n)
        source = "GEO"
        normalized = normalize_sample_table(rows, source)
        cols = STANDARD_COLUMNS + SOURCE_EXTRA_COLUMNS[source]
        csv_rows = [[r[c] for c in cols] for r in normalized]
        total = n + extra_rows if capped else n
        return {
            "total_samples": total,
            "shown_samples": n,
            "capped": capped,
            "columns": cols,
            "rows": csv_rows,
            "summary": {"primary": n, "normal": 0},
        }

    def test_writes_csv_file(self, tmp_path):
        table_data = self._make_table_data(3)
        filepath = write_sample_csv("GSE999", table_data, str(tmp_path))
        assert os.path.exists(filepath)

    def test_csv_has_header_row(self, tmp_path):
        table_data = self._make_table_data(3)
        filepath = write_sample_csv("GSE999", table_data, str(tmp_path))
        with open(filepath, newline="") as f:
            reader = csv.reader(f)
            header = next(reader)
        assert header == table_data["columns"]

    def test_csv_has_correct_data_rows(self, tmp_path):
        table_data = self._make_table_data(3)
        filepath = write_sample_csv("GSE999", table_data, str(tmp_path))
        with open(filepath, newline="") as f:
            reader = csv.reader(f)
            next(reader)  # skip header
            data_rows = [row for row in reader if not row[0].startswith("#")]
        assert len(data_rows) == 3

    def test_csv_returns_filepath(self, tmp_path):
        table_data = self._make_table_data(3)
        filepath = write_sample_csv("GSE999", table_data, str(tmp_path))
        assert isinstance(filepath, str)
        assert "GSE999" in filepath or str(tmp_path) in filepath

    def test_writes_cap_comment_when_capped(self, tmp_path):
        table_data = self._make_table_data(n=3, capped=True, extra_rows=100)
        filepath = write_sample_csv("GSE999", table_data, str(tmp_path))
        with open(filepath) as f:
            content = f.read()
        assert "additional samples not shown" in content
        assert "100" in content

    def test_no_cap_comment_when_not_capped(self, tmp_path):
        table_data = self._make_table_data(3)
        filepath = write_sample_csv("GSE999", table_data, str(tmp_path))
        with open(filepath) as f:
            content = f.read()
        assert "additional samples not shown" not in content

    def test_creates_output_dir_if_missing(self, tmp_path):
        table_data = self._make_table_data(2)
        nested_dir = str(tmp_path / "nested" / "output")
        filepath = write_sample_csv("GSE999", table_data, nested_dir)
        assert os.path.exists(filepath)


# ---------------------------------------------------------------------------
# TestSummarizeSamples
# ---------------------------------------------------------------------------

class TestSummarizeSamples:

    def _make_table_data_from_rows(self, rows, source="GEO"):
        normalized = normalize_sample_table(rows, source)
        cols = STANDARD_COLUMNS + SOURCE_EXTRA_COLUMNS.get(source, [])
        csv_rows = [[r[c] for c in cols] for r in normalized]
        summary_dict = _build_summary_dict(normalized)
        return {
            "total_samples": len(rows),
            "shown_samples": len(normalized),
            "capped": False,
            "columns": cols,
            "rows": csv_rows,
            "summary": summary_dict,
        }

    def test_summary_contains_sample_count(self):
        rows = _make_rows(12)
        table_data = self._make_table_data_from_rows(rows)
        summary = summarize_samples(table_data)
        assert "12" in summary

    def test_summary_contains_sample_type_counts(self):
        rows = _make_rows(8)
        for i in range(4):
            rows[i]["sample_type"] = "normal"
        table_data = self._make_table_data_from_rows(rows)
        summary = summarize_samples(table_data)
        assert "primary" in summary
        assert "normal" in summary

    def test_summary_contains_treatment_counts(self):
        rows = _make_rows(9)
        for i in range(3):
            rows[i]["treatment_status"] = "treated"
        table_data = self._make_table_data_from_rows(rows)
        summary = summarize_samples(table_data)
        assert "treated" in summary
        assert "untreated" in summary

    def test_summary_handles_all_na(self):
        rows = [{"sample_id": f"S{i:04d}"} for i in range(5)]
        table_data = self._make_table_data_from_rows(rows)
        summary = summarize_samples(table_data)
        # Should not crash; should contain sample count
        assert "5" in summary

    def test_summary_is_string(self):
        rows = _make_rows(3)
        table_data = self._make_table_data_from_rows(rows)
        summary = summarize_samples(table_data)
        assert isinstance(summary, str)

    def test_summary_format_matches_spec(self):
        # "Samples (12): 8 primary, 4 normal | untreated: 9, treated: 3"
        rows = _make_rows(12)
        for i in range(4):
            rows[i]["sample_type"] = "normal"
        for i in range(3):
            rows[i]["treatment_status"] = "treated"
        table_data = self._make_table_data_from_rows(rows)
        summary = summarize_samples(table_data)
        assert summary.startswith("Samples (")
        assert "|" in summary


# ---------------------------------------------------------------------------
# TestBuildSampleTable
# ---------------------------------------------------------------------------

class TestBuildSampleTable:

    def test_returns_dict_with_required_keys(self):
        rows = _make_rows(5)
        result = build_sample_table(rows, "GEO")
        required_keys = {"total_samples", "shown_samples", "capped", "columns", "rows", "summary"}
        assert required_keys <= set(result.keys())

    def test_total_samples_matches_row_count_by_default(self):
        rows = _make_rows(5)
        result = build_sample_table(rows, "GEO")
        assert result["total_samples"] == 5
        assert result["shown_samples"] == 5

    def test_total_samples_uses_total_count_when_provided(self):
        rows = _make_rows(5)
        result = build_sample_table(rows, "GEO", total_count=200)
        assert result["total_samples"] == 200
        assert result["shown_samples"] == 5

    def test_capped_false_when_under_cap(self):
        rows = _make_rows(10)
        result = build_sample_table(rows, "GEO")
        assert result["capped"] is False

    def test_capped_true_when_over_cap(self):
        rows = _make_rows(SAMPLE_CAP + 10)
        result = build_sample_table(rows, "GEO", total_count=SAMPLE_CAP + 10)
        assert result["capped"] is True

    def test_columns_include_standard_and_source_extras(self):
        rows = _make_rows(3)
        result = build_sample_table(rows, "GEO")
        for col in STANDARD_COLUMNS:
            assert col in result["columns"]
        for col in SOURCE_EXTRA_COLUMNS["GEO"]:
            assert col in result["columns"]

    def test_rows_are_lists(self):
        rows = _make_rows(3)
        result = build_sample_table(rows, "GEO")
        assert isinstance(result["rows"], list)
        assert all(isinstance(r, list) for r in result["rows"])

    def test_rows_length_matches_shown_samples(self):
        rows = _make_rows(5)
        result = build_sample_table(rows, "GEO")
        assert len(result["rows"]) == result["shown_samples"]

    def test_rows_values_align_with_columns(self):
        rows = [{"sample_id": "S0001", "tissue_site": "lung"}]
        result = build_sample_table(rows, "GEO")
        col_idx = result["columns"].index("sample_id")
        assert result["rows"][0][col_idx] == "S0001"
        col_idx = result["columns"].index("tissue_site")
        assert result["rows"][0][col_idx] == "lung"

    def test_summary_is_dict_in_table(self):
        rows = _make_rows(3)
        result = build_sample_table(rows, "GEO")
        assert isinstance(result["summary"], dict)

    def test_capped_at_sample_cap_rows(self):
        rows = _make_rows(SAMPLE_CAP + 20)
        result = build_sample_table(rows, "GEO", total_count=SAMPLE_CAP + 20)
        assert result["shown_samples"] == SAMPLE_CAP
        assert len(result["rows"]) == SAMPLE_CAP


# ---------------------------------------------------------------------------
# TestBuildSummaryDict
# ---------------------------------------------------------------------------

class TestBuildSummaryDict:

    def test_counts_sample_types(self):
        rows = [
            {"sample_type": "primary", "treatment_status": "untreated"},
            {"sample_type": "primary", "treatment_status": "untreated"},
            {"sample_type": "normal", "treatment_status": "untreated"},
        ]
        result = _build_summary_dict(rows)
        assert result.get("primary") == 2
        assert result.get("normal") == 1

    def test_counts_treatment_status(self):
        rows = [
            {"sample_type": "primary", "treatment_status": "treated"},
            {"sample_type": "primary", "treatment_status": "untreated"},
            {"sample_type": "primary", "treatment_status": "untreated"},
        ]
        result = _build_summary_dict(rows)
        assert result.get("treated") == 1
        assert result.get("untreated") == 2

    def test_ignores_na_values(self):
        rows = [
            {"sample_type": "N/A", "treatment_status": "N/A"},
            {"sample_type": "primary", "treatment_status": "untreated"},
        ]
        result = _build_summary_dict(rows)
        assert "N/A" not in result

    def test_empty_rows_returns_empty_dict(self):
        result = _build_summary_dict([])
        assert result == {}

"""Integration tests: verify scripts can be invoked via CLI and return valid JSON."""

import json
import subprocess

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

"""Unit tests for utility functions in app/utils.py."""
import pandas as pd
import pytest
from app.utils import getOverLapData, filterPeaksFile


class TestGetOverLapData:
    """getOverLapData should sort samples alphabetically and exclude Control."""

    def _make_samples(self, names):
        # Each sample has one replicate file ending in .txt
        return {name: {f"{name}_rep1.txt": None} for name in names}

    def test_sorts_alphabetically(self):
        samples = self._make_samples(["Zebra", "Apple", "Mango"])
        result = getOverLapData(samples)
        assert list(result.keys()) == ["Apple", "Mango", "Zebra"]

    def test_excludes_control(self):
        samples = self._make_samples(["SampleA", "Control", "SampleB"])
        result = getOverLapData(samples)
        assert "Control" not in result
        assert set(result.keys()) == {"SampleA", "SampleB"}

    def test_strips_extension_from_replicates(self):
        samples = {"Sample1": {"rep1.txt": None, "rep2.txt": None}}
        result = getOverLapData(samples)
        assert result["Sample1"] == ["rep1", "rep2"]

    def test_empty_input(self):
        assert getOverLapData({}) == {}

    def test_only_control_returns_empty(self):
        samples = self._make_samples(["Control"])
        assert getOverLapData(samples) == {}


class TestFilterPeaksFile:
    """filterPeaksFile should filter by length, remove PTMs, handle missing column."""

    def _df(self, peptides):
        return pd.DataFrame({"Peptide": peptides})

    def test_filters_by_length(self):
        samples = {"f.csv": self._df(["ACDE", "ACDEFGHIJK"])}
        result = filterPeaksFile(samples, minLen=5, maxLen=9)
        assert result["f.csv"].shape[0] == 0  # 4-mer excluded, 10-mer excluded

    def test_accepts_valid_peptides(self):
        samples = {"f.csv": self._df(["ACDEFGHI"])}
        result = filterPeaksFile(samples, minLen=1, maxLen=20)
        assert result["f.csv"].shape[0] == 1

    def test_raises_on_missing_peptide_column(self):
        df = pd.DataFrame({"NotPeptide": ["ACDE"]})
        with pytest.raises(KeyError):
            filterPeaksFile({"f.csv": df})

    def test_length_filter_returns_empty_df(self):
        samples = {"f.csv": self._df(["ACDE"])}
        # Length 4; filter to 5-9 → length step returns empty df (no exception raised)
        result = filterPeaksFile(samples, minLen=5, maxLen=9)
        assert result["f.csv"].shape[0] == 0

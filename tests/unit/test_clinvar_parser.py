"""Unit tests for ClinVar parser."""

import pytest
from refold.data.clinvar_parser import _parse_hgvs_protein


@pytest.mark.parametrize("hgvs,expected", [
    ("p.Gly12Val", ("G", 12, "V")),
    ("p.Arg506Gln", ("R", 506, "Q")),
    ("p.Tyr152Cys", ("Y", 152, "C")),
    ("p.Arg175His", ("R", 175, "H")),
])
def test_parse_hgvs_protein_valid(hgvs, expected):
    result = _parse_hgvs_protein(hgvs)
    assert result == expected


@pytest.mark.parametrize("hgvs", [
    "c.35G>A",
    "p.Gly12Ter",
    "p.Ter100Met",
    "not_hgvs",
    "",
    "p.Ala12Ala",
])
def test_parse_hgvs_protein_invalid(hgvs):
    result = _parse_hgvs_protein(hgvs)
    assert result is None

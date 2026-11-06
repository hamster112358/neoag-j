"""Regression tests for row-level evidence provenance (B-E review C1-01, C2-02).

C1-01: PRIME gap-filling fabricated pseudo-random scores for pairs the tool
never returned and wrote them into the real tool file with no marking; the
non-strict missing-tool fallback could also register a fully stubbed tool
output as real.
C2-02: immunogenicity stub detection only looked at the file-level
evidence_status column, so row-level placeholders were invisible.
"""

from pathlib import Path

from neoag.adapters.prime import parse_prime, write_prime_evidence
from neoag.evidence_provenance import (
    MODE_TOOL_RUN,
    SOURCE_LOCAL,
    STATUS_INVALID_FOR_PRODUCTION,
    STATUS_REAL,
    STATUS_STUB,
    ProvenanceRecord,
    attach_provenance,
)
from neoag.immunogenicity_composite import _is_stub_evidence, component_scores
from neoag.tools.runner import _prime_evidence_rows, _write_prime_stub
from neoag.utils import read_tsv

PROFILE = {"immunogenicity": {"enabled": True, "sources": ["prime"]}}


def _real_provenance() -> ProvenanceRecord:
    return ProvenanceRecord(tool="prime", source=SOURCE_LOCAL, status=STATUS_REAL, mode=MODE_TOOL_RUN, file="tool.tsv")


def test_missing_pairs_are_blank_placeholder_rows_not_fabricated_scores():
    by_pep = {"AATKWT": {"Peptide": "AATKWT", "Score_A0201": "0.8", "%Rank_A0201": "1.2"}}
    rows = _prime_evidence_rows(
        [("AATKWT", "HLA-A*02:01"), ("GGGXXX", "HLA-A*02:01")],
        by_pep,
        "S1",
        "tools/prime.tsv",
    )
    real, gap = rows
    assert real["prime_score"] == "0.8" and real["prime_rank"] == "1.2"
    assert real["evidence_row_status"] == STATUS_REAL
    assert gap["prime_score"] == "" and gap["prime_rank"] == ""
    assert gap["evidence_row_status"] == STATUS_STUB


def test_all_stub_file_overrides_a_real_provenance_argument(tmp_path: Path):
    tool_stub = tmp_path / "prime_tool.tsv"
    _write_prime_stub([("AAAAAAA", "HLA-A*02:01")], tool_stub, "S1")
    rows = parse_prime(tool_stub, "S1")
    assert all(row["source_file"] == "stub" for row in rows)  # marker survives parsing

    out = tmp_path / "prime_evidence.tsv"
    write_prime_evidence(out, rows, provenance=_real_provenance())
    written = read_tsv(out)
    assert all(row["evidence_row_status"] == STATUS_STUB for row in written)
    assert all(row["evidence_status"] == STATUS_INVALID_FOR_PRODUCTION for row in written)
    assert all(_is_stub_evidence(row) for row in written)


def test_sticky_row_marks_survive_a_real_file(tmp_path: Path):
    rows = [
        {
            "sample_id": "S1", "peptide": "AATKWT", "hla_allele": "HLA-A*02:01",
            "prime_score": "0.9", "prime_rank": "0.5", "source_file": "tools/prime.tsv",
            "evidence_row_status": STATUS_REAL,
        },
        {
            "sample_id": "S1", "peptide": "GGGXXX", "hla_allele": "HLA-A*02:01",
            "prime_score": "", "prime_rank": "", "source_file": "tools/prime.tsv",
            "evidence_row_status": STATUS_STUB,
        },
    ]
    out = tmp_path / "prime_evidence.tsv"
    write_prime_evidence(out, rows, provenance=_real_provenance())
    written = read_tsv(out)
    assert [row["evidence_row_status"] for row in written] == [STATUS_REAL, STATUS_STUB]
    # The file itself is real tool output; only the placeholder row is stub.
    assert all(row["evidence_status"] == STATUS_REAL for row in written)
    assert not _is_stub_evidence(written[0])
    assert _is_stub_evidence(written[1])


def test_is_stub_evidence_row_level_precedence_and_file_fallback():
    assert _is_stub_evidence({"evidence_row_status": STATUS_STUB, "evidence_status": STATUS_REAL})
    assert not _is_stub_evidence({"evidence_row_status": STATUS_REAL, "evidence_status": STATUS_STUB})
    assert _is_stub_evidence({"evidence_status": "stub"})
    assert _is_stub_evidence({"evidence_status": STATUS_INVALID_FOR_PRODUCTION})
    assert not _is_stub_evidence({"evidence_status": STATUS_REAL})
    assert not _is_stub_evidence({})


def test_component_scores_ignore_row_level_stub_rows():
    row = {
        "prime_score": "0.95",
        "prime_rank": "0.1",
        "evidence_row_status": STATUS_STUB,
        "evidence_status": STATUS_REAL,
    }
    assert component_scores(row, PROFILE) == {}
    row_real = dict(row, evidence_row_status=STATUS_REAL)
    assert "prime" in component_scores(row_real, PROFILE)


def test_attach_provenance_marks_file_level_stub_rows_even_when_status_is_invalid():
    rec = ProvenanceRecord(tool="prime", status=STATUS_INVALID_FOR_PRODUCTION, mode="stub", file="stub")
    attached = attach_provenance([{"peptide": "X", "source_file": "stub"}], rec)
    assert attached[0]["evidence_row_status"] == STATUS_STUB
    rec_real = ProvenanceRecord(tool="prime", status=STATUS_REAL, mode=MODE_TOOL_RUN, file="tool.tsv")
    attached_real = attach_provenance([{"peptide": "X", "source_file": "tool.tsv"}], rec_real)
    assert attached_real[0]["evidence_row_status"] == STATUS_REAL

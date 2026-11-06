from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace


def _module():
    path = Path.cwd() / "scripts/rna_allele_counts_pysam.py"
    spec = importlib.util.spec_from_file_location("rna_allele_counts_pysam", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_classifies_snv_and_indel_observations():
    module = _module()
    call = module.classify_observation
    assert call(query_sequence="AT", query_position=0, indel=0, ref="A", alts=["T"]) == -1
    assert call(query_sequence="TT", query_position=0, indel=0, ref="A", alts=["T"]) == 0
    assert call(query_sequence="ATG", query_position=0, indel=2, ref="A", alts=["ATG"]) == 0
    assert call(query_sequence="A", query_position=0, indel=-2, ref="ATG", alts=["A"]) == 0
    assert call(query_sequence="A", query_position=0, indel=0, ref="ATG", alts=["A"]) == -1


def test_multiallelic_reads_classify_to_their_own_alt():
    # Review C3-02: every ALT allele must be countable, not only ALT[0].
    module = _module()
    call = module.classify_observation
    assert call(query_sequence="TT", query_position=0, indel=0, ref="A", alts=["C", "T"]) == 1
    assert call(query_sequence="CC", query_position=0, indel=0, ref="A", alts=["C", "T"]) == 0
    assert call(query_sequence="GG", query_position=0, indel=0, ref="A", alts=["C", "T"]) is None


def test_complex_variant_is_explicitly_unassessed():
    module = _module()
    assert module.variant_type("AT", "GCAT") == "COMPLEX"
    # A complex ALT is never matched by read similarity; an anchor read with
    # the ref base and no indel counts as reference-anchored, and the
    # row-level status (checked in build_rows) keeps the site UNASSESSED_COMPLEX
    # with blank numerics regardless.
    assert module.classify_observation(query_sequence="AT", query_position=0, indel=0, ref="AT", alts=["GCAT"]) == -1


def test_build_rows_expands_per_alt_with_shared_depth():
    module = _module()
    sites = [("chr1", 100, "A", ["C", "T"])]
    counts = [{"status": "ASSESSED", "error": "", "ref": 5, "alts": [2, 3], "other": 1}]
    rows = module.build_rows(sites, counts)
    assert [r["alt"] for r in rows] == ["C", "T"]
    assert [r["rna_alt_reads"] for r in rows] == ["2", "3"]
    # Depth keeps every sequenced allele category (review C3-02).
    assert all(r["rna_depth"] == "11" and r["rna_other_reads"] == "1" for r in rows)
    assert rows[0]["rna_vaf"] == "0.1818"
    assert all(r["count_status"] == "ASSESSED" for r in rows)


def test_build_rows_leaves_errors_and_complex_blank_not_zero():
    module = _module()
    sites = [
        ("chr1", 100, "A", ["T"]),
        ("chr9", 200, "AT", ["GCAT"]),
        ("chrZ", 300, "A", ["T"]),
    ]
    counts = [
        {"status": "ERROR", "error": "OSError: truncated"},
        {"status": "ASSESSED", "error": "", "ref": 4, "alts": [0], "other": 0},
        {"status": "CONTIG_NOT_IN_BAM", "error": "contig chrZ missing"},
    ]
    rows = module.build_rows(sites, counts)
    by_key = {r["variant_key"]: r for r in rows}
    err = by_key["chr1:100A>T"]
    assert err["count_status"] == "ERROR" and err["rna_depth"] == "" and err["rna_alt_reads"] == ""
    assert err["rna_count_error"] == "OSError: truncated"
    cx = by_key["chr9:200AT>GCAT"]
    assert cx["count_status"] == "UNASSESSED_COMPLEX" and cx["rna_depth"] == ""
    miss = by_key["chrZ:300A>T"]
    assert miss["count_status"] == "CONTIG_NOT_IN_BAM" and miss["rna_ref_reads"] == ""


def test_count_site_alleles_marks_errors_instead_of_zeros():
    # Review C3-01: the old except -> (0,0,0) made a failed pileup look like
    # "no expression"; failures must be reported, never counted.
    module = _module()

    class BoomBam:
        def pileup(self, *a, **k):
            raise ValueError("unknown reference")

    result = module.count_site_alleles(
        BoomBam(), "chr1", 100, "A", ["T"], references={"chr1"}, min_mapq=20, min_baseq=13,
    )
    assert result["status"] == "ERROR" and "unknown reference" in result["error"]

    class EmptyBam:
        def pileup(self, *a, **k):
            return iter(())

    ok = module.count_site_alleles(
        EmptyBam(), "chr1", 100, "A", ["T"], references={"chr1"}, min_mapq=20, min_baseq=13,
    )
    assert ok == {"status": "ASSESSED", "error": "", "ref": 0, "alts": [0], "other": 0}

    missing = module.count_site_alleles(
        EmptyBam(), "chrZ", 100, "A", ["T"], references={"chr1"}, min_mapq=20, min_baseq=13,
    )
    assert missing["status"] == "CONTIG_NOT_IN_BAM"


def test_preflight_refuses_missing_index_and_cross_build():
    # Review C3-03/C3-01 fail-fast checks.
    module = _module()

    class NoIndexBam:
        references = ("chr1",)

        def check_index(self):
            raise ValueError("No baim/bai or csi index")

    import pytest

    with pytest.raises(SystemExit, match="no usable index"):
        module.preflight_checks(NoIndexBam(), {"chr1"})

    class IndexedBam:
        references = ("chr1", "chr2")

        def check_index(self):
            return None

    with pytest.raises(SystemExit, match="coordinate/build mismatch"):
        module.preflight_checks(IndexedBam(), {"NM_001"})
    # Same contigs, no refusal; chr-prefix tolerance both directions.
    module.preflight_checks(IndexedBam(), {"1", "chr2"})


def _fake_pysam(bam_factory):
    return SimpleNamespace(AlignmentFile=bam_factory)


def test_main_end_to_end_with_fake_pysam(tmp_path, monkeypatch, capsys):
    module = _module()

    class Pr:
        def __init__(self, seq, qp, indel=0):
            self.is_refskip = False
            self.is_del = False
            self.query_position = qp
            self.indel = indel
            self.alignment = SimpleNamespace(query_sequence=seq)

    class Col:
        def __init__(self, pos, prs):
            self.reference_pos = pos
            self.pileups = prs

    class FakeBam:
        references = ("chr1",)

        def check_index(self):
            return None

        def pileup(self, chrom, start, end, **kwargs):
            assert (chrom, start, end) == ("chr1", 99, 100)
            return [Col(99, [Pr("TT", 0), Pr("CC", 0), Pr("GG", 0), Pr("AA", 0)])]

        def close(self):
            return None

    monkeypatch.setitem(sys.modules, "pysam", _fake_pysam(lambda path, mode: FakeBam()))
    vcf = tmp_path / "somatic.vcf"
    vcf.write_text(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\t.\tA\tC,T\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    out = tmp_path / "counts.tsv"
    monkeypatch.setattr(
        sys, "argv",
        ["prog", "--somatic-vcf", str(vcf), "--rna-bam", str(tmp_path / "fake.bam"), "--output-tsv", str(out)],
    )
    (tmp_path / "fake.bam").write_bytes(b"x")  # existence check only
    module.main()
    lines = out.read_text(encoding="utf-8").strip().splitlines()
    header = lines[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in lines[1:]]
    assert [r["alt"] for r in rows] == ["C", "T"]
    # one T read, one C read, one other (G), one ref (A) -> depth 4 shared
    assert rows[0]["rna_alt_reads"] == "1" and rows[1]["rna_alt_reads"] == "1"
    assert rows[0]["rna_depth"] == "4" and rows[0]["rna_other_reads"] == "1" and rows[0]["rna_ref_reads"] == "1"
    stdout = capsys.readouterr().out
    assert "cross-check" in stdout

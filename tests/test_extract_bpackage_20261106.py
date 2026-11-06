"""Regression tests for the B1/B2 extract-package fixes (B-E review).

Covers: B2-01 (no fabricated WT), B2-02 (ALT-index attribution),
B2-03 (stop-gained drop visibility), B2-04 (no tumor-column guessing),
B2-05 (compound consequence terms), B2-06 (FILTER ledger + forwarded flag),
B2-07 (dead sort key), B1-04 (post-run plugin verification),
B1-05 (CSQ header sniff now samples records).
"""

import pytest

import neoag.vep.annotate as vep_annotate
from neoag.tools.upstream import vcf_has_csq_annotations
from neoag.vep.extract_peptides import (
    MULTI_AA_SINGLE,
    MULTI_AA_SUBSTITUTION,
    _csq_rank,
    classify_multi_aa_flag,
    extract_variant_peptides_from_vcf,
    pick_csq_transcript,
)

CSQ_FIELDS = [
    "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature",
    "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position",
    "Protein_position", "Amino_acids", "Codons", "Existing_variation", "DISTANCE",
    "STRAND", "FLAGS", "SYMBOL_SOURCE", "HGNC_ID", "CANONICAL", "MANE", "MANE_SELECT",
    "TSL", "HGVS_OFFSET", "FrameshiftSequence", "WildtypeProtein",
]
CSQ_HEADER = (
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
    "Format: " + "|".join(CSQ_FIELDS) + '">'
)
WT = "ACDEFGHIKLMNPQRSTVWYACDEFG"


def csq_entry(allele, consequence, *, pos, aa, wt=WT, fs="", canonical="", mane=""):
    row = {name: "" for name in CSQ_FIELDS}
    row.update({
        "Allele": allele, "Consequence": consequence, "IMPACT": "MODERATE",
        "SYMBOL": "GENE1", "Gene": "ENSG0001", "Feature_type": "Transcript",
        "Feature": "ENST00001.1", "BIOTYPE": "protein_coding",
        "Protein_position": pos, "Amino_acids": aa, "CANONICAL": canonical,
        "MANE_SELECT": mane, "FrameshiftSequence": fs, "WildtypeProtein": wt,
    })
    return "|".join(row[name] for name in CSQ_FIELDS)


def _write_vcf(path, *, header_extra, records):
    path.write_text(
        "\n".join([CSQ_HEADER, header_extra] + records) + "\n", encoding="utf-8"
    )


SINGLE = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR"


def test_b2_02_multiallelic_metrics_follow_the_csq_allele(tmp_path):
    entry_a = csq_entry("A", "synonymous_variant", pos="4", aa="K/K")
    entry_g = csq_entry("G", "missense_variant", pos="4", aa="K/N")
    vcf = tmp_path / "ma.vcf"
    _write_vcf(
        vcf,
        header_extra=SINGLE,
        records=[
            f"chr1\t100\t.\tC\tA,G\t.\tPASS\tCSQ={entry_a},{entry_g}\t"
            "GT:AD:AF\t0/2:100,1,8:0.005,0.073"
        ],
    )
    out = tmp_path / "pep.tsv"
    summary = extract_variant_peptides_from_vcf(vcf, out, sample_id="T1", lengths=(8,))
    assert summary["peptide_rows"] > 0
    header = out.read_text(encoding="utf-8").split("\n")[0].split("\t")
    row = out.read_text(encoding="utf-8").strip().split("\n")[1].split("\t")
    # The peptide came from ALT[1] (G); VAF/AD must be read at index 1,
    # not at ALT[0] (A): A had 1/109 reads, G had 8/109.
    assert row[header.index("vaf")] == "0.0730"
    assert row[header.index("tumor_alt_count")] == "8"


def test_b2_04_multi_sample_vcf_without_names_raises(tmp_path):
    entry = csq_entry("T", "missense_variant", pos="4", aa="K/N")
    vcf = tmp_path / "amb.vcf"
    _write_vcf(
        vcf,
        header_extra="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\tSAMPLE_B",
        records=[f"chr1\t100\t.\tA\tT\t.\tPASS\tCSQ={entry}\tGT:AD:AF\t0/1:20,5:0.25\t0/0:30,0:0.0"],
    )
    out = tmp_path / "pep.tsv"
    with pytest.raises(ValueError, match="multi-sample VCF requires an explicit tumor"):
        extract_variant_peptides_from_vcf(vcf, out, sample_id="T1", lengths=(8,))
    # naming either column resolves it; the unnamed case must never guess.
    summary = extract_variant_peptides_from_vcf(
        vcf, out, sample_id="T1", lengths=(8,), tumor_sample_name="SAMPLE_B"
    )
    assert summary["tumor_sample_column"] == "SAMPLE_B"


def test_b2_05_compound_terms_are_not_complex_but_true_multi_stays(tmp_path):
    assert classify_multi_aa_flag("missense_variant&splice_region_variant", "K/N", "10") == MULTI_AA_SINGLE
    assert classify_multi_aa_flag("stop_gained&splice_region_variant", "K/*", "10") == MULTI_AA_SINGLE
    assert classify_multi_aa_flag("missense_variant&splice_region_variant", "MG/IS", "1388-1389") == MULTI_AA_SUBSTITUTION
    assert classify_multi_aa_flag("stop_gained&frameshift_variant", "K/*", "10") == "complex"

    entry = csq_entry("T", "missense_variant&splice_region_variant", pos="4", aa="K/N")
    vcf = tmp_path / "compound.vcf"
    _write_vcf(vcf, header_extra=SINGLE, records=[f"chr1\t100\t.\tA\tT\t.\tPASS\tCSQ={entry}\tGT:AD:AF\t0/1:20,5:0.25"])
    out = tmp_path / "pep.tsv"
    summary = extract_variant_peptides_from_vcf(
        vcf, out, sample_id="T1", lengths=(8,), exclude_multi_aa=True
    )
    assert summary["peptide_rows"] > 0  # not silently excluded any more
    assert summary["variants_skipped_multi_aa"] == 0


def test_b2_06_filter_funnel_is_recorded_and_flag_is_forwarded(tmp_path):
    entry = csq_entry("T", "missense_variant", pos="4", aa="K/N")
    vcf = tmp_path / "filters.vcf"
    _write_vcf(
        vcf,
        header_extra=SINGLE,
        records=[
            f"chr1\t100\t.\tA\tT\t.\tPASS\tCSQ={entry}\tGT:AD:AF\t0/1:20,5:0.25",
            f"chr1\t200\t.\tA\tT\t.\tLOWQUAL\tCSQ={entry}\tGT:AD:AF\t0/1:20,5:0.25",
        ],
    )
    out = tmp_path / "pep.tsv"
    summary = extract_variant_peptides_from_vcf(vcf, out, sample_id="T1", lengths=(8,))
    assert summary["filter_status_counts"] == {"LOWQUAL": 1, "PASS": 1}
    assert summary["variants_parsed"] == 1  # PASS-only is still the default
    summary_all = extract_variant_peptides_from_vcf(
        vcf, out, sample_id="T1", lengths=(8,), pass_only=False
    )
    assert summary_all["variants_parsed"] == 2


def test_b2_03_stop_gained_stays_peptide_negative_but_is_counted(tmp_path):
    entry = csq_entry("T", "stop_gained", pos="4", aa="K/*")
    vcf = tmp_path / "stop.vcf"
    _write_vcf(vcf, header_extra=SINGLE, records=[f"chr1\t100\t.\tA\tT\t.\tPASS\tCSQ={entry}\tGT:AD:AF\t0/1:20,5:0.25"])
    out = tmp_path / "pep.tsv"
    summary = extract_variant_peptides_from_vcf(vcf, out, sample_id="T1", lengths=(8,))
    assert summary["peptide_rows"] == 0
    assert summary["variants_no_peptides_by_consequence"].get("stop_gained") == 1


def test_b2_01_rows_mark_wt_unavailable_when_no_plugin_context(tmp_path):
    entry = csq_entry("T", "frameshift_variant", pos="4", aa="K/MA*", wt="", fs="MAASAANPQRST")
    vcf = tmp_path / "fs.vcf"
    _write_vcf(vcf, header_extra=SINGLE, records=[f"chr1\t100\t.\tA\tTG\t.\tPASS\tCSQ={entry}\tGT:AD:AF\t0/1:20,5:0.25"])
    out = tmp_path / "pep.tsv"
    summary = extract_variant_peptides_from_vcf(vcf, out, sample_id="T1", lengths=(4,))
    assert summary["peptide_rows"] > 0
    text = out.read_text(encoding="utf-8")
    assert "WT_UNAVAILABLE" in text
    assert "ACDEFG" not in text  # no fabricated 42-repeat backbone anywhere


def test_b2_01_summary_exposes_drop_ledgers(tmp_path):
    entry_no_context = csq_entry("T", "missense_variant", pos="4", aa="K/N", wt="")
    vcf = tmp_path / "noproto.vcf"
    _write_vcf(vcf, header_extra=SINGLE, records=[f"chr1\t100\t.\tA\tT\t.\tPASS\tCSQ={entry_no_context}\tGT:AD:AF\t0/1:20,5:0.25"])
    out = tmp_path / "pep.tsv"
    summary = extract_variant_peptides_from_vcf(vcf, out, sample_id="T1", lengths=(8,))
    assert summary["variants_dropped_no_protein_context"] == 1
    assert summary["peptide_rows"] == 0
    assert summary["csq_wildtype_column"] == "yes"


def test_b2_07_csq_rank_has_two_slots_and_stable_ties():
    idx = {name: i for i, name in enumerate(CSQ_FIELDS)}
    entry = csq_entry("T", "missense_variant", pos="4", aa="K/N")
    rank = _csq_rank(entry.split("|"), idx)
    assert isinstance(rank, tuple) and len(rank) == 2
    first = csq_entry("T", "missense_variant", pos="4", aa="K/N")
    second = csq_entry("T", "missense_variant", pos="9", aa="S/A")
    picked = pick_csq_transcript([first.split("|"), second.split("|")], idx)
    assert picked is not None
    assert picked == first.split("|")  # equal keys keep CSQ entry order (deterministic)


def test_b1_04_csq_plugin_columns_detected_and_missing_plugins_raise(tmp_path, monkeypatch):
    vcf = tmp_path / "plug.vcf"
    _write_vcf(vcf, header_extra=SINGLE, records=[])
    assert vep_annotate.csq_plugin_columns(vcf) == {"WildtypeProtein": True, "FrameshiftSequence": True}

    plain_header = '##INFO=<ID=CSQ,Number=.,Type=String,Description="x. Format: Allele|Consequence|Amino_acids">'
    noplugin = tmp_path / "noplugin.vcf"
    out_vcf = tmp_path / "annot.vcf"

    def fake_run(cmd, **kwargs):
        out_vcf.write_text("\n".join([plain_header, SINGLE]) + "\n", encoding="utf-8")

    monkeypatch.setattr(vep_annotate.subprocess, "run", fake_run)
    with pytest.raises(RuntimeError, match="plugin failed to load"):
        vep_annotate.run_vep_pvacseq_annotate(
            input_vcf="in.vcf", output_vcf=str(out_vcf), reference_fasta="ref.fa",
            plugins_dir=str(tmp_path / "plugins"), index_vcf=False,
        )


def test_b1_04_columns_recorded_when_no_plugins_requested(tmp_path, monkeypatch):
    plain_header = '##INFO=<ID=CSQ,Number=.,Type=String,Description="x. Format: Allele|Consequence|Amino_acids">'
    out_vcf = tmp_path / "annot2.vcf"

    def fake_run(cmd, **kwargs):
        out_vcf.write_text("\n".join([plain_header, SINGLE]) + "\n", encoding="utf-8")

    monkeypatch.setattr(vep_annotate.subprocess, "run", fake_run)
    info = vep_annotate.run_vep_pvacseq_annotate(
        input_vcf="in.vcf", output_vcf=str(out_vcf), reference_fasta="ref.fa", index_vcf=False,
    )
    assert info["vep_plugin_columns"] == "WildtypeProtein=no;FrameshiftSequence=no"


def test_b1_05_header_only_csf_sniff_now_samples_records(tmp_path):
    declared_only = tmp_path / "declared.vcf"
    declared_only.write_text(
        "\n".join([CSQ_HEADER, SINGLE, "chr1\t100\t.\tA\tT\t.\tPASS\tDP=10\tGT:AD\t0/1:20,5"]) + "\n",
        encoding="utf-8",
    )
    assert not vcf_has_csq_annotations(declared_only)

    annotated = tmp_path / "annotated.vcf"
    entry = csq_entry("T", "missense_variant", pos="4", aa="K/N")
    annotated.write_text(
        "\n".join([CSQ_HEADER, SINGLE, f"chr1\t100\t.\tA\tT\t.\tPASS\tCSQ={entry}\tGT:AD\t0/1:20,5"]) + "\n",
        encoding="utf-8",
    )
    assert vcf_has_csq_annotations(annotated)

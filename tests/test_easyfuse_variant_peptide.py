"""Tests for unified EasyFuse catalog builder."""

from pathlib import Path

from neoag.adapters.easyfuse_variant_peptide import (
    EasyFusePeptideConfig,
    build_easyfuse_catalog,
    build_fusion_centered_minigene,
    easyfuse_to_variant_peptide_rows,
    fusion_generation_fingerprint,
    sliding_fusion_neo_peptides,
    write_easyfuse_qc_tables,
)
from neoag.adapters.variant_peptide_adapter import (
    _catalog_rows_to_raw_peptides,
    run_variant_peptide_upstream,
)
from neoag.vep.extract_peptides import (
    build_proteome_kmer_index,
    read_fasta_sequences,
)

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "data/fixtures/easyfuse_fusions.pass.tsv"

CSQ_HEADER = (
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
    "Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|"
    "HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|"
    "Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|MANE_SELECT|"
    'TSL|HGVS_OFFSET|FrameshiftSequence|WildtypeProtein">'
)


def _write_mini_vcf(path: Path) -> None:
    path.write_text(
        "\n".join([
            CSQ_HEADER,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "chr1\t100\t.\tA\tT\t.\tPASS\tCSQ=T|missense_variant|MODERATE|GENE1|ENSG0001|"
            "Transcript|ENST00001.1|protein_coding|1/2||ENST00001.1:c.10A>T|"
            "ENSP00001.1:p.Lys4Asn|10|10|4|K/N|||1||HGNC|HGNC:1|YES||MANE_Select|"
            "NM_0001.1|1|||ACDEFGHIKLMNPQRSTVWYACDEFG",
        ]) + "\n",
        encoding="utf-8",
    )


def test_sliding_fusion_neo_peptides_dedupes_windows():
    windows = sliding_fusion_neo_peptides("ABCDEFGHIJK", (8, 9), bp_pos=5)
    peptides = {w["mutant_peptide"] for w in windows}
    assert "ABCDEHGI" not in peptides
    assert len(peptides) == len(windows)
    assert all(8 <= len(p) <= 9 for p in peptides)


def test_sliding_fusion_windows_record_both_junction_sides():
    windows = sliding_fusion_neo_peptides("ACDEFGHIKLMN", (9,), bp_pos=5)
    spanning = [row for row in windows if row["crosses_junction"] == "yes"]
    assert spanning
    assert all(row["fusion_left_peptide"] and row["fusion_right_peptide"] for row in spanning)
    assert all(
        row["fusion_junction_display"]
        == row["fusion_left_peptide"] + "|" + row["fusion_right_peptide"]
        for row in spanning
    )




def test_build_fusion_centered_minigene():
    mg = build_fusion_centered_minigene("ABCDEFGHIJKLMN", peptide_start_aa=4, peptide_end_aa=11, mini_len=2)
    assert mg == "BC|DEFGHIJK|LM"
    assert mg.split("|")[1] == "DEFGHIJK"


def test_easyfuse_catalog_minigene_centers_each_fusion_peptide():
    result = build_easyfuse_catalog(FIXTURE, "S1", "default", lengths=(8,), peptide_cfg=EasyFusePeptideConfig(junction_only=False, mini_len=3))
    assert result.catalog_rows
    assert all(row["minigene"].split("|")[1] == row["mutant_peptide"] for row in result.catalog_rows)

def test_build_easyfuse_catalog_strict_defaults():
    result = build_easyfuse_catalog(FIXTURE, "S1", "default", lengths=(8, 9))
    assert len(result.events) == 2
    assert result.catalog_rows
    assert all(r["multi_aa_flag"] == "fusion_neo" for r in result.catalog_rows)
    assert all(r["peptide_source"] == "easyfuse" for r in result.catalog_rows)
    assert all(
        r["fusion_generation_method"].startswith("easyfuse_neo_peptide_sliding_window")
        for r in result.catalog_rows
    )
    assert result.summary_qc["junction_only"] == "yes"
    assert result.summary_qc["dedup_per_event"] == "yes"
    assert len(result.filter_qc) == 4
    assert len(result.collapse_qc) == 2


def test_easyfuse_to_variant_peptide_rows_wrapper_has_no_legacy_peptides():
    rows, events, legacy = easyfuse_to_variant_peptide_rows(FIXTURE, "S1", "default", lengths=(8, 9))
    assert events
    assert rows
    assert legacy == []


def test_isoform_strategy_max_junction_reads(tmp_path):
    dup = tmp_path / "dup.tsv"
    dup.write_text(
        "\t".join([
            "BPID", "FTID", "Fusion_Gene", "Breakpoint1", "Breakpoint2", "type", "frame",
            "neo_peptide_sequence", "neo_peptide_sequence_bp", "prediction_prob",
            "prediction_class", "ft_junc_cnt", "ft_anch_cnt",
        ]) + "\n"
        + "\t".join([
            "bp1", "FT_short", "GENE1_GENE2", "chr1:1:+", "chr2:2:+", "trans", "in_frame",
            "SHORTSEQ", "3", "0.9", "positive", "5", "12",
        ]) + "\n"
        + "\t".join([
            "bp1", "FT_long", "GENE1_GENE2", "chr1:1:+", "chr2:2:+", "trans", "in_frame",
            "SHORTSEQEXTRA", "3", "0.7", "positive", "20", "12",
        ]) + "\n",
        encoding="utf-8",
    )
    result = build_easyfuse_catalog(
        dup,
        "S1",
        "default",
        lengths=(8,),
        peptide_cfg=EasyFusePeptideConfig(
            junction_only=False,
            dedup_per_event=True,
            isoform_strategy="max_junction_reads",
        ),
    )
    assert result.collapse_qc[0]["selected_ftid"] == "FT_long"


def test_write_easyfuse_qc_tables(tmp_path):
    result = build_easyfuse_catalog(FIXTURE, "S1", "default", lengths=(8,))
    paths = write_easyfuse_qc_tables(result, tmp_path)
    assert Path(paths["filter_qc"]).is_file()
    assert Path(paths["collapse_qc"]).is_file()
    assert Path(paths["summary_qc"]).is_file()


def test_catalog_rows_to_raw_peptides_links_fusion_events():
    result = build_easyfuse_catalog(FIXTURE, "S1", "default", lengths=(8,))
    event_map = {e["event_id"]: e for e in result.events}
    peptides = _catalog_rows_to_raw_peptides(
        result.catalog_rows[:2],
        event_map,
        sample_id="S1",
        hla_alleles=["HLA-A*02:01"],
        source_tool="EasyFuse",
    )
    assert peptides
    assert peptides[0]["event_type"] == "Fusion"
    assert peptides[0]["source_tool"] == "EasyFuse"
    assert peptides[0]["crosses_junction"] == "yes"


def test_easyfuse_catalog_supports_optional_12mer():
    windows = sliding_fusion_neo_peptides("ACDEFGHIKLMNPQRST", (12,), bp_pos=8)
    assert windows and all(len(row["mutant_peptide"]) == 12 for row in windows)
    result = build_easyfuse_catalog(FIXTURE, "S1", "default", lengths=(8, 9, 10, 11, 12))
    assert result.summary_qc["peptide_lengths"] == "8,9,10,11,12"


def test_run_variant_peptide_upstream_merges_easyfuse(tmp_path):
    vcf = tmp_path / "mini.vcf"
    _write_mini_vcf(vcf)
    cfg = {
        "inputs": {
            "variant_peptide_extraction": True,
            "variant_peptide_lengths": "8",
            "variant_peptide_mini_len": 8,
            "variant_peptide_filter_normal_proteome": False,
            "easyfuse_pass_csv": str(FIXTURE),
            "tumor_sample_name": "TUMOR",
        }
    }
    parsed = tmp_path / "parsed"
    tools = tmp_path / "tools"
    parsed.mkdir()
    tools.mkdir()

    outs = run_variant_peptide_upstream(
        cfg,
        variants_vcf=vcf,
        parsed_dir=parsed,
        tools_dir=tools,
        sample_id="S1",
        profile_name="default",
        hla_alleles=["HLA-A*02:01"],
    )
    assert int(outs["variant_peptide_rows_vcf"]) > 0
    assert int(outs["variant_peptide_rows_easyfuse"]) > 0
    assert Path(outs["easyfuse_variant_peptides"]).is_file()
    assert Path(outs["filter_qc"]).is_file()
    assert Path(outs["summary_qc"]).is_file()
    catalog = (tmp_path / "tools" / "variant_peptides.tsv").read_text(encoding="utf-8")
    assert "fusion_neo" in catalog
    assert "\teasyfuse\t" in catalog or catalog.endswith("easyfuse\n")


# ---------------------------------------------------------------------------
# 2026-11 fusion audit regressions: junction eligibility, normal-parent
# exact-match flag, coordinate ambiguity, provenance and cache fingerprint.
# ---------------------------------------------------------------------------

_EF_COLS = [
    "BPID", "FTID", "Fusion_Gene", "Breakpoint1", "Breakpoint2", "type", "frame",
    "neo_peptide_sequence", "neo_peptide_sequence_bp",
    "fusion_protein_sequence", "fusion_protein_sequence_bp",
    "prediction_prob", "prediction_class", "ft_junc_cnt", "ft_anch_cnt",
]


def _fusion_table(tmp_path, **overrides):
    row = {
        "BPID": "chr1:100:+_chr2:200:-", "FTID": "FT_A_A::B", "Fusion_Gene": "A::B",
        "Breakpoint1": "chr1:100:+", "Breakpoint2": "chr2:200:-",
        "type": "trans", "frame": "in_frame",
        "neo_peptide_sequence": "", "neo_peptide_sequence_bp": "",
        "fusion_protein_sequence": "", "fusion_protein_sequence_bp": "",
        "prediction_prob": "0.9", "prediction_class": "positive",
        "ft_junc_cnt": "10", "ft_anch_cnt": "30",
    }
    row.update(overrides)
    path = tmp_path / "ef_fusion.tsv"
    path.write_text(
        "\t".join(_EF_COLS) + "\n"
        + "\t".join(str(row[c]) for c in _EF_COLS) + "\n",
        encoding="utf-8",
    )
    return path


# 60-residue upstream segment + distinct downstream segment; the junction
# sits at fusion-protein position 60, far from any neo-scale coordinate, so
# windows anchored at the N-terminus (the deployed defect) cannot fake them.
# Only letters in the sanitizer whitelist ACDEFGHIKLMNPQRSTVWY may be used.
_DEEP_LEFT = "ACDEFGHIKLMNPQRSTVWY" * 3
_DEEP_RIGHT = "KMSEQWAPRTVYHCGDNFVC"
_DEEP_FUSION = _DEEP_LEFT + _DEEP_RIGHT
_DEEP_NEO = _DEEP_LEFT[-8:] + _DEEP_RIGHT[:5]  # 13 aa, junction after res 8


def _deep_row(**overrides):
    row = {
        "fusion_protein_sequence": _DEEP_FUSION,
        "fusion_protein_sequence_bp": "60",
        "neo_peptide_sequence": _DEEP_NEO,
        "neo_peptide_sequence_bp": "8",  # implied neo start 60-8=52 >= 0
    }
    row.update(overrides)
    return row


def test_fusion_protein_windows_straddle_reported_breakpoint(tmp_path):
    # Regression: windows must straddle the VALIDATED fusion-protein
    # breakpoint (residue 60), never a neo-scale coordinate near the N-term.
    table = _fusion_table(tmp_path, **_deep_row())
    result = build_easyfuse_catalog(table, "S1", "default", lengths=(8,))
    rows = result.catalog_rows
    assert rows
    expected_starts = {s + 1 for s in range(53, 60)}
    assert {int(r["peptide_start_aa"]) for r in rows} == expected_starts
    for r in rows:
        assert int(r["peptide_start_aa"]) <= 60 < int(r["peptide_end_aa"])
        assert r["mutant_peptide"] == _DEEP_FUSION[
            int(r["peptide_start_aa"]) - 1: int(r["peptide_end_aa"])
        ]
        assert r["crosses_junction"] == "yes"
        assert r["peptide_context"] == "fusion_protein"
        assert r["peptide_context_status"] == "CONTEXT_FUSION_PROTEIN"
        assert r["peptide_context_ambiguous_coordinates"] == "no"
        assert r["junction_eligibility"] == "ELIGIBLE_JUNCTION_OVERLAP"
        assert r["junction_protein_position"] == "60"


def test_native_parent_windows_flagged_and_novel_kept(tmp_path):
    # Junction context "MNACDEFGHIKLMNP", junction after residue 8.
    # Proteome carries the upstream parent (first 13 residues) and a
    # downstream parent starting at residue 8: windows matching either
    # parent exactly are flagged native; the genuinely chimeric one is not.
    fp = "MNACDEFGHIKLMNP"
    table = _fusion_table(
        tmp_path,
        fusion_protein_sequence=fp,
        fusion_protein_sequence_bp="8",
        neo_peptide_sequence="ACDEFGHI",
        neo_peptide_sequence_bp="8",
    )
    fasta = tmp_path / "normal_proteins.fa"
    fasta.write_text(
        ">left_parent\nMNACDEFGHIKL\n>right_parent\nGHIKLMNPQRS\n",
        encoding="utf-8",
    )
    index = build_proteome_kmer_index(read_fasta_sequences(fasta), (8,))
    result = build_easyfuse_catalog(
        table, "S1", "default", lengths=(8,),
        peptide_cfg=EasyFusePeptideConfig(
            normal_proteome_index=index, normal_proteome_source="unit-test",
        ),
    )
    rows = {int(r["peptide_start_aa"]): r for r in result.catalog_rows}
    assert set(rows) == set(range(2, 9))  # 0-based starts 1..7
    # starts 2-5 exact-match the upstream parent, start 8 the downstream
    # parent; starts 6-7 are the genuinely chimeric windows.
    expected_native = {2, 3, 4, 5, 8}
    for start, r in rows.items():
        want = "yes" if start in expected_native else "no"
        assert r["mutant_peptide"] == fp[start - 1:start + 7]
        assert r["in_normal_proteome"] == want, (start, r["mutant_peptide"])
        assert r["contains_novel_aa"] == ("no" if want == "yes" else "yes")
        assert r["normal_proteome_source"] == "unit-test"
        assert r["junction_eligibility"] == "ELIGIBLE_JUNCTION_OVERLAP"


def test_missing_reference_proteome_stays_unassessed(tmp_path):
    table = _fusion_table(
        tmp_path,
        fusion_protein_sequence="MNACDEFGHIKLMNP",
        fusion_protein_sequence_bp="8",
        neo_peptide_sequence="ACDEFGHI",
        neo_peptide_sequence_bp="8",
    )
    result = build_easyfuse_catalog(table, "S1", "default", lengths=(8,))
    assert result.catalog_rows
    for r in result.catalog_rows:
        assert r["in_normal_proteome"] == "UNASSESSED"
        assert r["contains_novel_aa"] == "UNASSESSED"
        assert r["normal_proteome_source"] == "not_provided"


def test_contradictory_coordinates_yield_unassessed_eligibility(tmp_path):
    # Both spaces carry a breakpoint but the implied neo start is negative:
    # the spaces disagree about where the junction sits, so no window may
    # claim junction eligibility.
    table = _fusion_table(
        tmp_path, **_deep_row(neo_peptide_sequence_bp="65")
    )
    result = build_easyfuse_catalog(table, "S1", "default", lengths=(8,))
    rows = result.catalog_rows
    assert rows
    for r in rows:
        assert r["peptide_context_status"] == "CONTEXT_FUSION_PROTEIN"
        assert r["peptide_context_ambiguous_coordinates"] == "yes"
        assert r["junction_eligibility"] == "UNASSESSED_AMBIGUOUS_COORDINATES"
    ctx_row = next(
        c for c in result.coordinate_qc
        if c["context_status"] == "CONTEXT_FUSION_PROTEIN"
    )
    assert "implied neo start" in ctx_row["context_detail"]


def test_context_invalid_row_emits_no_peptides(tmp_path):
    # A breakpoint coordinate that is valid in NEITHER space must produce no
    # windows at all and be flagged, never windows in a foreign coordinate
    # space (the deployed defect class).
    table = _fusion_table(
        tmp_path,
        neo_peptide_sequence="ACDEFGHIKLMNP",
        neo_peptide_sequence_bp="99",  # invalid for the 13-aa neo sequence
    )
    result = build_easyfuse_catalog(table, "S1", "default", lengths=(8, 9))
    assert result.catalog_rows == []
    invalid = [
        c for c in result.coordinate_qc if c["context_status"] == "CONTEXT_INVALID"
    ]
    assert len(invalid) == 1
    assert invalid[0]["junction_window_rows"] == "0"
    assert result.summary_qc["peptide_context_invalid_rows"] == "1"


def test_neo_fallback_rows_keep_eligibility_provenance(tmp_path):
    table = _fusion_table(
        tmp_path,
        neo_peptide_sequence="ACDEFGHIKLMNP",
        neo_peptide_sequence_bp="8",
    )
    result = build_easyfuse_catalog(table, "S1", "default", lengths=(8, 9))
    rows = result.catalog_rows
    assert rows
    for r in rows:
        assert r["peptide_context_status"] == "CONTEXT_NEO_FALLBACK"
        assert r["peptide_context_ambiguous_coordinates"] == "no"
        assert r["junction_eligibility"] == "ELIGIBLE_JUNCTION_OVERLAP"
        assert r["junction_protein_position"] == r["peptide_context_bp"] == "8"


def test_out_of_frame_fusion_keeps_provenance_and_eligibility(tmp_path):
    # The frame label is provenance; it must not fabricate or erase the
    # coordinate-based eligibility verdict.
    table = _fusion_table(tmp_path, **_deep_row(frame="out_of_frame"))
    result = build_easyfuse_catalog(table, "S1", "default", lengths=(8,))
    rows = result.catalog_rows
    assert rows
    for r in rows:
        assert int(r["peptide_start_aa"]) <= 60 < int(r["peptide_end_aa"])
        assert r["junction_eligibility"] == "ELIGIBLE_JUNCTION_OVERLAP"
    events = {e["event_id"]: e for e in result.events}
    assert events


def test_alternative_context_choice_prefers_validated_fusion_protein(tmp_path):
    # Both an assembled protein coordinate and a neo coordinate are valid:
    # the assembled (alternative-ORF-resolved) protein wins, and its
    # coordinate space is recorded per peptide.
    table = _fusion_table(tmp_path, **_deep_row())
    result = build_easyfuse_catalog(table, "S1", "default", lengths=(8,))
    for r in result.catalog_rows:
        assert r["peptide_context"] == "fusion_protein"
        assert r["peptide_context_bp"] == "60"


def test_generation_fingerprint_is_deterministic_and_visible(tmp_path):
    fp_a = fusion_generation_fingerprint(
        (8, 9), junction_only=True, dedup_per_event=True,
        isoform_strategy="max_junction_reads",
    )
    fp_b = fusion_generation_fingerprint(
        (8, 9), junction_only=True, dedup_per_event=True,
        isoform_strategy="max_junction_reads",
    )
    assert fp_a == fp_b and fp_a.startswith("EFGEN1-")
    assert fp_a != fusion_generation_fingerprint(
        (8, 9), junction_only=False, dedup_per_event=True,
        isoform_strategy="max_junction_reads",
    )
    assert fp_a != fusion_generation_fingerprint(
        (8, 9), junction_only=True, dedup_per_event=True,
        isoform_strategy="max_junction_reads",
        normal_proteome_provided=True, normal_proteome_source="pinned.fa",
    )

    r1 = build_easyfuse_catalog(FIXTURE, "S1", "default", lengths=(8,))
    r2 = build_easyfuse_catalog(FIXTURE, "S1", "default", lengths=(8,))
    stamp = r1.summary_qc["generation_fingerprint"]
    assert stamp == r2.summary_qc["generation_fingerprint"]
    paths = write_easyfuse_qc_tables(r1, tmp_path)
    header = Path(paths["summary_qc"]).read_text(encoding="utf-8").splitlines()[0]
    assert "generation_fingerprint" in header.split("\t")
    assert stamp in Path(paths["summary_qc"]).read_text(encoding="utf-8")


def test_upstream_run_carries_fusion_provenance_into_raw_tables(tmp_path):
    vcf = tmp_path / "mini.vcf"
    _write_mini_vcf(vcf)
    cfg = {
        "inputs": {
            "variant_peptide_extraction": True,
            "variant_peptide_lengths": "8",
            "variant_peptide_mini_len": 8,
            "easyfuse_pass_csv": str(FIXTURE),
            "tumor_sample_name": "TUMOR",
        }
    }
    parsed = tmp_path / "parsed"
    tools = tmp_path / "tools"
    parsed.mkdir()
    tools.mkdir()
    outs = run_variant_peptide_upstream(
        cfg,
        variants_vcf=vcf,
        parsed_dir=parsed,
        tools_dir=tools,
        sample_id="S1",
        profile_name="default",
        hla_alleles=["HLA-A*02:01"],
    )
    catalog_text = Path(outs["variant_peptides"]).read_text(encoding="utf-8")
    catalog_header = catalog_text.splitlines()[0].split("\t")
    for col in (
        "junction_eligibility", "in_normal_proteome",
        "peptide_context_ambiguous_coordinates", "normal_proteome_source",
    ):
        assert col in catalog_header
    assert "ELIGIBLE_JUNCTION_OVERLAP" in catalog_text

    raw_text = Path(outs["raw_peptides"]).read_text(encoding="utf-8")
    raw_header = raw_text.splitlines()[0].split("\t")
    assert "junction_eligibility" in raw_header
    assert "in_normal_proteome" in raw_header
    # No proteome FASTA was configured: every row must say UNASSESSED,
    # never an implied "no".
    raw_rows = [dict(zip(raw_header, line.split("\t")))
                for line in raw_text.splitlines()[1:]]
    assert raw_rows
    for row in raw_rows:
        assert row["in_normal_proteome"] == "UNASSESSED"


def test_snv_only_run_has_no_fusion_columns(tmp_path):
    vcf = tmp_path / "mini.vcf"
    _write_mini_vcf(vcf)
    cfg = {
        "inputs": {
            "variant_peptide_extraction": True,
            "variant_peptide_lengths": "8",
            "variant_peptide_mini_len": 8,
            "tumor_sample_name": "TUMOR",
        }
    }
    parsed = tmp_path / "parsed"
    tools = tmp_path / "tools"
    parsed.mkdir()
    tools.mkdir()
    run_variant_peptide_upstream(
        cfg,
        variants_vcf=vcf,
        parsed_dir=parsed,
        tools_dir=tools,
        sample_id="S1",
        profile_name="default",
        hla_alleles=["HLA-A*02:01"],
    )
    catalog_header = (tools / "variant_peptides.tsv").read_text(
        encoding="utf-8"
    ).splitlines()[0].split("\t")
    assert "junction_eligibility" not in catalog_header
    assert "in_normal_proteome" in catalog_header
    rows = [
        dict(zip(catalog_header, line.split("\t")))
        for line in (tools / "variant_peptides.tsv").read_text(
            encoding="utf-8"
        ).splitlines()[1:]
    ]
    assert rows
    # VCF path without a reference proteome: unassessed, not fabricated "no".
    assert all(r["in_normal_proteome"] == "UNASSESSED" for r in rows)


# ---------------------------------------------------------------------------
# 2026-11-10 follow-up investigation regressions: which inputs can produce
# peptides made ENTIRELY of downstream-partner residues (downstream_probe.py
# cases A-F). Each test pins the corrected tree's behaviour and, where the
# corrected tree deliberately keeps a documented residual gap, pins that too.
# Sequences use only sanitizer-whitelist letters ACDEFGHIKLMNPQRSTVWY.
# ---------------------------------------------------------------------------

_HEAD22 = "MNACDEFGHIKLMNPQRSTVWY"          # 22 aa
_TAIL20 = "KMSEQWAPRTVYHCGDNFVC"            # 20 aa


def test_downstream_first_context_stays_eligible_documents_residual_gap(tmp_path):
    # Residual gap (documented, fix requires parent-order verification): a
    # context that starts with the downstream partner but self-reports a
    # valid internal junction is consumed unvalidated. The coordinate fix
    # guarantees windows straddle the ASSEMBLED junction (so no window can
    # lie entirely in the N-terminal head), but eligibility cannot be
    # denied from sequence alone.
    bad = _TAIL20 + _HEAD22  # label says head gene first; sequence says otherwise
    table = _fusion_table(
        tmp_path,
        fusion_protein_sequence=bad,
        fusion_protein_sequence_bp="20",
    )
    result = build_easyfuse_catalog(table, "S1", "default", lengths=(8,))
    rows = result.catalog_rows
    assert rows
    for r in rows:
        assert int(r["peptide_start_aa"]) <= 20 < int(r["peptide_end_aa"])
        # nothing emitted from the pure-downstream N-terminal head band
        assert int(r["peptide_end_aa"]) > 20
        assert r["junction_eligibility"] == "ELIGIBLE_JUNCTION_OVERLAP"
    assert min(int(r["peptide_start_aa"]) for r in rows) == 14  # 20-8+2


def test_internal_start_assembly_never_emits_head_natives(tmp_path):
    # Deployed defect: windows anchored at constant position 13 emitted the
    # assembled head verbatim (EWSR1::WT1 RTHTGKTSEK = WT1 403-412 pattern).
    # Corrected: windows sit at the reported junction 22, so no pure-head
    # substring is emitted at all - even when the head IS a normal protein.
    fusion = _HEAD22 + _TAIL20
    table = _fusion_table(
        tmp_path,
        fusion_protein_sequence=fusion,
        fusion_protein_sequence_bp="22",
    )
    fasta = tmp_path / "head.fa"
    fasta.write_text(f">head_parent\n{_HEAD22}\n", encoding="utf-8")
    index = build_proteome_kmer_index(read_fasta_sequences(fasta), (8,))
    result = build_easyfuse_catalog(
        table, "S1", "default", lengths=(8,),
        peptide_cfg=EasyFusePeptideConfig(
            normal_proteome_index=index, normal_proteome_source="unit-test",
        ),
    )
    rows = result.catalog_rows
    assert rows
    emitted = {r["mutant_peptide"] for r in rows}
    assert not any(pep in _HEAD22 for pep in emitted)
    head_band = {fusion[i:i + 8] for i in range(2, 14)}  # deployed N-term band
    assert emitted.isdisjoint(head_band)
    for r in rows:
        assert r["in_normal_proteome"] == "no"
        assert int(r["peptide_start_aa"]) <= 22 < int(r["peptide_end_aa"])


def test_inconsistent_neo_bp13_cannot_reach_downstream_only_windows(tmp_path):
    # An inconsistent explicit neo bp (13 while the fusion junction is 10)
    # used to open the fully-downstream channel (probe B2: 12 native
    # downstream peptides). Corrected: the validated fusion coordinate wins,
    # the disagreement is flagged, and no window starts past the junction.
    fusion = "ACDEFGHIKL" + _TAIL20  # junction after residue 10
    table = _fusion_table(
        tmp_path,
        fusion_protein_sequence=fusion,
        fusion_protein_sequence_bp="10",
        neo_peptide_sequence=fusion[:26],
        neo_peptide_sequence_bp="13",
    )
    fasta = tmp_path / "down.fa"
    fasta.write_text(f">down_parent\n{_TAIL20}\n", encoding="utf-8")
    index = build_proteome_kmer_index(read_fasta_sequences(fasta), (8,))
    result = build_easyfuse_catalog(
        table, "S1", "default", lengths=(8,),
        peptide_cfg=EasyFusePeptideConfig(
            normal_proteome_index=index, normal_proteome_source="unit-test",
        ),
    )
    rows = result.catalog_rows
    assert rows
    assert _TAIL20[:8] not in {r["mutant_peptide"] for r in rows}
    for r in rows:
        assert int(r["peptide_start_aa"]) <= 10  # never downstream-only
        assert int(r["peptide_end_aa"]) > 10
        assert r["peptide_context_ambiguous_coordinates"] == "yes"
        assert r["junction_eligibility"] == "UNASSESSED_AMBIGUOUS_COORDINATES"
        assert r["in_normal_proteome"] == "no"  # every window is chimeric


def test_multi_ftid_collapse_selects_reads_regardless_of_context(tmp_path):
    # Probe E: one BPID, two FTID hypotheses. Selection (max_junction_reads)
    # is by read support only - a downstream-first hypothesis with more
    # reads wins and supplies all peptides under the shared event label.
    # Pinned as a documented residual gap of the current design.
    bad = _TAIL20 + _DEEP_LEFT[:20]
    path = tmp_path / "ef_two.tsv"
    row_tpl = {
        "BPID": "BP_X", "Fusion_Gene": "A::B",
        "Breakpoint1": "chr1:100:+", "Breakpoint2": "chr2:200:-",
        "type": "trans", "frame": "in_frame",
        "neo_peptide_sequence": "", "neo_peptide_sequence_bp": "",
        "prediction_prob": "0.9", "prediction_class": "positive",
        "ft_anch_cnt": "30",
    }
    lines = ["\t".join(_EF_COLS)]
    for ftid, protein, bp, junc in (
            ("FT_good", _DEEP_FUSION, "60", "5"),
            ("FT_bad", bad, "20", "20")):
        row = dict(row_tpl, FTID=ftid, fusion_protein_sequence=protein,
                   fusion_protein_sequence_bp=bp, ft_junc_cnt=junc)
        lines.append("\t".join(str(row[c]) for c in _EF_COLS))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    result = build_easyfuse_catalog(path, "S1", "default", lengths=(8,))
    assert result.collapse_qc
    qc = result.collapse_qc[0]
    assert qc["selected_ftid"] == "FT_bad"
    assert qc["rejected_ftids"] == "FT_good"
    rows = result.catalog_rows
    assert rows and all(r["fusion_transcript_id"] == "FT_bad" for r in rows)
    for r in rows:
        assert int(r["peptide_start_aa"]) <= 20 < int(r["peptide_end_aa"])
        assert r["gene"] == "A::B"  # label from the event, protein from FT_bad


def test_shared_segment_junction_window_is_native_and_still_eligible(tmp_path):
    # Probe F: a window straddling a GENUINE junction can exact-match a
    # normal protein when the parents share sequence. The native flag must
    # not silently revoke junction eligibility (orthogonal columns), and
    # junction eligibility must not pretend the window is novel.
    fusion = _HEAD22[:20] + _TAIL20  # junction after residue 20
    window = fusion[13:21]  # residues 14..21, straddles the junction
    table = _fusion_table(
        tmp_path,
        fusion_protein_sequence=fusion,
        fusion_protein_sequence_bp="20",
    )
    fasta = tmp_path / "shared.fa"
    fasta.write_text(f">shared_protein\nM{window}P\n", encoding="utf-8")
    index = build_proteome_kmer_index(read_fasta_sequences(fasta), (8,))
    result = build_easyfuse_catalog(
        table, "S1", "default", lengths=(8,),
        peptide_cfg=EasyFusePeptideConfig(
            normal_proteome_index=index, normal_proteome_source="unit-test",
        ),
    )
    hit = [r for r in result.catalog_rows if r["mutant_peptide"] == window]
    assert hit
    for r in hit:
        assert r["in_normal_proteome"] == "yes"
        assert r["contains_novel_aa"] == "no"
        assert r["junction_eligibility"] == "ELIGIBLE_JUNCTION_OVERLAP"
        assert r["crosses_junction"] == "yes"

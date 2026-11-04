"""Tests for the SV construction audit fixes (2026-11-04, items E1/E2).

E1: intronic DNA breakpoints now project onto bracketing exon boundaries so a
    spliced RNA junction can support the event (PROJECTED_SPLICED_JUNCTION).
E2 path A: a QC-pass RNA junction promotes a DNA-only ORF hypothesis to
    rna_supported_hypothesis and lets it reach peptide generation.
E2 path B: expressed products are classified honestly (same-gene products are
    junctions/frameshifts, not forced fusions).
"""

from pathlib import Path

from neoag.sv.exact_evidence import load_expressed_products
from neoag.sv.phase1 import build_sv_phase1_raw
from neoag.utils import read_tsv

ROOT = Path(__file__).resolve().parents[1]
FX = ROOT / "data" / "fixtures_sv"
FXI = ROOT / "data" / "fixtures_sv_intron"


def _run(tmp_path, *, junctions=True, products=None):
    return build_sv_phase1_raw(
        sample_id="SVINTRON",
        sv_vcfs=[FXI / "intron_del.vcf"],
        callers=["GRIDSS2"],
        tumor_sample_name="TUMOR",
        normal_sample_name="NORMAL",
        reference_fasta=FXI / "mini_ref.fa",
        gencode_gtf=FXI / "mini.gtf",
        hla=FX / "hla.txt",
        outdir=tmp_path,
        expression_tsv=FX / "expression.tsv",
        rna_junction_tsv=(FXI / "rna_junctions_projected.tsv") if junctions else None,
        expressed_products_tsv=products,
        normal_expression_tsv=FX / "normal_expression.tsv",
        normal_hla_ligands_tsv=FX / "normal_hla_ligands.tsv",
        peptide_lengths=(8, 9, 10, 11),
    )


def test_intronic_breakpoint_projected_junction_supports_event(tmp_path):
    out = _run(tmp_path)
    sv_events = read_tsv(out["sv_events_full"])
    assert len(sv_events) == 1
    assert sv_events[0]["rna_evidence_match"] == "PROJECTED_SPLICED_JUNCTION"
    assert sv_events[0]["rna_evidence_qc"] == "PASS"
    assert sv_events[0]["rna_support_status"] == "RNA_JUNCTION_SUPPORTED"
    assert sv_events[0]["effect_class"] == "SV_Frameshift"


def test_rna_supported_dna_hypothesis_reaches_peptides(tmp_path):
    out = _run(tmp_path)
    raw_events = read_tsv(out["raw_events"])
    raw_peptides = read_tsv(out["raw_peptides"])
    proteins = read_tsv(out["sv_protein_reconstruction"])
    # Path A was previously disabled for peptides; the RNA-validated
    # hypothesis must now produce candidates.
    assert len(raw_events) == 1
    assert raw_events[0]["event_type"] == "SV_Frameshift"
    assert len(raw_peptides) > 0
    assert all(p["contains_novel_aa"] == "yes" for p in raw_peptides)
    assert proteins[0]["reconstruction_method"].startswith("RNA_SUPPORTED_DNA_HYPOTHESIS:")
    assert proteins[0]["reconstruction_confidence"] == "medium"


def test_same_event_without_rna_evidence_stays_hypothesis(tmp_path):
    out = _run(tmp_path, junctions=False)
    assert read_tsv(out["raw_peptides"]) == []
    assert read_tsv(out["raw_events"]) == []
    sv_events = read_tsv(out["sv_events_full"])
    assert sv_events[0]["rna_evidence_match"] == "NO_EXACT_MATCH"
    assert sv_events[0]["reconstruction_status"] == "hypothesis_only"


def test_intragenic_expressed_product_is_not_forced_to_fusion(tmp_path):
    out = _run(tmp_path, products=FXI / "expressed_products_intragenic.tsv")
    proteins = read_tsv(out["sv_protein_reconstruction"])
    assert len(proteins) == 1
    # An out-of-frame same-gene product is a frameshift, and its whole novel
    # tail is peptide-eligible - not only the junction-crossing windows.
    assert proteins[0]["protein_type"] == "SV_Frameshift"
    assert read_tsv(out["raw_events"])[0]["event_type"] == "SV_Frameshift"


def _product_row(tmp_path, row: dict) -> str:
    cols = [
        "genome_build", "chrom1", "pos1", "strand1", "chrom2", "pos2", "strand2",
        "gene1", "gene2", "transcript1", "transcript2", "protein_sequence",
        "wildtype_protein_sequence", "junction_aa_position", "in_frame",
        "orf_status", "source_tool", "source_record_id", "effect_class",
    ]
    path = tmp_path / "products.tsv"
    values = [str(row.get(col, "")) for col in cols]
    path.write_text("\t".join(cols) + "\n" + "\t".join(values) + "\n", encoding="utf-8")
    return str(path)


def _base_product(tmp_path, **extra):
    row = {
        "genome_build": "GRCh38", "chrom1": "chr1", "pos1": 250, "strand1": ".",
        "chrom2": "chr1", "pos2": 500, "strand2": ".",
        "gene1": "GENE1", "gene2": "GENE1", "transcript1": "TX1", "transcript2": "TX1",
        "protein_sequence": "A" * 30 + "R" * 30, "wildtype_protein_sequence": "A" * 60,
        "junction_aa_position": 30, "in_frame": "no", "orf_status": "CONFIRMED",
        "source_tool": "FixTool", "source_record_id": "P1",
    }
    row.update(extra)
    return _product_row(tmp_path, row)


def test_expressed_product_effect_class_rules(tmp_path):
    derived_frameshift = next(iter(load_expressed_products(_base_product(tmp_path), default_build="GRCh38").values()))
    assert derived_frameshift.effect_class == "SV_Frameshift"

    derived_junction = next(iter(load_expressed_products(_base_product(tmp_path, in_frame="yes"), default_build="GRCh38").values()))
    assert derived_junction.effect_class == "SV_Junction"

    cross_gene = next(iter(load_expressed_products(_base_product(tmp_path, gene2="GENE2"), default_build="GRCh38").values()))
    assert cross_gene.effect_class == "SV_Fusion"

    explicit = next(iter(load_expressed_products(_base_product(tmp_path, effect_class="SV_Fusion"), default_build="GRCh38").values()))
    assert explicit.effect_class == "SV_Fusion"
    assert explicit.to_reconstruction("E1", "S1").protein_type == "SV_Fusion"

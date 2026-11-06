from __future__ import annotations

from pathlib import Path

import pytest

from neoag.sv.phase1 import build_sv_phase1_raw
from neoag.sv.score_pipeline import run_sv_score
from neoag.utils import read_tsv

ROOT = Path(__file__).resolve().parents[1]
FX = ROOT / "data" / "fixtures_sv"


def test_sv_score_binding_and_score_chain(tmp_path):
    adapter = tmp_path / "adapter"
    build_sv_phase1_raw(
        sample_id="SVMINI",
        sv_vcfs=[FX / "mini_sv.vcf"],
        callers=["GRIDSS2"],
        tumor_sample_name="TUMOR",
        normal_sample_name="NORMAL",
        reference_fasta=FX / "mini_ref.fa",
        gencode_gtf=FX / "mini.gtf",
        hla=FX / "hla.txt",
        outdir=adapter,
        expression_tsv=FX / "expression.tsv",
        rna_junction_tsv=FX / "rna_junctions_exact.tsv",
        expressed_products_tsv=FX / "expressed_products.tsv",
        normal_expression_tsv=FX / "normal_expression.tsv",
        normal_hla_ligands_tsv=FX / "normal_hla_ligands.tsv",
    )
    # The fixed stub binding asset cannot join SV candidate peptide-HLA keys,
    # which the core-predictor coverage validator (C1-08) rejects by design.
    # Provide run-matched binding evidence files instead of the stub copy so
    # this integration test exercises the real presentation join.
    cands = read_tsv(adapter / "parsed" / "raw_peptides.tsv")
    net_file = tmp_path / "net.xls"
    net_file.write_text(
        "Pos\tPeptide\tHLA\tScore_EL\t%Rank_EL\tScore_BA\t%Rank_BA\tBindLevel\tAlleleSupport\n"
        + "".join(
            f"1\t{p['peptide']}\t{p['hla_allele']}\t1.0\t1.0\t0.5\t2.0\tNB\ttrained\n"
            for p in cands
        ),
        encoding="utf-8",
    )
    mhc_file = tmp_path / "mhc.csv"
    mhc_file.write_text(
        "peptide,allele,mhcflurry_affinity,mhcflurry_affinity_percentile,mhcflurry_processing_score,mhcflurry_presentation_score\n"
        + "".join(f"{p['peptide']},{p['hla_allele']},50.0,1.5,0.6,0.7\n" for p in cands),
        encoding="utf-8",
    )
    out = run_sv_score(
        outdir=tmp_path / "scored",
        profile_name_or_path="sv_wgs_phase1",
        sample_id="SVMINI",
        raw_events=adapter / "parsed" / "raw_events.tsv",
        raw_peptides=adapter / "parsed" / "raw_peptides.tsv",
        netmhcpan=net_file,
        mhcflurry=mhc_file,
        expression=FX / "expression.tsv",
        normal_expression=FX / "normal_expression.tsv",
        normal_hla_ligands=FX / "normal_hla_ligands.tsv",
        immunogenicity_stub=True,
    )
    assert "ccf_2" in out
    assert Path(out["ccf_2"]).exists()
    assert Path(out["ccf_lite"]).exists()
    ccf_header = Path(out["ccf_2"]).read_text().splitlines()[0]
    assert "clonality_confidence" in ccf_header
    assert "ccf_resolution" in ccf_header
    peptides = read_tsv(out["ranked_peptides"])
    events = read_tsv(out["ranked_events"])
    assert len(events) >= 1
    assert len(peptides) >= 1
    assert peptides[0].get("efficacy_score")
    assert Path(out["evidence_report"]).is_file()

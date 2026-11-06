"""Tests for the 2026-11-06 coordinate-space fixes.

Fusion side (EasyFuse parser): peptides must be enumerated on one sequence
with that same sequence's own breakpoint, validated before generation;
inconsistent rows are flagged in the coordinate QC table, never guessed.

Splice side (event catalog): junction events without an upstream gene name
get their symbol from a coordinate-to-symbol overlap lookup when a GTF is
provided, with the lookup basis recorded in ``gene_symbol_source``.
"""

from pathlib import Path

from neoag.adapters.easyfuse_variant_peptide import (
    EasyFusePeptideConfig,
    build_easyfuse_catalog,
    easyfuse_row_to_variant_peptide_rows,
    select_fusion_peptide_context,
    write_easyfuse_qc_tables,
)
from neoag.adapters.event_catalog import (
    load_gene_interval_index,
    lookup_gene_symbol,
    parse_splice_catalog,
)
from neoag.splice.normalization import normalize_splice_sources

ROOT = Path(__file__).resolve().parents[1]


def _write_fusion_table(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.write_text(
        "\t".join(header) + "\n"
        + "\n".join("\t".join(row) for row in rows) + "\n",
        encoding="utf-8",
    )


FUSION_HEADER = [
    "BPID", "FTID", "Fusion_Gene", "Breakpoint1", "Breakpoint2", "type", "frame",
    "neo_peptide_sequence", "neo_peptide_sequence_bp",
    "fusion_protein_sequence", "fusion_protein_sequence_bp",
    "prediction_prob", "prediction_class", "ft_junc_cnt", "ft_anch_cnt",
]


def _fusion_row(*, neo="NEONEOPEP", neo_bp="5", protein="AAAAAAAABBBBBBBBCCCCCCCC", protein_bp="9", frame="in_frame"):
    return [
        "bp1", "FT1", "GENE1_GENE2", "chr1:1:+", "chr2:2:-", "trans", frame,
        neo, neo_bp, protein, protein_bp, "0.9", "positive", "10", "15",
    ]


def _catalog(tmp_path, rows, **cfg_kwargs):
    table = tmp_path / "fusions.tsv"
    _write_fusion_table(table, FUSION_HEADER, rows)
    cfg = EasyFusePeptideConfig(junction_only=True, dedup_per_event=True, **cfg_kwargs)
    return build_easyfuse_catalog(table, "S1", "default", lengths=(8, 9), peptide_cfg=cfg)


def test_fusion_context_prefers_protein_coordinate(tmp_path):
    # neo-space bp (5) and protein-space bp (9) point at different places;
    # windows must be built on the protein at its own bp, not on the protein
    # at the neo-space bp (the live-run defect).
    result = _catalog(tmp_path, [_fusion_row()])
    assert result.catalog_rows, "valid protein coordinates must yield peptides"
    for row in result.catalog_rows:
        assert row["peptide_context"] == "fusion_protein"
        assert row["peptide_context_bp"] == "9"
        assert row["peptide_context_status"] == "CONTEXT_FUSION_PROTEIN"
        start0 = int(row["peptide_start_aa"]) - 1
        assert row["junction_position_in_peptide_1based"] == str(9 - start0)
    assert result.summary_qc["peptide_context_fusion_protein_rows"] == "1"
    assert result.summary_qc["peptide_context_neo_fallback_rows"] == "0"
    assert result.summary_qc["peptide_context_invalid_rows"] == "0"


def test_fusion_context_neo_fallback_requires_valid_neo_bp(tmp_path):
    result = _catalog(tmp_path, [_fusion_row(protein_bp="", protein="")])
    assert result.catalog_rows
    for row in result.catalog_rows:
        assert row["peptide_context"] == "neo_peptide"
        assert row["peptide_context_status"] == "CONTEXT_NEO_FALLBACK"
        assert row["peptide_context_bp"] == "5"
    assert result.summary_qc["peptide_context_neo_fallback_rows"] == "1"
    # A neo bp outside the neo window is inconsistent -> no peptides, flagged.
    bad_dir = tmp_path / "bad"
    bad_dir.mkdir()
    bad = _catalog(bad_dir, [_fusion_row(neo_bp="42", protein_bp="", protein="")])
    assert bad.catalog_rows == []
    assert [qc["context_status"] for qc in bad.coordinate_qc] == ["CONTEXT_INVALID"]
    assert "bp=42" in bad.coordinate_qc[0]["context_detail"]
    assert bad.summary_qc["peptide_context_invalid_rows"] == "1"


def test_fusion_context_protein_bp_validated_against_protein_length(tmp_path):
    # protein bp beyond the protein length must not leak into the neo fallback
    # with an unrelated neo bp; the row is flagged instead.
    row = _fusion_row(protein_bp="99", neo_bp="42")
    result = _catalog(tmp_path, [row])
    assert result.catalog_rows == []
    qc = result.coordinate_qc[0]
    assert qc["context_status"] == "CONTEXT_INVALID"
    assert qc["junction_window_rows"] == "0"
    event_row = easyfuse_row_to_variant_peptide_rows(
        dict(zip(FUSION_HEADER, row)),
        {"event_id": "E1"},
        sample_id="S1",
        lengths=(8,),
        row_index=1,
    )
    assert event_row == []


def test_select_fusion_peptide_context_branches():
    fp = select_fusion_peptide_context({
        "neo_peptide_sequence": "ABCDEFGH", "neo_peptide_sequence_bp": "3",
        "fusion_protein_sequence": "A" * 40, "fusion_protein_sequence_bp": "35.6",
    })
    assert (fp["context"], fp["bp_pos"], fp["status"]) == ("fusion_protein", 35, "CONTEXT_FUSION_PROTEIN")
    neo = select_fusion_peptide_context({
        "neo_peptide_sequence": "ABCDEFGH", "neo_peptide_sequence_bp": "3",
    })
    assert (neo["context"], neo["bp_pos"], neo["status"]) == ("neo_peptide", 3, "CONTEXT_NEO_FALLBACK")
    bad = select_fusion_peptide_context({
        "neo_peptide_sequence": "ABCDEFGH", "neo_peptide_sequence_bp": "9",
    })
    assert bad["status"] == "CONTEXT_INVALID"
    assert bad["sequence"] == "" and bad["bp_pos"] == 0


def test_coordinate_qc_file_written(tmp_path):
    result = _catalog(tmp_path, [_fusion_row()])
    paths = write_easyfuse_qc_tables(result, tmp_path / "qc")
    qc_path = Path(paths["coordinate_qc"])
    assert qc_path.is_file()
    text = qc_path.read_text(encoding="utf-8")
    assert "CONTEXT_FUSION_PROTEIN" in text
    assert text.splitlines()[0].split("\t") == [
        "sample_id", "row_index", "event_id", "fusion_gene",
        "neo_peptide_length", "neo_peptide_sequence_bp",
        "fusion_protein_length", "fusion_protein_sequence_bp",
        "peptide_context", "peptide_context_bp", "context_status",
        "context_detail", "junction_window_rows",
    ]


GENE_GTF = "\n".join([
    'chr1\ttest\tgene\t1001\t5000\t.\t+\t.\tgene_id "ENSG1"; gene_name "COVERIN";',
    'chr1\ttest\tgene\t1001\t1500\t.\t+\t.\tgene_id "ENSG2"; gene_name "TINYGENE";',
    'chr1\ttest\tgene\t6001\t9000\t.\t-\t.\tgene_id "ENSG3"; gene_name "NEARGENE";',
]) + "\n"


def _write_regtools(path: Path, rows: list[list[str]]) -> None:
    header = ["chrom", "start", "end", "name", "score", "strand", "splice_site"]
    path.write_text(
        "\t".join(header) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n",
        encoding="utf-8",
    )


def test_gene_interval_lookup_prefers_strand_matched_containment(tmp_path):
    # fresh file per run so the (path, mtime) cache cannot serve a stale index
    gtf = tmp_path / "fresh_genes.gtf"
    gtf.write_text(GENE_GTF, encoding="utf-8")
    index = load_gene_interval_index(gtf)
    contained = lookup_gene_symbol(index, "chr1", 1200, 1400, "+")
    assert contained == ("TINYGENE", "contained")  # smaller containing gene wins
    assert lookup_gene_symbol(index, "chrX1", 1, 2, "+") == ("", "")
    assert lookup_gene_symbol(index, "1", 1200, 1400, "+") == contained  # chr prefix stripped
    overlap = lookup_gene_symbol(index, "chr1", 4800, 6500, "-")
    assert overlap[1] == "overlap"


def test_splice_catalog_coordinate_to_symbol_lookup(tmp_path):
    table = tmp_path / "junctions.tsv"
    _write_regtools(table, [
        ["chr1", "1199", "1400", "J_RESOLVED", "20", "+", "GT-AG"],   # intron inside both chr1 genes
        ["chr1", "9500", "9600", "J_EMPTY", "5", "+", "GT-AG"],       # overlaps nothing
    ])
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(GENE_GTF, encoding="utf-8")

    without = parse_splice_catalog(table, "S1", "default")
    by_name = {e["event_name"]: e for e in without}
    assert by_name["J_RESOLVED"]["gene"] == "chr1:1199"  # old honest fallback preserved
    assert by_name["J_RESOLVED"]["gene_symbol_source"] == "coordinate_fallback_no_reference"

    with_ref = parse_splice_catalog(table, "S1", "default", gene_gtf=gtf)
    by_name = {e["event_name"]: e for e in with_ref}
    assert by_name["J_RESOLVED"]["gene"] == "TINYGENE"
    assert by_name["J_RESOLVED"]["gene_symbol_source"] == "coordinate_lookup_contained"
    assert by_name["J_EMPTY"]["gene"] == "chr1:9500"  # still honest when nothing overlaps
    assert by_name["J_EMPTY"]["gene_symbol_source"] == "coordinate_lookup_no_overlap"


def test_splice_catalog_never_overrides_provided_gene(tmp_path):
    table = tmp_path / "named.tsv"
    table.write_text(
        "chrom\tstart\tend\tname\tscore\tstrand\tsplice_site\tgene_names\n"
        "chr1\t1199\t1400\tJ_NAMED\t9\t+\tGT-AG\tGIVENNAMED\n",
        encoding="utf-8",
    )
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(GENE_GTF, encoding="utf-8")
    events = parse_splice_catalog(table, "S1", "default", gene_gtf=gtf)
    assert events[0]["gene"] == "GIVENNAMED"
    assert events[0]["gene_symbol_source"] == "source_provided"


def _norm_events(outputs):
    from neoag.utils import read_tsv
    return read_tsv(outputs["raw_events"])


def test_normalize_splice_sources_names_events_by_containment(tmp_path):
    # Live STAR SJ route: intron 1200-1400 on chr1(+) sits fully inside both
    # TINYGENE (1001-1500) and COVERIN (1001-5000); the smaller must win.
    star = tmp_path / "SJ.out.tab"
    star.write_text("1\t1200\t1400\t2\t1\t2\t20\t18\t30\n", encoding="utf-8")
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(GENE_GTF, encoding="utf-8")

    without = normalize_splice_sources(
        sample_id="S1", junctions=star, junction_tool="STAR",
        junction_coordinate_system="star_sj", outdir=tmp_path / "n0",
    )
    ev0 = _norm_events(without)
    assert ev0 and ev0[0]["gene"].startswith("chr1:")  # honest coordinate fallback
    assert ev0[0]["gene_symbol_source"] == "coordinate_fallback_no_reference"

    with_gtf = normalize_splice_sources(
        sample_id="S1", junctions=star, junction_tool="STAR",
        junction_coordinate_system="star_sj", outdir=tmp_path / "n1",
        annotation_gtf=gtf,
    )
    ev1 = _norm_events(with_gtf)
    assert ev1[0]["gene"] == "TINYGENE"
    assert ev1[0]["gene_symbol_source"] == "coordinate_lookup_contained"
    # evidence is unchanged: only the label moves, reads stay on the junction
    assert ev0[0]["verified_rna_junction_reads"] == ev1[0]["verified_rna_junction_reads"]


def test_normalize_splice_sources_containment_only_never_guesses(tmp_path):
    # intron 4800-6500 only partially overlaps COVERIN and straddles NEARGENE:
    # containment-only must leave it on the coordinate fallback, never a name.
    star = tmp_path / "SJ.out.tab"
    star.write_text("1\t4800\t6500\t2\t1\t2\t30\t28\t40\n", encoding="utf-8")
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(GENE_GTF, encoding="utf-8")
    outputs = normalize_splice_sources(
        sample_id="S1", junctions=star, junction_tool="STAR",
        junction_coordinate_system="star_sj", outdir=tmp_path / "n2",
        annotation_gtf=gtf,
    )
    ev = _norm_events(outputs)
    assert ev[0]["gene"] == "chr1:4800"
    assert ev[0]["gene_symbol_source"] == "coordinate_lookup_no_containing_gene"

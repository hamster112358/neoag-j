"""Round-3 review B1-01 (FURTHER_SUGGESTIONS WS7): opt-in multitranscript VEP shadow mode.

Pinned properties:
- the single-transcript baseline is UNCHANGED: same rows, same header, no
  extra columns, no new summary keys (byte-identity against the pre-change
  code was additionally verified out-of-band on this fixture);
- multi mode enumerates every acceptable transcript, keeps distinct mutant
  peptides as separate hypotheses, and merges identical mutant peptides into
  the best-ranked transcript's row with transcript provenance;
- ALT -> CSQ -> transcript association is per-entry: FORMAT metrics are read
  at each entry's own allele index (multiallelic);
- shadow only: invalid modes raise, defaults never change output.
"""

from neoag.vep.extract_peptides import (
    MULTI_TRANSCRIPT_PROVENANCE_FIELD,
    OUTPUT_FIELDS,
    extract_variant_peptides_from_vcf,
)

CSQ_HEADER = (
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
    "Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|"
    "HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|"
    "Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|MANE_SELECT|"
    'TSL|HGVS_OFFSET|FrameshiftSequence|WildtypeProtein">'
)


def _entry(allele, cons, gene, ensg, tx, mane, canonical, pos, aa, wt, fs="", hgvsp=""):
    return "|".join([
        allele, cons, "MODERATE", gene, ensg, "Transcript", tx, "protein_coding",
        "1/2", "", "", hgvsp, "100", "9", pos, aa, "cTG/cTG", "", "", "1", "", "", "",
        canonical, mane, mane and "MANE" or "", "1", "", fs, wt,
    ])


def _write_vcf(path):
    # variant A: same protein change annotated on MANE ENST4 and non-MANE ENST5
    # variant B: multiallelic - ALT G on MANE ENST1, ALT T on canonical ENST2,
    #            plus an intron entry on ENST3 (not protein-altering)
    rec_a = (
        "chr2\t200\t.\tC\tT\t.\tPASS\tCSQ="
        + _entry("T", "missense_variant", "GENE2", "ENSG2", "ENST4", "MANE_Select", "YES",
                 "5", "A/G", "AAAAAAGAAAAA", hgvsp="p.A5G")
        + ","
        + _entry("T", "missense_variant", "GENE2", "ENSG2", "ENST5", "", "",
                 "5", "A/G", "AAAAAAGAAAAA", hgvsp="p.A5G")
        + "\tGT:AD:DP\t0/1:5,15:20\n"
    )
    rec_b = (
        "chr1\t100\t.\tA\tG,T\t.\tPASS\tCSQ="
        + _entry("G", "missense_variant", "GENE1", "ENSG1", "ENST1", "MANE_Select", "YES",
                 "3", "L/V", "MKLLVVVVVVVVVV")
        + ","
        + _entry("T", "missense_variant", "GENE1", "ENSG1", "ENST2", "", "YES",
                 "3", "L/I", "MKLLVVVVVVVVVV")
        + ","
        + ("T|intron_variant|LOW|GENE1|ENSG1|Transcript|ENST3|protein_coding"
           "||||||||||||||||||||")
        + "\tGT:AD:DP\t0/1:20,10,30:60\n"
    )
    path.write_text(
        "##fileformat=VCFv4.2\n" + CSQ_HEADER + "\n"
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n"
        + rec_b + rec_a,
        encoding="utf-8",
    )


def _extract(tmp_path, mode):
    vcf = tmp_path / "case.vcf"
    _write_vcf(vcf)
    out = tmp_path / f"peptides.{mode}.tsv"
    summary = extract_variant_peptides_from_vcf(
        vcf, out, lengths=(8,), sample_id="B1", transcript_selection=mode,
    )
    rows = []
    header = None
    for i, line in enumerate(out.read_text(encoding="utf-8").splitlines()):
        cols = line.split("\t")
        if i == 0:
            header = cols
        else:
            rows.append(dict(zip(header, cols)))
    return summary, header, rows


def test_b1_01_single_mode_baseline_header_and_transcripts_unchanged(tmp_path):
    summary, header, rows = _extract(tmp_path, "single")
    assert header == list(OUTPUT_FIELDS)  # no new columns leak into the baseline
    assert MULTI_TRANSCRIPT_PROVENANCE_FIELD not in header
    # B1-01 shadow keys never appear in single-mode summaries
    for key in ("transcript_selection", "variants_with_alternate_transcripts",
                "multitranscript_deduplicated_peptides"):
        assert key not in summary
    # every variant contributed only its best-ranked transcript
    assert {r["transcript_id"] for r in rows} == {"ENST1", "ENST4"}
    assert "transcript_selection=single" not in summary["generation_method"]


def test_b1_01_multi_mode_keeps_distinct_peptides_as_separate_hypotheses(tmp_path):
    _summary, header, rows = _extract(tmp_path, "multi")
    assert header == list(OUTPUT_FIELDS) + [MULTI_TRANSCRIPT_PROVENANCE_FIELD]
    b_rows = [r for r in rows if r["gene"] == "GENE1"]
    tx = {r["transcript_id"] for r in b_rows}
    assert tx == {"ENST1", "ENST2"}  # intron entry ENST3 never contributes
    # each row's provenance names its own transcript first
    for r in b_rows:
        assert r[MULTI_TRANSCRIPT_PROVENANCE_FIELD].startswith(f"TX={r['transcript_id']}|")
    # the two transcripts produce disjoint peptide sets (different mutant proteins)
    p1 = {r["mutant_peptide"] for r in b_rows if r["transcript_id"] == "ENST1"}
    p2 = {r["mutant_peptide"] for r in b_rows if r["transcript_id"] == "ENST2"}
    assert p1 and p2 and not (p1 & p2)


def test_b1_01_multi_mode_dedups_identical_peptides_with_transcript_provenance(tmp_path):
    summary, _header, rows = _extract(tmp_path, "multi")
    a_rows = [r for r in rows if r["gene"] == "GENE2"]
    # identical mutation on ENST4 and ENST5: one row per peptide, not two
    assert len(a_rows) == len({r["mutant_peptide"] for r in a_rows})
    assert all(r["transcript_id"] == "ENST4" for r in a_rows)  # best rank kept
    for r in a_rows:
        prov = r[MULTI_TRANSCRIPT_PROVENANCE_FIELD]
        assert "TX=ENST4|HGVSp=p.A5G" in prov and "TX=ENST5|HGVSp=p.A5G" in prov
    assert summary["multitranscript_deduplicated_peptides"] == len(a_rows)
    assert summary["variants_with_alternate_transcripts"] == 2
    assert "transcript_selection=multi" in summary["generation_method"]


def test_b1_01_multi_mode_reads_metrics_at_each_entrys_own_allele_index(tmp_path):
    _summary, _header, rows = _extract(tmp_path, "multi")
    b1 = {r["transcript_id"]: r for r in rows if r["gene"] == "GENE1"}
    # AD = 20(ref),10(ALT G),30(ALT T): ENST1 annotates ALT G, ENST2 annotates ALT T
    assert b1["ENST1"]["tumor_alt_count"] == "10"
    assert b1["ENST2"]["tumor_alt_count"] == "30"
    assert b1["ENST1"]["vaf"] == "0.1667"
    assert b1["ENST2"]["vaf"] == "0.5000"
    assert b1["ENST1"]["tumor_depth"] == b1["ENST2"]["tumor_depth"] == "60"


def test_b1_01_multi_mode_is_superset_of_baseline_for_same_vcf(tmp_path):
    _s1, _h1, single_rows = _extract(tmp_path, "single")
    _s2, _h2, multi_rows = _extract(tmp_path, "multi")
    single_keys = {(r["transcript_id"], r["mutant_peptide"]) for r in single_rows}
    multi_keys = {(r["transcript_id"], r["mutant_peptide"]) for r in multi_rows}
    assert single_keys < multi_keys  # shadow mode only ADDS hypotheses
    single_peps = {r["mutant_peptide"] for r in single_rows}
    multi_peps = {r["mutant_peptide"] for r in multi_rows}
    assert single_peps < multi_peps  # recall gain is visible at peptide level


def test_b1_01_invalid_selection_raises(tmp_path):
    import pytest

    vcf = tmp_path / "case.vcf"
    _write_vcf(vcf)
    with pytest.raises(ValueError, match="unknown transcript_selection"):
        extract_variant_peptides_from_vcf(
            vcf, tmp_path / "x.tsv", transcript_selection="both",
        )

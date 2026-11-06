"""Package G-a tests: runtime config, sample roles, genome-build keys.

Review items locked down:

B1-02  the runtime TOML is emitted from typed data and round-trip
       verified; a quote or backslash inside a path can no longer change
       the document's meaning.
B1-03  a multi-sample VCF whose tumor/normal roles cannot be resolved
       fails closed instead of falling through to column-order guessing.
B5-03  tool-consensus splice keys use the declared run build instead of a
       hardcoded GRCh38 label.
C6-06  cross-site RNA event keys embed each row's genome build (or the
       declared run build), so same-coordinate events from different
       builds do not silently compare.
"""

from __future__ import annotations

import importlib.util
import tomllib
from pathlib import Path

import pytest

from neoag.cross_site_rna import build_cross_site_rna_evidence, event_keys
from neoag.open_neo.tool_consensus import _write_splice
from neoag.utils import read_tsv, write_toml

ROOT = Path(__file__).resolve().parents[1]


def _script_module():
    path = ROOT / "scripts" / "run_candidate_upstream.py"
    spec = importlib.util.spec_from_file_location("run_candidate_upstream", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# --------------------------------------------------------------------------- B1-02


def test_write_toml_roundtrips_hostile_paths(tmp_path):
    nasty = "/data/case's \"quoted\"\\dir/ref.fa"
    cfg = {
        "sample": {"id": "S \"1\"", "profile": "default"},
        "tools": {"enabled": []},
        "inputs": {
            "entry_mode": "snv_indel",
            "variants_vcf": nasty,
            "variant_peptide_extraction": True,
            "hla_alleles": ["HLA-A*02:01", "HLA-B*07:02"],
            "reference_fasta": "/data/space dir/ref.fa",
            "tumor_sample_name": "TUMOR",
            "normal_sample_name": "",
        },
    }
    path = write_toml(tmp_path / "run.toml", cfg)
    parsed = tomllib.loads(path.read_text(encoding="utf-8"))
    assert parsed["inputs"]["variants_vcf"] == nasty
    assert parsed["sample"]["id"] == "S \"1\""
    assert parsed["inputs"]["variant_peptide_extraction"] is True
    assert parsed["inputs"]["hla_alleles"] == ["HLA-A*02:01", "HLA-B*07:02"]
    assert parsed["inputs"]["normal_sample_name"] == ""


def test_write_toml_omits_none_values(tmp_path):
    path = write_toml(tmp_path / "run.toml", {"inputs": {"a": 1, "b": None}})
    parsed = tomllib.loads(path.read_text(encoding="utf-8"))
    assert parsed == {"inputs": {"a": 1}}


def test_upstream_script_uses_typed_toml(tmp_path):
    # The snv branch previously built the config by string interpolation;
    # pin that the script now imports and calls the typed writer.
    source = (ROOT / "scripts" / "run_candidate_upstream.py").read_text(encoding="utf-8")
    assert "write_toml(" in source
    assert 'f\'id = "{args.sample_id}"\'' not in source


# --------------------------------------------------------------------------- B1-03


def _vcf(tmp_path, samples: list[str]) -> Path:
    path = tmp_path / f"case_{len(samples)}.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(samples)
        + "\n"
        + "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t"
        + "\t".join("./1" for _ in samples)
        + "\n",
        encoding="utf-8",
    )
    return path


def test_ambiguous_two_sample_roles_fail_closed(tmp_path):
    module = _script_module()
    vcf = _vcf(tmp_path, ["SAMPLE1", "SAMPLE2"])
    with pytest.raises(SystemExit, match="cannot resolve tumor/normal roles"):
        module.resolve_vcf_sample_roles(vcf)


def test_explicit_roles_follow_names_not_column_order(tmp_path):
    module = _script_module()
    # Swapped column order: NORMAL first.  Explicit names must win.
    vcf = _vcf(tmp_path, ["NORMAL", "TUMOR"])
    tumor, normal = module.resolve_vcf_sample_roles(vcf)
    assert (tumor, normal) == ("TUMOR", "NORMAL")
    vcf2 = _vcf(tmp_path, ["S1", "S2"])
    tumor, normal = module.resolve_vcf_sample_roles(vcf2, tumor_sample_name="S2", normal_sample_name="S1")
    assert (tumor, normal) == ("S2", "S1")


def test_explicit_single_role_completes_pair(tmp_path):
    module = _script_module()
    vcf = _vcf(tmp_path, ["ALPHA", "BETA"])
    tumor, normal = module.resolve_vcf_sample_roles(vcf, normal_sample_name="ALPHA")
    assert (tumor, normal) == ("BETA", "ALPHA")


def test_single_sample_vcf_is_the_tumor(tmp_path):
    module = _script_module()
    vcf = _vcf(tmp_path, ["ONLY"])
    tumor, normal = module.resolve_vcf_sample_roles(vcf)
    assert tumor == "ONLY"
    assert normal == ""


# --------------------------------------------------------------------------- B5-03


def _rna_tsv(tmp_path) -> Path:
    path = tmp_path / "rna.tsv"
    path.write_text(
        "chrom\tstart\tend\tname\tscore\tstrand\tgene_names\n"
        "chr1\t99\t200\tJ1\t10\t+\tGENE1\n",
        encoding="utf-8",
    )
    return path


def _declared(rna: Path) -> dict:
    return {"splice_dna": {}, "splice_rna": {"regtools": str(rna)}, "splice_neoantigen": {}}


def test_splice_consensus_uses_declared_build(tmp_path):
    rows = _write_splice(_declared(_rna_tsv(tmp_path)), tmp_path / "c37", genome_build="GRCh37")[1]
    canonical = [row["canonical_junction_id"] for row in rows if row["canonical_junction_id"].startswith("SJ|")]
    assert canonical and all(cid.startswith("SJ|GRCh37|") for cid in canonical)


def test_splice_consensus_default_build_unchanged(tmp_path):
    rows = _write_splice(_declared(_rna_tsv(tmp_path)), tmp_path / "c38")[1]
    canonical = [row["canonical_junction_id"] for row in rows if row["canonical_junction_id"].startswith("SJ|")]
    assert canonical and all(cid.startswith("SJ|GRCh38|") for cid in canonical)


def test_row_level_build_column_overrides_default(tmp_path):
    rna = tmp_path / "rna37.tsv"
    rna.write_text(
        "chrom\tstart\tend\tname\tscore\tstrand\tgene_names\tgenome_build\n"
        "chr1\t99\t200\tJ1\t10\t+\tGENE1\tGRCh37\n",
        encoding="utf-8",
    )
    rows = _write_splice(_declared(rna), tmp_path / "c38row", genome_build="GRCh38")[1]
    canonical = [row["canonical_junction_id"] for row in rows if row["canonical_junction_id"].startswith("SJ|")]
    assert canonical and all(cid.startswith("SJ|GRCh37|") for cid in canonical)


# --------------------------------------------------------------------------- C6-06


def _event(event_id: str, *, build: str = "") -> dict:
    return {
        "event_id": event_id,
        "event_type": "SNV",
        "chrom": "chr1",
        "pos": "100",
        "ref": "A",
        "alt": "T",
        "genome_build": build,
        "rna_alt_reads": "5",
        "rna_depth": "40",
    }


def test_snv_keys_use_row_build_not_hardcoded_label():
    (specificity, method, key) = event_keys(_event("E1", build="GRCh37"))[0]
    assert key == "GRCh37|chr1|100|A|T"
    (specificity, method, key) = event_keys(_event("E2"), default_build="GRCh37")[0]
    assert key == "GRCh37|chr1|100|A|T"
    (specificity, method, key) = event_keys(_event("E3"))[0]
    assert key == "GRCh38|chr1|100|A|T"


def test_cross_build_events_do_not_match(tmp_path):
    def write_events(path: Path, rows: list[dict]) -> None:
        cols = sorted({key for row in rows for key in row})
        path.write_text(
            "\t".join(cols) + "\n"
            + "".join("\t".join(row.get(col, "") for col in cols) + "\n" for row in rows),
            encoding="utf-8",
        )

    primary = tmp_path / "primary.tsv"
    secondary = tmp_path / "secondary.tsv"
    write_events(primary, [_event("P1", build="GRCh38")])
    # Same coordinates, different build in the secondary table: without
    # liftOver provenance this must NOT be counted as the same event.
    write_events(secondary, [_event("S1", build="GRCh37")])
    out = tmp_path / "cross.tsv"
    summary = build_cross_site_rna_evidence(primary, secondary, out, secondary_sample_id="S", genome_build="GRCh38")
    rows = read_tsv(out)
    assert rows[0]["cross_site_status"] == "PRIMARY_ONLY"
    assert summary["default_genome_build"] == "GRCh38"
    # And when both carry the same declared build they do match.
    write_events(secondary, [_event("S1", build="GRCh38")])
    summary2 = build_cross_site_rna_evidence(primary, secondary, out, secondary_sample_id="S", identity_status="CONFIRMED")
    assert read_tsv(out)[0]["cross_site_status"] == "EXACT_SHARED"

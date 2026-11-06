from __future__ import annotations
from pathlib import Path
from ..utils import MISSING, first, read_csv, read_tsv, safe_id, to_float
from ..evidence_provenance import ProvenanceRecord, provenance_from_file, without_provenance, write_evidence_tsv
from ..schemas import MHCFLURRY_EVIDENCE_FIELDS


def _number_or_blank(row: dict[str, str], keys: list[str]) -> str:
    """Review C1-09: a missing or unreadable metric is "" (typed missingness).

    The old defaults turned absence into percentile 99 / score 0 - i.e. a
    prediction that was never made was reported as a failed prediction. A
    measured 0 still parses, so legitimate zero measurements are preserved.
    """
    raw = first(row, keys, "")
    text = str(raw).strip()
    if not text or text in MISSING:
        return ""
    value = to_float(text, float("nan"))
    if value != value:  # NaN check: unparseable input is absence, not 0
        return ""
    return str(value)


def parse_mhcflurry(path: str | Path, sample_id: str = "") -> list[dict[str, str]]:
    p = Path(path)
    rows_in = read_csv(p) if p.suffix.lower() == ".csv" else read_tsv(p)
    out = []
    for r in rows_in:
        pep = first(r, ["peptide", "Peptide", "sequence"], "")
        allele = first(r, ["allele", "Allele", "hla", "HLA"], "")
        if pep and allele:
            out.append({
                "sample_id": sample_id,
                "peptide": pep,
                "hla_allele": allele,
                "peptide_hla_key": safe_id(f"{pep}_{allele}"),
                "mhcflurry_affinity": _number_or_blank(r, ["mhcflurry_affinity", "affinity", "prediction", "ic50"]),
                "mhcflurry_affinity_percentile": _number_or_blank(r, ["mhcflurry_affinity_percentile", "affinity_percentile", "percentile_rank", "percentile"]),
                "mhcflurry_processing_score": _number_or_blank(r, ["mhcflurry_processing_score", "processing_score"]),
                "mhcflurry_presentation_score": _number_or_blank(r, ["mhcflurry_presentation_score", "presentation_score"]),
                # Review C1-11: MHCflurry's per-allele training coverage is
                # not reported in its output, so every row records the
                # coverage as unassessed rather than implying full support.
                "mhcflurry_allele_support_status": "MODEL_COVERAGE_UNASSESSED",
                "source_file": str(p),
            })
    return out

def write_mhcflurry_evidence(path, rows, provenance: ProvenanceRecord | None = None):
    src = rows[0].get("source_file") if rows else path
    prov = provenance or provenance_from_file("mhcflurry", src, mode="passthrough")
    write_evidence_tsv(path, rows, without_provenance(MHCFLURRY_EVIDENCE_FIELDS), prov)

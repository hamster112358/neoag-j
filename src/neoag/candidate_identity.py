"""Stable candidate identities used by ranking and report de-duplication."""

from __future__ import annotations

import hashlib
import re
from typing import Any, Mapping

from .status_words import first_text


IDENTITY_FIELDS = (
    "event_identity_id",
    "protein_change_identity_id",
    "peptide_sequence_id",
    "peptide_hla_id",
)


def _text(row: Mapping[str, Any], *fields: str) -> str:
    # Review D4-05: the null vocabulary ("NA", "NONE", ".", "UNASSESSED") is
    # centralized in status_words.first_text so every identity layer agrees on
    # which values mean "no information".
    return first_text(row, *fields)


def _stable_id(prefix: str, *parts: str) -> str:
    payload = "|".join(part.strip() for part in parts)
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:24]
    return f"{prefix}_{digest}"


def normalize_peptide(value: Any) -> str:
    return re.sub(r"\s+", "", str(value or "")).upper()


def normalize_hla(value: Any) -> str:
    text = str(value or "").strip().upper().replace("_", "-")
    text = re.sub(r"^HLA-?", "", text)
    match = re.match(r"^([A-Z0-9]+)\*?([0-9]{2,3}):?([0-9]{2,3})(?::?([0-9]{2,3}))?", text)
    if not match:
        return text
    locus, first, second, third = match.groups()
    suffix = f":{third}" if third else ""
    return f"HLA-{locus}*{first}:{second}{suffix}"


def candidate_identity(row: Mapping[str, Any]) -> dict[str, str]:
    """Return event, protein-change, peptide and peptide-HLA identities.

    The event identity remains event-specific. A shared peptide-HLA identity can
    therefore collapse duplicate ranking positions while preserving links to
    distinct breakpoints, transcript hypotheses or ORFs.

    event_identity_status states how the event identity was derived
    (review D4-04) so consumers can honour an UNRESOLVED_ROW_FALLBACK
    without string-sniffing the id prefix:
      EXACT_PRIMARY_ID       - a caller-supplied event id was present
      COORDINATE_DERIVED_ID  - synthesised from event_type/gene/coordinates
      UNRESOLVED_ROW_FALLBACK- row digest; NOT a biological event identity
    """
    event_source = _text(
        row,
        "canonical_event_id",
        "event_group_id",
        "event_id",
        "source_event_id",
        "event_name",
    )
    if event_source:
        identity_status = "EXACT_PRIMARY_ID"
    else:
        coordinate_parts = [
            _text(row, field)
            for field in ("event_type", "gene", "chrom", "pos", "ref", "alt")
        ]
        if any(coordinate_parts):
            event_source = "|".join(coordinate_parts)
            identity_status = "COORDINATE_DERIVED_ID"
        else:
            # Review D4-04: with no primary id AND no coordinates the old code
            # hashed the literal string "|||||", so every unrelated unknown row
            # collided on one event id and lost identity in deduplication. The
            # fallback is now row-specific (digest of the full raw row) and
            # explicitly marked UNRESOLVED so it cannot be mistaken for a
            # biological event identity.
            payload = "|".join(
                f"{key}={row.get(key, '')}" for key in sorted(row, key=str)
            )
            event_source = "UNRESOLVED_EVENT_ROW:" + hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]
            identity_status = "UNRESOLVED_ROW_FALLBACK"
    event_id = _stable_id("EVT", event_source)

    protein_change = _text(
        row,
        "protein_change_id",
        "orf_id",
        "transcript_hypothesis_id",
        "combined_protein_change",
        "protein_change",
        "hgvsp",
        "amino_acid_change",
    )
    protein_id = _stable_id("PROT", event_id, protein_change or "UNRESOLVED_PROTEIN_CHANGE")

    peptide = normalize_peptide(_text(row, "peptide", "mutant_peptide", "mt_peptide"))
    peptide_id = _stable_id("PEP", peptide or "UNRESOLVED_PEPTIDE")
    hla = normalize_hla(_text(row, "hla_allele", "hla", "allele", "restricting_hla"))
    peptide_hla_id = _stable_id("PHLA", peptide or "UNRESOLVED_PEPTIDE", hla or "UNRESOLVED_HLA")
    return {
        "event_identity_id": event_id,
        "protein_change_identity_id": protein_id,
        "peptide_sequence_id": peptide_id,
        "peptide_hla_id": peptide_hla_id,
        "event_identity_status": identity_status,
    }


def identity_value(row: Mapping[str, Any], field: str) -> str:
    value = str(row.get(field) or "").strip()
    return value or candidate_identity(row)[field]

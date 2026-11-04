"""Strict RNA-junction and expressed-product inputs for the SV pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ..utils import first, read_tsv, to_float
from .identity import canonical_breakpoint_key, normalize_contig
from .protein_reconstruct import ProteinReconstruction
from .reference import TranscriptModel


def _required(row: dict[str, str], names: list[str], label: str) -> str:
    value = first(row, names, "").strip()
    if not value:
        raise ValueError(f"Exact SV evidence row is missing required column/value: {label}")
    return value


def _key(row: dict[str, str], default_build: str) -> str:
    build = first(row, ["genome_build", "build", "assembly"], default_build).strip()
    # Strand is optional evidence: canonical_breakpoint_key() normalises any
    # token outside {+,-} to "." anyway, and DEL/END-style DNA records carry
    # no strand at all, so requiring a non-"." value here previously made it
    # impossible to ever join a product/junction table to a DEL adjacency.
    strand1 = str(row.get("strand1") or row.get("orientation1") or ".").strip()
    strand2 = str(row.get("strand2") or row.get("orientation2") or ".").strip()
    return canonical_breakpoint_key(
        build,
        _required(row, ["chrom1", "chr1"], "chrom1"),
        int(float(_required(row, ["pos1", "breakpoint1"], "pos1"))),
        strand1,
        _required(row, ["chrom2", "chr2"], "chrom2"),
        int(float(_required(row, ["pos2", "breakpoint2"], "pos2"))),
        strand2,
    )


@dataclass(frozen=True)
class ExactJunctionEvidence:
    adjacency_key: str
    split_reads: int
    spanning_reads: int
    unique_start_count: int
    min_anchor_bp: int
    min_mapq: int
    source_tool: str
    source_record_id: str

    @property
    def support_reads(self) -> int:
        # Split reads establish the exact junction; spanning reads are auxiliary.
        return self.split_reads

    @property
    def qc_pass(self) -> bool:
        return (
            self.split_reads >= 3
            and self.unique_start_count >= 2
            and self.min_anchor_bp >= 10
            and self.min_mapq >= 20
        )


def load_exact_junction_evidence(
    path: str | Path | None,
    *,
    default_build: str,
) -> dict[str, ExactJunctionEvidence]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Missing RNA junction evidence: {p}")
    out: dict[str, ExactJunctionEvidence] = {}
    for row in read_tsv(p):
        key = _key(row, default_build)
        evidence = ExactJunctionEvidence(
            adjacency_key=key,
            split_reads=int(to_float(first(row, ["split_reads", "junction_reads"], "0"), 0)),
            spanning_reads=int(to_float(first(row, ["spanning_reads"], "0"), 0)),
            unique_start_count=int(to_float(first(row, ["unique_start_count", "unique_starts"], "0"), 0)),
            min_anchor_bp=int(to_float(first(row, ["min_anchor_bp", "anchor_bp"], "0"), 0)),
            min_mapq=int(to_float(first(row, ["min_mapq", "mapq"], "0"), 0)),
            source_tool=first(row, ["source_tool", "caller"], ""),
            source_record_id=first(row, ["source_record_id", "record_id"], ""),
        )
        previous = out.get(key)
        if previous is None or (evidence.qc_pass, evidence.support_reads) > (previous.qc_pass, previous.support_reads):
            out[key] = evidence
    return out


# ---------------------------------------------------------------------------
# Projected exon-boundary matching (audit item SV E1).
#
# The strict tier joins a DNA adjacency to an RNA junction only when the
# breakpoint coordinates are byte-identical.  That holds for precise BNDs but
# fails for the common case the auditor flagged: an intragenic deletion whose
# DNA breakpoints sit in introns while spliced RNA reports the junction at the
# flanking exon boundaries (1-based: last base of the retained upstream exon,
# first base of the retained downstream exon).  The projected tier keeps exact
# coordinate identity - it only moves the lookup to annotated exon boundaries
# that bracket the DNA breakpoints - so a match still requires the RNA table
# to report that precise adjacency, and it never fires for a position that is
# itself a breakpoint (the strict tier already covers those).
# ---------------------------------------------------------------------------

PositionalKey = tuple[str, str, int, str, int]


def _canonical_pair(chrom1: str, pos1: int, chrom2: str, pos2: int) -> tuple[tuple[str, int], tuple[str, int]]:
    left = (normalize_contig(chrom1), int(pos1))
    right = (normalize_contig(chrom2), int(pos2))
    return tuple(sorted((left, right)))  # type: ignore[return-value]


def junction_positional_index(
    evidence: dict[str, ExactJunctionEvidence],
) -> dict[PositionalKey, ExactJunctionEvidence]:
    """Index junction evidence by (build, chrom, pos, chrom, pos), strand-insensitive.

    Coordinates come from the canonical key, so they are already contig- and
    position-sorted the same way _canonical_pair sorts them.
    """

    out: dict[PositionalKey, ExactJunctionEvidence] = {}
    for key, item in evidence.items():
        build, first_token, second_token = key.split("|")
        c1, p1, _s1 = first_token.split(":")
        c2, p2, _s2 = second_token.split(":")
        pos_key: PositionalKey = (build, c1, int(p1), c2, int(p2))
        previous = out.get(pos_key)
        if previous is None or (item.qc_pass, item.support_reads) > (previous.qc_pass, previous.support_reads):
            out[pos_key] = item
    return out


def _lookup_positional(
    positional: dict[PositionalKey, ExactJunctionEvidence],
    build: str,
    chrom1: str,
    pos1: int,
    chrom2: str,
    pos2: int,
) -> ExactJunctionEvidence | None:
    (c1, p1), (c2, p2) = _canonical_pair(chrom1, pos1, chrom2, pos2)
    hit = positional.get((build.strip().lower(), c1, p1, c2, p2))
    return hit if hit and hit.qc_pass else None


def _bracketing_exon_boundaries(tx: TranscriptModel, pos_lo: int, pos_hi: int) -> tuple[int, int] | None:
    """For an intragenic deletion: last exon end below the event, first exon start above it."""

    ends = [f.end for f in tx.exons if f.end < pos_lo]
    starts = [f.start for f in tx.exons if f.start > pos_hi]
    if not ends or not starts:
        return None
    return max(ends), min(starts)


def _retained_boundary(tx: TranscriptModel, pos: int, side: str) -> int | None:
    """Nearest annotated exon boundary on the retained side of a breakend.

    side="donor" is the boundary ending the transcript prefix kept on the 5'
    side of the breakpoint; side="acceptor" starts the retained 3' suffix.
    Transcript-strand aware (on '-' strands the 5' side is at higher
    coordinates).
    """

    if side == "donor":
        sides = (
            ([f.end for f in tx.exons if f.end <= pos], max)
            if tx.strand == "+"
            else ([f.start for f in tx.exons if f.start >= pos], min)
        )
    else:
        sides = (
            ([f.start for f in tx.exons if f.start >= pos], min)
            if tx.strand == "+"
            else ([f.end for f in tx.exons if f.end <= pos], max)
        )
    values, pick = sides
    return pick(values) if values else None


def match_projected_junction(
    positional: dict[PositionalKey, ExactJunctionEvidence],
    *,
    build: str,
    chrom1: str,
    pos1: int,
    chrom2: str,
    pos2: int,
    svtype: str,
    same_gene: bool,
    tx1: TranscriptModel | None,
    tx2: TranscriptModel | None,
) -> tuple[str, ExactJunctionEvidence] | None:
    """Return (match_kind, evidence) when an exon-boundary projection matches.

    Scope is deliberately narrow and honest:
      - same-gene same-chrom DEL: one candidate junction bracketing the
        deleted interval (both boundaries strictly outside it, so the wild-type
        splice junction of a neighbouring exon pair can never be hit);
      - different genes: the two pairings the DNA model itself assumes
        (partner1-donor x partner2-acceptor, and the reverse).
    Everything else (DUP/INV/unknown gene context) stays unmatched rather than
    being fuzzy-matched.
    """

    if not positional:
        return None
    same_chrom = normalize_contig(chrom1) == normalize_contig(chrom2)
    if same_gene and same_chrom:
        if svtype.upper() != "DEL":
            return None
        tx = tx1 or tx2
        if tx is None:
            return None
        bounds = _bracketing_exon_boundaries(tx, min(pos1, pos2), max(pos1, pos2))
        if bounds is None:
            return None
        hit = _lookup_positional(positional, build, chrom1, bounds[0], chrom2, bounds[1])
        return ("PROJECTED_SPLICED_JUNCTION", hit) if hit else None
    if tx1 is not None and tx2 is not None:
        for donor_tx, donor_pos, acceptor_tx, acceptor_pos, kind in (
            (tx1, pos1, tx2, pos2, "PROJECTED_PARTNER_JUNCTION"),
            (tx2, pos2, tx1, pos1, "PROJECTED_PARTNER_JUNCTION_REV"),
        ):
            donor = _retained_boundary(donor_tx, donor_pos, "donor")
            acceptor = _retained_boundary(acceptor_tx, acceptor_pos, "acceptor")
            if donor is None or acceptor is None:
                continue
            hit = _lookup_positional(positional, build, donor_tx.chrom, donor, acceptor_tx.chrom, acceptor)
            if hit:
                return kind, hit
    return None


@dataclass(frozen=True)
class ExpressedProduct:
    adjacency_key: str
    gene1: str
    gene2: str
    transcript1: str
    transcript2: str
    protein_sequence: str
    wildtype_protein_sequence: str
    junction_aa_position: int
    in_frame: str
    orf_status: str
    source_tool: str
    source_record_id: str
    effect_class: str = "SV_Fusion"

    def to_reconstruction(self, event_id: str, sample_id: str) -> ProteinReconstruction:
        sequence = self.protein_sequence.upper().replace("*", "")
        junction = self.junction_aa_position
        return ProteinReconstruction(
            protein_sequence_id=f"PROT_{event_id}",
            event_id=event_id,
            sample_id=sample_id,
            gene=(f"{self.gene1}::{self.gene2}" if self.effect_class == "SV_Fusion" and self.gene2 != self.gene1 else self.gene1),
            transcript_id=f"{self.transcript1}::{self.transcript2}",
            protein_type=self.effect_class,
            protein_sequence=sequence,
            wt_protein_sequence=self.wildtype_protein_sequence.upper().replace("*", ""),
            wt_prefix_aa=self.wildtype_protein_sequence.upper().replace("*", "")[:junction],
            novel_aa=sequence[junction : junction + 80],
            junction_aa_position=str(junction),
            novel_start_aa=str(junction),
            frameshift_start_aa="" if self.in_frame == "yes" else str(junction),
            in_frame=self.in_frame,
            reconstruction_method=f"external_expressed_transcript:{self.source_tool or 'unspecified'}",
            reconstruction_confidence="high",
            reconstruction_reason=(
                f"exact_adjacency={self.adjacency_key};orf_status={self.orf_status};"
                f"source_record_id={self.source_record_id}"
            ),
        )


def load_expressed_products(
    path: str | Path | None,
    *,
    default_build: str,
) -> dict[str, ExpressedProduct]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Missing expressed-product table: {p}")
    out: dict[str, ExpressedProduct] = {}
    for row in read_tsv(p):
        key = _key(row, default_build)
        sequence = _required(row, ["protein_sequence", "fusion_protein"], "protein_sequence").upper().replace("*", "")
        junction = int(float(_required(row, ["junction_aa_position", "junction_aa"], "junction_aa_position")))
        if not (0 < junction < len(sequence)):
            raise ValueError(f"Invalid junction_aa_position={junction} for product {key}")
        orf_status = _required(row, ["orf_status", "translation_status"], "orf_status").upper()
        if orf_status not in {"CONFIRMED", "COMPLETE", "TRANSLATABLE"}:
            raise ValueError(f"Expressed product {key} is not a confirmed translatable ORF: {orf_status}")
        # Audit item SV E2 (path B): do not force every product to SV_Fusion.
        # Honour an explicit type column; otherwise derive it honestly from the
        # product itself - same-gene products are junctions/frameshifts, only
        # cross-gene products are fusions.  Fusion-classified products are the
        # only ones peptide_builder will restrict to junction-crossing windows.
        raw_type = first(row, ["effect_class", "protein_type", "consequence_class"], "").strip()
        type_alias = {
            "fusion": "SV_Fusion", "sv_fusion": "SV_Fusion",
            "frameshift": "SV_Frameshift", "sv_frameshift": "SV_Frameshift",
            "insertion": "SV_Insertion", "sv_insertion": "SV_Insertion",
            "junction": "SV_Junction", "sv_junction": "SV_Junction",
        }
        in_frame = "yes" if first(row, ["in_frame", "frame_status"], "").lower() in {"yes", "true", "1", "in_frame", "inframe"} else "no"
        if raw_type:
            effect_class = type_alias.get(raw_type.lower().replace("-", "_")) or raw_type
        elif _required(row, ["gene1", "left_gene"], "gene1") == _required(row, ["gene2", "right_gene"], "gene2"):
            effect_class = "SV_Junction" if in_frame == "yes" else "SV_Frameshift"
        else:
            effect_class = "SV_Fusion"
        out[key] = ExpressedProduct(
            adjacency_key=key,
            gene1=_required(row, ["gene1", "left_gene"], "gene1"),
            gene2=_required(row, ["gene2", "right_gene"], "gene2"),
            transcript1=_required(row, ["transcript1", "left_transcript"], "transcript1"),
            transcript2=_required(row, ["transcript2", "right_transcript"], "transcript2"),
            protein_sequence=sequence,
            wildtype_protein_sequence=first(row, ["wildtype_protein_sequence", "wt_protein_sequence"], ""),
            junction_aa_position=junction,
            in_frame=in_frame,
            orf_status=orf_status,
            source_tool=first(row, ["source_tool", "caller"], ""),
            source_record_id=first(row, ["source_record_id", "record_id"], ""),
            effect_class=effect_class,
        )
    return out


def product_meta(product: ExpressedProduct) -> dict[str, Any]:
    return {
        "gene1": product.gene1,
        "gene2": product.gene2,
        "transcript1": product.transcript1,
        "transcript2": product.transcript2,
        "effect_class": product.effect_class,
        "fusion_in_frame": product.in_frame if product.effect_class == "SV_Fusion" else "unknown",
        "frameshift": "no" if product.in_frame == "yes" else "yes",
        "junction_aa_position": str(product.junction_aa_position),
        "reconstruction_status": "confirmed_expressed_product",
        "reconstruction_reason": (
            f"exact RNA product from {product.source_tool or 'external tool'}; "
            f"record={product.source_record_id or 'not supplied'}"
        ),
    }

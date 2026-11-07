
"""Peptide-level safety evidence resolver.

The gate is deliberately conservative and explainable. It writes a sidecar
`peptide_safety.tsv` that score can consume, instead of hiding safety logic
inside a single scalar ranking score.
"""

from __future__ import annotations

import csv
import gzip
from pathlib import Path
import re
from typing import Any, Mapping

from .utils import read_tsv, write_tsv, to_float, first
from .safety import apply_event_safety, load_normal_expression
from .safety import combine_safety_reasons, worst_safety_status
from .adapters.peptide_input import normalize_hla_allele

PEPTIDE_SAFETY_FIELDS = [
    'peptide_id','event_id','sample_id','event_type','mutation_source','peptide_consequence','gene',
    'peptide','hla_allele','mhc_class','matched_normal_status','normal_alt_reads','normal_vaf','tumor_only_flag',
    'normal_support_evidence_basis','normal_support_integrity','safety_decision_basis',
    'reference_proteome_exact_match','reference_proteome_match_type','reference_match_gene','reference_match_protein',
    'reference_match_position','reference_match_length','normal_hla_ligand_exact_match','normal_ligand_tissue',
    'normal_ligand_hla','normal_ligand_source_protein','normal_ligand_class','normal_tissue_max_tpm',
    'normal_tissue_max_tissue','critical_tissue_max_tpm','normal_hspc_tpm','normal_hspc_unit','critical_tissue_hit','critical_tissue_name','normal_junction_seen',
    'normal_junction_match_level',
    'normal_junction_source','normal_junction_max_reads','normal_junction_tissue','wildtype_peptide','mt_binding_rank',
    'wt_binding_rank','mt_wt_fold_change','mutation_position_in_peptide','mutation_anchor_only','anchor_risk_status',
    'closest_self_peptide','closest_self_gene','closest_self_similarity','closest_self_hla_binding_rank',
    'closest_self_normal_expression_tpm','safety_tier','safety_status','safety_reason','safety_multiplier','review_required',
    'normal_expression_status','normal_hspc_status','reference_proteome_status','normal_ligandome_status',
    'anchor_assessment_status','normal_junction_assessment_status','safety_evidence_completeness',
    'normal_proteome_exact_match_status','normal_transcript_junction_match_status',
    'normal_immunopeptidome_match_status','similar_peptide_cross_reactivity_status',
    'source_gene_expression_context_status','critical_organ_expression_context_status',
    'hematopoietic_expression_context_status','tcr_contact_anchor_context_status',
    'final_safety_conclusion','safety_missing_layers','safety_priority_cap',
]

EVENT_SAFETY_FIELDS = [
    'event_id','sample_id','gene','event_type','mutation_source','normal_expression_status','normal_junction_status',
    'normal_junction_match_level',
    'matched_normal_status','event_safety_status','event_safety_reason','gene_tpm_policy',
    'normal_hspc_status','reference_proteome_status','normal_ligandome_status','anchor_assessment_status',
    'normal_tissue_max_tpm','normal_tissue_max_tissue','critical_tissue_max_tpm','critical_tissue_name',
    'normal_hspc_tpm','normal_hspc_unit','safety_evidence_completeness','safety_missing_layers'
]


def _read_fasta_sequences(path: str | Path | None) -> dict[str, str]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f'Missing reference proteome FASTA: {p}')
    seqs: dict[str, str] = {}
    name = None; chunks: list[str] = []
    with p.open('r', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(chunks).upper().replace('*','')
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            seqs[name] = ''.join(chunks).upper().replace('*','')
    return seqs


def build_reference_index(
    path: str | Path | None,
    lengths: tuple[int, ...] = (8,9,10,11),
    target_peptides: set[str] | None = None,
) -> dict[str, list[tuple[str,int]]]:
    seqs = _read_fasta_sequences(path)
    idx: dict[str, list[tuple[str,int]]] = {}
    for prot, seq in seqs.items():
        for L in lengths:
            if len(seq) < L:
                continue
            for i in range(0, len(seq)-L+1):
                pep = seq[i:i+L]
                if 'X' in pep or (target_peptides is not None and pep not in target_peptides):
                    continue
                idx.setdefault(pep, []).append((prot, i+1))
    return idx


def _normal_ligand_rows(
    path: str | Path | None,
    target_peptides: set[str] | None = None,
) -> dict[str, list[dict[str, str]]]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    out: dict[str, list[dict[str, str]]] = {}
    for r in read_tsv(p):
        pep = first(r, ['peptide','Peptide','sequence','Sequence'], '').upper().strip()
        if not pep or (target_peptides is not None and pep not in target_peptides):
            continue
        out.setdefault(pep, []).append({
            'peptide': pep,
            'hla_allele': first(r, ['hla_allele','HLA','allele'], ''),
            'tissue': first(r, ['tissue','source_tissue','normal_tissue','source'], ''),
            'source_protein': first(r, ['source_protein','protein','gene','source_gene'], ''),
            'hla_class': first(r, ['hla_class','mhc_class','class'], ''),
        })
    return out


def _normalize_chromosome(value: Any) -> str:
    chrom = str(value or '').strip()
    if not chrom:
        return ''
    if chrom.lower().startswith('chr'):
        chrom = chrom[3:]
    if chrom.upper() in {'M', 'MT'}:
        chrom = 'M'
    return f'chr{chrom}'


def _junction_coordinate_keys(row: Mapping[str, Any]) -> list[str]:
    chrom = first(row, ['chromosome', 'chrom', 'chr', 'seqname', 'seqnames'], '')
    start = first(row, ['junction_start', 'start', 'donor', 'intron_start'], '')
    end = first(row, ['junction_end', 'end', 'acceptor', 'intron_end'], '')
    strand = first(row, ['strand', 'junction_strand'], '')
    if not (chrom and str(start).strip() and str(end).strip()):
        embedded = first(row, ['junction_id', 'event_id'], '')
        match = re.search(r'(chr)?([0-9]+|X|Y|M|MT):(\d+)-(\d+)(?::([+?.-]))?', str(embedded), re.I)
        if match:
            chrom, start, end, embedded_strand = match.group(2), match.group(3), match.group(4), match.group(5)
            strand = strand or embedded_strand or ''
    chrom = _normalize_chromosome(chrom)
    try:
        start_i, end_i = int(str(start)), int(str(end))
    except (TypeError, ValueError):
        return []
    if not chrom or start_i < 0 or end_i < 0:
        return []
    if start_i > end_i:
        start_i, end_i = end_i, start_i
    base = f'{chrom}:{start_i}-{end_i}'
    strand = str(strand or '').strip()
    return [f'{base}:{strand}', base] if strand in {'+', '-'} else [base]


def _normal_junctions(
    path: str | Path | None,
    target_keys: set[str] | None = None,
) -> dict[str, dict[str, str]]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    out={}
    handle = gzip.open(p, 'rt', encoding='utf-8', newline='') if p.suffix == '.gz' else p.open('r', encoding='utf-8', newline='')
    with handle:
        header_line = handle.readline()
        if not header_line:
            return out
        fieldnames = header_line.rstrip('\r\n').split('\t')
        field_index = {name: idx for idx, name in enumerate(fieldnames)}

        def value(parts: list[str], *names: str) -> str:
            for name in names:
                idx = field_index.get(name)
                if idx is not None and idx < len(parts):
                    found = parts[idx].strip()
                    if found:
                        return found
            return ''

        has_explicit_coordinates = bool(
            set(field_index).intersection({'chromosome', 'chrom', 'chr', 'seqname', 'seqnames'})
            and set(field_index).intersection({'junction_start', 'start', 'donor', 'intron_start'})
            and set(field_index).intersection({'junction_end', 'end', 'acceptor', 'intron_end'})
        )
        for line in handle:
            parts = line.rstrip('\r\n').split('\t')
            normalized=[]
            for name in ('event_id', 'junction_id', 'gene_pair', 'fusion'):
                raw = value(parts, name)
                if not raw:
                    continue
                canonical = raw.replace('--', '::')
                normalized.extend([raw, canonical])
                if '::' in canonical:
                    left, right = canonical.split('::', 1)
                    normalized.append(f'{right}::{left}')
            gene1, gene2 = value(parts, 'gene1'), value(parts, 'gene2')
            if gene1 and gene2:
                normalized.extend([f'{gene1}::{gene2}', f'{gene2}::{gene1}'])

            if has_explicit_coordinates:
                chrom = _normalize_chromosome(value(parts, 'chromosome', 'chrom', 'chr', 'seqname', 'seqnames'))
                start = value(parts, 'junction_start', 'start', 'donor', 'intron_start')
                end = value(parts, 'junction_end', 'end', 'acceptor', 'intron_end')
                strand = value(parts, 'strand', 'junction_strand')
                try:
                    start_i, end_i = int(start), int(end)
                except (TypeError, ValueError):
                    start_i = end_i = -1
                if chrom and start_i >= 0 and end_i >= 0:
                    if start_i > end_i:
                        start_i, end_i = end_i, start_i
                    base = f'{chrom}:{start_i}-{end_i}'
                    normalized.append(base)
                    if strand in {'+', '-'}:
                        normalized.append(f'{base}:{strand}')
            elif target_keys is None:
                row = dict(zip(fieldnames, parts))
                normalized.extend(_junction_coordinate_keys(row))

            if target_keys is not None:
                normalized = [key for key in normalized if key in target_keys]
                if not normalized:
                    continue
            row = dict(zip(fieldnames, parts))
            for key in dict.fromkeys(normalized):
                out[key] = row
    return out


def _normal_junction_scopes(path: str | Path | None) -> set[str]:
    if not path or not Path(path).is_file():
        return set()
    p = Path(path)
    handle = gzip.open(p, 'rt', encoding='utf-8', newline='') if p.suffix == '.gz' else p.open('r', encoding='utf-8', newline='')
    with handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if 'junction_class' not in (reader.fieldnames or []):
            return {'all_junctions'}
        scopes = {
            str(row.get('junction_class') or '').strip()
            for row in reader
            if str(row.get('junction_class') or '').strip()
        }
    # Backward compatibility for user-provided files predating scoped refs.
    return scopes or {'all_junctions'}


def _normal_junction_kind(
    peptide: Mapping[str, Any],
    event: Mapping[str, Any],
) -> str:
    consequence = str(peptide.get('peptide_consequence') or '').lower()
    event_type = str(peptide.get('event_type') or event.get('event_type') or '').lower()
    mutation_source = str(peptide.get('mutation_source') or event.get('mutation_source') or '').lower()
    if consequence == 'exon_deletion_junction':
        return 'not_applicable'
    # VEP can label an SNV/InDel as splice_region_variant while the peptide
    # catalog labels its protein product as splice_junction. Without an
    # observed RNA exon-exon junction this is still a sequence variant, not a
    # junction that can be looked up in GTEx.
    if event_type in {'snv', 'indel'} and mutation_source in {'snv', 'indel', 'wes', 'wgs'}:
        return 'not_applicable'
    if consequence == 'splice_junction' or event_type in {'splice', 'splice_junction'}:
        return 'splice'
    if consequence in {'fusion', 'sv_junction'} or event_type in {'fusion', 'sv_junction'}:
        return 'fusion'
    # crosses_junction also marks some variant-generated InDel peptides. Those
    # are not normal RNA splice/fusion junctions and must not require this layer.
    if peptide.get('crosses_junction') == 'yes' and mutation_source in {'sv', 'fusion'}:
        return 'fusion'
    if peptide.get('crosses_junction') == 'yes' and mutation_source == 'splice':
        return 'splice'
    return 'not_applicable'


_JUNCTION_COORD_KEY_RE = re.compile(r'^(chr[0-9A-Za-z]+:\d+-\d+)(?::([+?.-]))?$')


def _junction_row_strand(row: Mapping[str, Any]) -> str:
    """Strand of a normal-junction row: explicit column, then an embedded
    :+ / :- suffix in the junction identifier (review B3b-08)."""
    strand = str(first(row, ['strand', 'junction_strand'], '') or '').strip()
    if strand in {'+', '-'}:
        return strand
    embedded = str(first(row, ['junction_id', 'event_id'], '') or '')
    m = re.search(r':([+-])(?::|$)', embedded)
    return m.group(1) if m else ''


def _candidate_junction_keys(
    peptide: Mapping[str, Any],
    event: Mapping[str, Any],
) -> list[str]:
    keys: list[str] = []
    for row in (peptide, event):
        for field in ('junction_id', 'event_id', 'event_name', 'gene_pair', 'fusion'):
            value = str(row.get(field) or '').strip()
            if value and value not in keys:
                keys.append(value)
        for value in _junction_coordinate_keys(row):
            if value not in keys:
                keys.append(value)
    return keys


def _split_genes(gene: str) -> list[str]:
    raw = str(gene or '')
    for sep in ('::',';','|',','):
        raw = raw.replace(sep, ' ')
    genes: list[str] = []
    for token in (x.strip() for x in raw.split() if x.strip()):
        candidates: list[str] = []
        if '(' in token and token.endswith(')'):
            base, _, rest = token.partition('(')
            inner = rest[:-1].strip()
            if base.strip():
                candidates.append(base.strip())
            if inner and not inner.isdigit():
                candidates.append(inner)
        else:
            candidates.append(token)
        for item in candidates:
            if item and item not in genes:
                genes.append(item)
    return genes


def _normal_expression_aliases(
    peptide: Mapping[str, Any],
    event: Mapping[str, Any],
    normal_expr: Mapping[str, Mapping[str, Any]],
) -> str:
    """Return only exact gene aliases represented in the normal reference.

    Splice normalizers can retain a coordinate in ``gene`` while preserving
    the Ensembl gene ID in source provenance.  Safety assessment must recover
    that exact ID before expression/HSPC lookup; coordinate or nearby-gene
    inference is deliberately not allowed here.
    """
    candidates: list[str] = []
    fields = (
        'gene', 'gene_symbol', 'gene_name', 'ensembl_gene_id', 'gene_id',
        'gene_pair', 'source_event_id', 'source_record_id', 'source_records',
        'event_name',
    )
    for row in (peptide, event):
        for field_name in fields:
            value = str(row.get(field_name) or '').strip()
            if not value:
                continue
            candidates.extend(_split_genes(value))
            candidates.extend(re.findall(r'ENSG\d+(?:\.\d+)?', value))
    normalized: list[str] = []
    for candidate in candidates:
        for alias in (candidate, re.sub(r'^(ENSG\d+)\.\d+$', r'\1', candidate)):
            if alias in normal_expr and alias not in normalized:
                normalized.append(alias)
    if normalized:
        return ';'.join(normalized)
    return str(peptide.get('gene') or event.get('gene') or '')


def _normal_expr_for_genes(gene_value: str, normal_expr: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    nt=0.0; nh=0.0; crit=False; found=False; critical_tpm=0.0
    nt_tissue=''; critical_tissue=''; hspc_unit=''; expr_assessed=True; hspc_assessed=True
    genes = _split_genes(gene_value)
    for g in genes:
        row = normal_expr.get(g) or {}
        if row:
            found=True
            row_nt=to_float(row.get('normal_tissue_max_tpm'), 0.0)
            if row_nt >= nt: nt=row_nt; nt_tissue=row.get('normal_tissue_max_tissue','')
            row_hspc=to_float(row.get('normal_hspc_tpm'), 0.0)
            if row_hspc >= nh: nh=row_hspc; hspc_unit=row.get('normal_hspc_unit','TPM')
            row_critical=to_float(row.get('critical_tissue_max_tpm'), row_nt)
            if row_critical >= critical_tpm: critical_tpm=row_critical; critical_tissue=row.get('critical_tissue_name','')
            crit=crit or to_float(row.get('critical_tissue_hit'), 0.0) > 0
            expr_assessed=expr_assessed and row.get('normal_expression_status','ASSESSED') == 'ASSESSED'
            hspc_assessed=hspc_assessed and row.get('normal_hspc_status','ASSESSED') == 'ASSESSED'
        else:
            expr_assessed=False; hspc_assessed=False
    return {
        'found': found, 'normal_tissue_max_tpm': nt, 'normal_tissue_max_tissue': nt_tissue,
        'critical_tissue_max_tpm': critical_tpm, 'critical_tissue_name': critical_tissue,
        'normal_hspc_tpm': nh, 'normal_hspc_unit': hspc_unit or 'TPM', 'critical_tissue_hit': crit,
        'normal_expression_status': 'ASSESSED' if found and expr_assessed else 'UNASSESSED',
        'normal_hspc_status': 'ASSESSED' if found and hspc_assessed else 'UNASSESSED',
    }


def _find_normal_ligand(
    peptide: str,
    hla: str,
    ligs: Mapping[str, list[dict[str, str]]],
) -> dict[str, str] | None:
    h_norm = normalize_hla_allele(hla)
    exact = ligs.get(peptide, [])
    if not exact:
        return None
    for r in exact:
        if normalize_hla_allele(r.get('hla_allele','')) == h_norm and h_norm:
            return r
    return exact[0]


def _event_normal_support(
    event: Mapping[str, Any],
    cfg: Mapping[str, Any] | None = None,
) -> tuple[bool, str, str, str, bool, str, str]:
    """Classify matched-normal ALT evidence (review C6-04, round 3).

    The old reader coerced a MISSING normal VAF to 0.0, so an ALT-positive
    record with no VAF was automatically classified as normal-positive, and
    a thin normal depth read as a measured normal-negative.  Missing,
    measured, thin and contradictory inputs are now distinct states:

    - not_assessed               no normal evidence column is present at all
    - normal_alt_support         reliable matched-normal support: ALT reads at
                                 or above the read threshold with a reported
                                 or explicitly recalculated VAF at or above
                                 the VAF threshold (same allele's counts)
    - NORMAL_SUPPORT_REVIEW      ALT-positive record whose VAF is not
                                 available and cannot be recalculated, or
                                 whose measured VAF contradicts the read
                                 count (>= threshold reads, sub-threshold
                                 VAF) - neither a normal negative nor a
                                 confirmed normal positive until validated
    - NORMAL_VAF_CONTRADICTION   reported VAF and counts cannot both be true
                                 (explicit integrity warning, never resolved
                                 silently in either direction)
    - NORMAL_COVERAGE_INADEQUATE normal depth present but below the minimum:
                                 no claim either way (missing evidence is not
                                 a measured negative)
    - no_clear_normal_support    measured and below the positive threshold

    Returns (support, alt, vaf, status, review, integrity, evidence_basis).
    """
    cfg = cfg or {}
    alt_raw = str(event.get('normal_alt_count') or event.get('normal_alt_support') or '').strip()
    depth_raw = str(event.get('normal_depth') or '').strip()
    vaf_raw = str(event.get('normal_vaf') or event.get('normal_alt_vaf') or '').strip()
    if not (alt_raw or depth_raw or vaf_raw):
        return False, '', '', 'not_assessed', False, '', 'NO_NORMAL_EVIDENCE_COLUMNS'
    n_alt = to_float(alt_raw, 0.0) if alt_raw else None
    depth = to_float(depth_raw, 0.0) if depth_raw else None
    vaf = to_float(vaf_raw, 0.0) if vaf_raw else None
    pos_reads = float(cfg.get('normal_alt_positive_reads', 3))
    pos_vaf = float(cfg.get('normal_vaf_positive', 0.02))
    min_depth = float(cfg.get('normal_support_min_depth', 10))
    alt_str = str(int(n_alt)) if n_alt is not None else ''
    basis_parts = ['REPORTED' if vaf_raw else 'VAF_UNAVAILABLE']
    integrity = ''
    recalculated = False
    if vaf is None and n_alt is not None and depth is not None and depth > 0:
        # Explicit recalculation from the same allele's counts and depth.
        vaf = n_alt / depth
        basis_parts = [f'VAF_RECALCULATED_ALT{int(n_alt)}_DEPTH{int(depth)}']
        recalculated = True
    # Contradiction detection between reported counts and VAF (the VAF was
    # reported, not recalculated here, so any disagreement is the data's).
    if not recalculated and vaf is not None and n_alt is not None and depth is not None and depth > 0:
        recomputed = n_alt / depth
        gap = abs(recomputed - vaf)
        if gap > 0.1 and gap > 0.5 * max(recomputed, vaf):
            integrity = (
                f'normal_vaf_contradiction: reported VAF {vaf:.4f} is incompatible with '
                f'{int(n_alt)} ALT reads at depth {int(depth)} (recalculated {recomputed:.4f})'
            )
    if not integrity and vaf is not None and vaf >= pos_vaf and (n_alt or 0) <= 0:
        integrity = (
            f'normal_vaf_contradiction: reported VAF {vaf:.4f} with zero normal ALT reads'
        )
    status = ''
    review = False
    support = False
    if integrity:
        status = 'NORMAL_VAF_CONTRADICTION'
        review = True
    elif (n_alt or 0) >= pos_reads:
        if vaf is not None and vaf >= pos_vaf:
            support = True
            status = 'normal_alt_support'
        else:
            # Positive ALT reads with an unavailable or sub-threshold VAF:
            # needs validation; never auto-pass and never auto-reject.
            status = 'NORMAL_SUPPORT_REVIEW'
            review = True
    elif (n_alt or 0) > 0 and vaf is None:
        status = 'NORMAL_SUPPORT_REVIEW'
        review = True
    elif depth is not None and depth < min_depth:
        status = 'NORMAL_COVERAGE_INADEQUATE'
    else:
        status = 'no_clear_normal_support'
    vaf_str = f'{vaf:.4f}' if vaf is not None else ''
    return support, alt_str, vaf_str, status, review, integrity, ';'.join(basis_parts)


def _anchor_only(peptide: str, wt_peptide: str, hla: str) -> tuple[str, str, str]:
    if not peptide or not wt_peptide or len(peptide) != len(wt_peptide):
        return 'unknown', '', 'unassessed'
    diffs = [i+1 for i,(a,b) in enumerate(zip(peptide, wt_peptide)) if a != b]
    if not diffs:
        return 'no', '', 'no_difference'
    L = len(peptide)
    # Approximate MHC-I anchor positions; class/allele-specific anchor maps can replace this later.
    anchors = {2, L}
    only = all(d in anchors for d in diffs)
    return 'yes' if only else 'no', ','.join(map(str,diffs)), 'anchor_only' if only else 'tcr_facing_change_present'


def safety_multiplier_for(status: str, profile: Mapping[str, Any] | None = None) -> float:
    cfg = (profile or {}).get('safety', {})
    s = str(status or '').upper()
    if s in {'SAFETY_PASS','PASS'}:
        return 1.0
    if s in {'SAFETY_REVIEW','REVIEW','CAUTION'}:
        return float(cfg.get('caution_multiplier', 0.55))
    if s in {'SAFETY_PARTIAL','PARTIAL'}:
        return float(cfg.get('partial_multiplier', 0.75))
    if s in {'SAFETY_HIGH_RISK','HIGH_RISK'}:
        return 0.20
    if s in {'SAFETY_REJECT','REJECT','FAIL'}:
        return 0.0
    return 0.60


def build_peptide_safety_gate(
    *,
    raw_events: str | Path,
    raw_peptides: str | Path,
    out_peptide_safety: str | Path,
    out_event_safety: str | Path | None = None,
    profile: Mapping[str, Any] | None = None,
    normal_expression: str | Path | None = None,
    normal_hla_ligands: str | Path | None = None,
    reference_proteome: str | Path | None = None,
    normal_junctions: str | Path | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    events = {r.get('event_id',''): r for r in read_tsv(raw_events)}
    peptides = read_tsv(raw_peptides)
    norm_expr = load_normal_expression(normal_expression)
    target_peptides = {(row.get('peptide') or '').upper().strip() for row in peptides if row.get('peptide')}
    ligs = _normal_ligand_rows(normal_hla_ligands, target_peptides=target_peptides)
    ref_idx = build_reference_index(reference_proteome, target_peptides=target_peptides)
    target_junction_keys: set[str] = set()
    for peptide_row in peptides:
        event_row = events.get(peptide_row.get('event_id', ''), {})
        if _normal_junction_kind(peptide_row, event_row) == 'not_applicable':
            continue
        target_junction_keys.update(_candidate_junction_keys(peptide_row, event_row))
        gene_pair = str(event_row.get('gene') or '').strip()
        if gene_pair:
            target_junction_keys.add(gene_pair)
    junctions = (
        _normal_junctions(normal_junctions, target_keys=target_junction_keys)
        if target_junction_keys else {}
    )
    junction_scopes = _normal_junction_scopes(normal_junctions) if target_junction_keys else set()
    ref_assessed = bool(reference_proteome and Path(reference_proteome).is_file())
    ligand_assessed = bool(normal_hla_ligands and Path(normal_hla_ligands).is_file())
    safety_cfg = (profile or {}).get('safety', {})
    ref_reject = str(safety_cfg.get('reference_proteome_exact_match_action','FAIL')).upper() in {'FAIL','REJECT'}
    normal_lig_action = str(safety_cfg.get('normal_hla_ligand_overlap_action','CAUTION')).upper()
    # Review C6-03: a normal junction keyed only by the candidate's gene pair
    # is context, not proof of the candidate's exact junction. "reject" (the
    # default) keeps the historical production behavior but labels the match
    # level honestly; "review" caps gene-pair-only matches at CAUTION/review.
    # Review C6-03 (round 3, FURTHER_SUGGESTIONS WS8): the operator decision
    # makes review-level context the DEFAULT - a recurrent normal fusion of
    # the same gene pair at a different breakpoint does not prove THIS
    # candidate's peptide exists in normal tissue.  "reject" remains available
    # as an explicit conservative override.
    gene_pair_action = str(safety_cfg.get('normal_junction_gene_pair_only_action','review')).lower()
    rows=[]; event_rows=[]
    for p in peptides:
        e = events.get(p.get('event_id',''), {})
        peptide = (p.get('peptide') or '').upper().strip()
        hla = p.get('hla_allele','')
        reasons=[]; status='SAFETY_PASS'; review='no'
        normal_support, normal_alt, normal_vaf, normal_status, normal_review, normal_integrity, normal_basis = _event_normal_support(e, safety_cfg)
        if normal_support:
            status='SAFETY_REJECT'; reasons.append('matched_normal_support')
        elif normal_review:
            # Review C6-04: ALT-positive / VAF-unavailable (or contradictory)
            # normal evidence goes to review, not to a silent pass or a
            # confirmed normal-positive call.
            if status == 'SAFETY_PASS': status='SAFETY_REVIEW'
            reasons.append('matched_normal_support_review')
            if normal_integrity:
                reasons.append('normal_vaf_integrity_conflict')
        ref_hits = ref_idx.get(peptide, []) if peptide else []
        if ref_hits and ref_reject:
            status='SAFETY_REJECT'; reasons.append('reference_proteome_exact_match')
        lig = _find_normal_ligand(peptide, hla, ligs) if peptide else None
        if lig:
            if normal_lig_action in {'FAIL','REJECT'} or normalize_hla_allele(lig.get('hla_allele','')) == normalize_hla_allele(hla):
                if status != 'SAFETY_REJECT':
                    status='SAFETY_HIGH_RISK'
                reasons.append('normal_hla_ligand_exact_match')
            else:
                if status == 'SAFETY_PASS': status='SAFETY_REVIEW'
                reasons.append('normal_hla_ligand_overlap')
        gene_pair = e.get('gene','')
        candidate_junction_keys = _candidate_junction_keys(p, e)
        jrow = None
        junction_match_level = ''
        # Review B3b-08: a verified normal-junction match requires the full
        # strand-aware key.  The candidate's strand is resolved when its own
        # annotation (junction identifier suffix, transcript/exon structure
        # upstream, or an explicit strand column) states it uniquely; an
        # unresolved strand never earns an exact match, in either direction.
        cand_strands = {
            m.group(2) for key in candidate_junction_keys
            if (m := _JUNCTION_COORD_KEY_RE.fullmatch(key)) and m.group(2) in {'+', '-'}
        }
        cand_strand = next(iter(cand_strands)) if len(cand_strands) == 1 else ''
        for key in candidate_junction_keys:
            hit = junctions.get(key)
            if not hit:
                continue
            coord = _JUNCTION_COORD_KEY_RE.fullmatch(key)
            if coord is None:
                # Non-coordinate key (junction/fusion identifier or gene-pair
                # string): the identifier itself matched exactly.
                jrow = hit
                junction_match_level = 'EXACT_NORMAL_JUNCTION_MATCH'
                break
            normal_strand = _junction_row_strand(hit) or (coord.group(2) if coord.group(2) in {'+', '-'} else '')
            if cand_strand and normal_strand in {'+', '-'} and normal_strand != cand_strand:
                # Both strands resolved and OPPOSITE: coordinates equal is not
                # this junction. Not even a possible match.
                continue
            if cand_strand and normal_strand == cand_strand:
                jrow = hit
                junction_match_level = 'EXACT_NORMAL_JUNCTION_MATCH'
                break
            # Coordinates match with a strand unknown on at least one side.
            jrow = hit
            junction_match_level = 'POSSIBLE_NORMAL_JUNCTION_STRAND_AMBIGUOUS'
            break
        if not jrow and gene_pair and junctions.get(gene_pair):
            # Review C6-03: keyed only by the gene pair, this normal junction
            # cannot prove the candidate's specific breakpoint occurs in
            # normal tissue - it only says the pair is junction-capable.
            jrow = junctions[gene_pair]
            junction_match_level = 'GENE_PAIR_ONLY_CONTEXT'
        junction_kind = _normal_junction_kind(p, e)
        is_junction_type = junction_kind != 'not_applicable'
        if junction_kind == 'fusion':
            junction_assessed = bool({'normal_recurrent_fusion', 'all_junctions'} & junction_scopes)
        elif junction_kind == 'splice':
            junction_assessed = bool({'normal_splice_junction', 'all_junctions'} & junction_scopes)
        else:
            junction_assessed = False
        if jrow and is_junction_type:
            max_reads = to_float(first(jrow, ['normal_reads','normal_sample_count','junction_reads','reads'], '0'), 0.0)
            gene_pair_only = junction_match_level == 'GENE_PAIR_ONLY_CONTEXT'
            strand_ambiguous = junction_match_level == 'POSSIBLE_NORMAL_JUNCTION_STRAND_AMBIGUOUS'
            if strand_ambiguous:
                # Review B3b-08: an unresolved strand makes the match neither
                # an exact positive nor a definitive negative - review level,
                # no matter how many reads the catalog row carries.
                if status == 'SAFETY_PASS': status='SAFETY_REVIEW'
                reasons.append('possible_normal_junction_strand_ambiguous')
            elif gene_pair_only and gene_pair_action == 'review':
                if status == 'SAFETY_PASS': status='SAFETY_REVIEW'
                reasons.append('gene_pair_only_normal_junction_context')
            elif max_reads >= float(safety_cfg.get('normal_junction_reject_reads', 3)):
                status='SAFETY_REJECT'
                reasons.append('normal_junction_seen_gene_pair_only' if gene_pair_only else 'normal_junction_seen')
            else:
                if status == 'SAFETY_PASS': status='SAFETY_REVIEW'
                reasons.append('low_level_normal_junction_gene_pair_only' if gene_pair_only else 'low_level_normal_junction_seen')
        expression_aliases = _normal_expression_aliases(p, e, norm_expr)
        expr = _normal_expr_for_genes(expression_aliases, norm_expr)
        wt_rank = to_float(p.get('wildtype_binding_rank') or p.get('netmhcpan_wt_rank_el') or p.get('netmhcpan_wt_rank_ba'), 99.0)
        mt_rank = to_float(p.get('netmhcpan_el_rank') or p.get('netmhcpan_mt_rank_el') or p.get('binding_rank'), 99.0)
        anchor_only, mut_pos, anchor_status = _anchor_only(peptide, p.get('wildtype_peptide',''), hla)
        anchor_heuristic_hit = anchor_only == 'yes' and wt_rank <= float(safety_cfg.get('wildtype_strong_binding_rank', 0.5))
        if anchor_heuristic_hit:
            # Review C6-05: the position-2/C-terminal rule is an approximate
            # cross-allele heuristic, not an allele-specific structural or
            # functional validation.  It stays a caution signal by default;
            # only an explicit, separately justified profile knob may let it
            # alone drive a hard rejection, and "record_only" keeps it purely
            # documentary.
            anchor_action = str(safety_cfg.get('anchor_only_action', 'review')).lower()
            reasons.append('anchor_only_mutation_with_wt_binding')
            if anchor_action == 'reject':
                status = 'SAFETY_REJECT'
                reasons.append('anchor_only_reject_explicitly_enabled')
            elif anchor_action != 'record_only' and status == 'SAFETY_PASS':
                status = 'SAFETY_REVIEW'
        sim_raw = str(p.get('self_similarity_score') or '').strip()
        sim = to_float(sim_raw, 0.0)
        similarity_cutoff = float(safety_cfg.get('self_similarity_high_risk', safety_cfg.get('self_similarity_caution', 0.85)))
        if sim_raw and sim >= similarity_cutoff and expr.get('critical_tissue_hit'):
            if status not in {'SAFETY_REJECT'}: status='SAFETY_HIGH_RISK'
            reasons.append('high_self_similarity_critical_tissue')
        expression_status = expr.get('normal_expression_status', 'UNASSESSED')
        hspc_status = expr.get('normal_hspc_status', 'UNASSESSED')
        expression_reference_available = bool(normal_expression and Path(normal_expression).is_file())
        if not expr.get('found'):
            if expression_reference_available:
                reasons.extend([
                    'normal_expression_gene_not_in_reference',
                    'normal_hspc_gene_not_in_reference',
                ])
            else:
                reasons.extend([
                    'normal_expression_reference_not_provided',
                    'normal_hspc_reference_not_provided',
                ])
        reference_status = 'ASSESSED' if ref_assessed else 'UNASSESSED'
        ligandome_status = 'ASSESSED' if ligand_assessed else 'UNASSESSED'
        junction_type = is_junction_type
        junction_status = ('ASSESSED' if junction_assessed else 'UNASSESSED') if junction_type else 'NOT_APPLICABLE'
        if anchor_status != 'unassessed':
            anchor_assessment = 'ASSESSED'
        elif junction_type or str(p.get('peptide_consequence','')).lower() in {'frameshift','insertion'}:
            anchor_assessment = 'NOT_APPLICABLE'
        else:
            anchor_assessment = 'UNASSESSED'
        direct_safety_layers = {
            'reference_proteome': reference_status,
            'normal_ligandome': ligandome_status,
            'normal_junction': junction_status,
            'similar_self_peptide': 'ASSESSED' if sim_raw else 'UNASSESSED',
        }
        required_layers = dict(direct_safety_layers)
        if not junction_type:
            required_layers.update({
                'normal_expression': expression_status,
                'normal_hspc': hspc_status,
                'anchor_only': anchor_assessment,
            })
        assessed = [name for name, layer_status in required_layers.items() if layer_status != 'NOT_APPLICABLE']
        missing = [name for name in assessed if required_layers[name] != 'ASSESSED']
        completeness = (len(assessed) - len(missing)) / max(len(assessed), 1)
        if missing:
            status = worst_safety_status(status, 'SAFETY_PARTIAL')
            reasons.append('safety_evidence_incomplete')
        if status in {'SAFETY_REVIEW','SAFETY_HIGH_RISK','SAFETY_PARTIAL'}:
            review='yes'
        if not reasons:
            reasons.append('no_major_signal')
        direct_reason_tokens = {
            'reference_proteome_exact_match', 'normal_hla_ligand_exact_match',
            'normal_hla_ligand_overlap', 'normal_junction_seen',
            'low_level_normal_junction_seen', 'high_self_similarity_critical_tissue',
            'matched_normal_support',
            # Review C6-04: review/constraint calls that rest on directly
            # counted normal evidence.
            'matched_normal_support_review',
            'normal_vaf_integrity_conflict',
            # Review C6-03 labels for gene-pair-only matches.
            'normal_junction_seen_gene_pair_only',
            'low_level_normal_junction_gene_pair_only',
            'gene_pair_only_normal_junction_context',
            # Review B3b-08: coordinate-direct match with unresolved strand.
            'possible_normal_junction_strand_ambiguous',
        }
        # Review C6-05 (round 3): every decision names its epistemic basis -
        # direct molecular evidence, an approximation (heuristic), or neither.
        heuristic_reason_tokens = {
            'anchor_only_mutation_with_wt_binding', 'anchor_only_reject_explicitly_enabled',
        }
        if any(reason in direct_reason_tokens for reason in reasons):
            safety_decision_basis = 'DIRECT_ASSESSED_EVIDENCE'
        elif any(reason in heuristic_reason_tokens for reason in reasons):
            safety_decision_basis = 'HEURISTIC_ONLY'
        elif missing:
            safety_decision_basis = 'EVIDENCE_COMPLETENESS_ONLY'
        else:
            safety_decision_basis = 'NO_SAFETY_SIGNAL'
        if status == 'SAFETY_REJECT':
            final_safety_conclusion = 'REJECT_DIRECT_SAFETY_EVIDENCE'
        elif any(reason in direct_reason_tokens for reason in reasons):
            final_safety_conclusion = 'REVIEW_DIRECT_SAFETY_EVIDENCE'
        elif missing:
            final_safety_conclusion = 'PARTIAL_DIRECT_SAFETY_EVIDENCE'
        else:
            final_safety_conclusion = 'NO_DIRECT_SAFETY_SIGNAL_DETECTED'
        row={
            'peptide_id': p.get('peptide_id',''), 'event_id': p.get('event_id',''), 'sample_id': p.get('sample_id',''),
            'event_type': p.get('event_type',''), 'mutation_source': p.get('mutation_source',''),
            'peptide_consequence': p.get('peptide_consequence',''), 'gene': p.get('gene',''), 'peptide': peptide,
            'hla_allele': hla, 'mhc_class': p.get('mhc_class',''), 'matched_normal_status': normal_status,
            'normal_alt_reads': normal_alt, 'normal_vaf': normal_vaf, 'tumor_only_flag': 'yes' if not e else 'no',
            # Review C6-04: how the normal VAF was obtained (reported vs
            # recalculated from the same allele's counts) and any integrity
            # conflict; review C6-05: direct-vs-heuristic decision basis.
            'normal_support_evidence_basis': normal_basis,
            'normal_support_integrity': normal_integrity,
            'safety_decision_basis': safety_decision_basis,
            'reference_proteome_exact_match': 'yes' if ref_hits else 'no' if ref_assessed else 'not_assessed',
            'reference_proteome_match_type': 'exact' if ref_hits else '',
            'reference_match_gene': '', 'reference_match_protein': ref_hits[0][0] if ref_hits else '',
            'reference_match_position': str(ref_hits[0][1]) if ref_hits else '', 'reference_match_length': str(len(peptide)) if ref_hits else '',
            'normal_hla_ligand_exact_match': 'yes' if lig else 'no' if ligand_assessed else 'not_assessed',
            'normal_ligand_tissue': (lig or {}).get('tissue',''), 'normal_ligand_hla': (lig or {}).get('hla_allele',''),
            'normal_ligand_source_protein': (lig or {}).get('source_protein',''), 'normal_ligand_class': (lig or {}).get('hla_class',''),
            'normal_tissue_max_tpm': f"{expr.get('normal_tissue_max_tpm',0.0):.4f}" if expr.get('found') else '',
            'normal_tissue_max_tissue': expr.get('normal_tissue_max_tissue',''),
            'critical_tissue_max_tpm': f"{expr.get('critical_tissue_max_tpm',0.0):.4f}" if expr.get('found') else '',
            'normal_hspc_tpm': f"{expr.get('normal_hspc_tpm',0.0):.4f}" if expr.get('found') else '',
            'normal_hspc_unit': expr.get('normal_hspc_unit','') if expr.get('found') else '',
            'critical_tissue_hit': ('yes' if expr.get('critical_tissue_hit') else 'no') if expr.get('found') else 'not_assessed',
            'critical_tissue_name': expr.get('critical_tissue_name',''),
            'normal_junction_seen': (
                'yes' if jrow and is_junction_type else
                'no' if junction_assessed else
                'not_applicable' if junction_kind == 'not_applicable' else
                'not_assessed'
            ),
            'normal_junction_source': first(jrow or {}, ['source','dataset'], '') if is_junction_type else '',
            'normal_junction_max_reads': first(jrow or {}, ['normal_reads','normal_sample_count','junction_reads','reads'], ''),
            'normal_junction_tissue': first(jrow or {}, ['tissue','normal_tissue'], ''),
            'normal_junction_match_level': junction_match_level if (jrow and is_junction_type) else '',
            'wildtype_peptide': p.get('wildtype_peptide',''), 'mt_binding_rank': f'{mt_rank:.4f}', 'wt_binding_rank': f'{wt_rank:.4f}',
            'mt_wt_fold_change': '', 'mutation_position_in_peptide': mut_pos, 'mutation_anchor_only': anchor_only,
            'anchor_risk_status': anchor_status, 'closest_self_peptide': '', 'closest_self_gene': '',
            'closest_self_similarity': p.get('self_similarity_score',''), 'closest_self_hla_binding_rank': '',
            'closest_self_normal_expression_tpm': '', 'safety_tier': status, 'safety_status': 'FAIL' if status=='SAFETY_REJECT' else ('CAUTION' if status in {'SAFETY_REVIEW','SAFETY_HIGH_RISK'} else ('SAFETY_PARTIAL' if status=='SAFETY_PARTIAL' else 'PASS')),
            'safety_reason': ';'.join(dict.fromkeys(reasons)), 'safety_multiplier': f'{safety_multiplier_for(status, profile):.4f}',
            'review_required': review,
            'normal_expression_status': expression_status, 'normal_hspc_status': hspc_status,
            'reference_proteome_status': reference_status, 'normal_ligandome_status': ligandome_status,
            'anchor_assessment_status': anchor_assessment, 'normal_junction_assessment_status': junction_status,
            'safety_evidence_completeness': f'{completeness:.4f}', 'safety_missing_layers': ';'.join(missing),
            'safety_priority_cap': str(safety_cfg.get('partial_priority_cap','C_CAUTION')) if missing else '',
            'normal_proteome_exact_match_status': 'DETECTED' if ref_hits else 'NOT_DETECTED' if ref_assessed else 'UNASSESSED',
            # Review WS8: DETECTED now means EXACT strand-aware detection
            # only.  Gene-pair context and strand-ambiguous coordinate hits
            # keep their own vocabulary so no downstream consumer can read a
            # contextual signal as exact junction detection.
            'normal_transcript_junction_match_status': (
                'DETECTED' if jrow and is_junction_type and junction_match_level == 'EXACT_NORMAL_JUNCTION_MATCH' else
                'POSSIBLE_STRAND_AMBIGUOUS' if jrow and is_junction_type and junction_match_level == 'POSSIBLE_NORMAL_JUNCTION_STRAND_AMBIGUOUS' else
                'GENE_PAIR_CONTEXT_NOT_EXACT' if jrow and is_junction_type else
                'NOT_DETECTED' if junction_assessed else
                'NOT_APPLICABLE' if not is_junction_type else 'UNASSESSED'
            ),
            'normal_immunopeptidome_match_status': 'DETECTED' if lig else 'NOT_DETECTED' if ligand_assessed else 'UNASSESSED',
            'similar_peptide_cross_reactivity_status': (
                'HIGH_SIMILARITY_REVIEW' if sim_raw and sim >= similarity_cutoff else
                'LOW_SIMILARITY' if sim_raw else 'UNASSESSED'
            ),
            'source_gene_expression_context_status': expression_status,
            'critical_organ_expression_context_status': (
                'DETECTED_AUXILIARY' if expr.get('found') and expr.get('critical_tissue_hit') else
                'NOT_DETECTED_AUXILIARY' if expr.get('found') else 'UNASSESSED'
            ),
            'hematopoietic_expression_context_status': (
                'DETECTED_AUXILIARY' if expr.get('found') and to_float(expr.get('normal_hspc_tpm'), 0.0) > 0 else
                'NOT_DETECTED_AUXILIARY' if expr.get('found') else 'UNASSESSED'
            ),
            'tcr_contact_anchor_context_status': anchor_status.upper(),
            'final_safety_conclusion': final_safety_conclusion,
        }
        rows.append(row)
    # event-level rollup
    by_event: dict[str, list[dict[str, Any]]] = {}
    for r in rows:
        by_event.setdefault(r['event_id'], []).append(r)
    for eid, rs in by_event.items():
        e = events.get(eid, {})
        event_for_safety = dict(e)
        event_for_safety['gene'] = _normal_expression_aliases({}, e, norm_expr)
        event_assessment = apply_event_safety(event_for_safety, profile or {}, norm_expr)
        worst = event_assessment.get('safety_status', 'SAFETY_PARTIAL')
        reasons = [event_assessment.get('safety_reason', '')]
        matched_normal = any(r.get('matched_normal_status') == 'normal_alt_support' for r in rs)
        if matched_normal:
            worst = worst_safety_status(worst, 'FAIL')
            reasons.append('matched_normal_support')
        elif any(
            r.get('matched_normal_status') in {'NORMAL_SUPPORT_REVIEW', 'NORMAL_VAF_CONTRADICTION'}
            for r in rs
        ):
            # Review C6-04: the event inherits the review, not the rejection.
            worst = worst_safety_status(worst, 'CAUTION')
            reasons.append('matched_normal_support_review')
        junction_rows = [r for r in rs if r.get('normal_junction_seen') == 'yes']
        # Review C6-03/B3b-08: the rollup treats each match LEVEL by its own
        # evidence strength - exact strand-aware reads can reject, while
        # gene-pair context and strand-ambiguous coordinate hits stay at
        # review unless the operator explicitly keeps the old reject policy.
        exact_junction_rows = [
            r for r in junction_rows
            if r.get('normal_junction_match_level') == 'EXACT_NORMAL_JUNCTION_MATCH'
        ]
        ambiguous_junction_rows = [
            r for r in junction_rows
            if r.get('normal_junction_match_level') == 'POSSIBLE_NORMAL_JUNCTION_STRAND_AMBIGUOUS'
        ]
        gene_pair_junction_rows = [
            r for r in junction_rows
            if r.get('normal_junction_match_level') == 'GENE_PAIR_ONLY_CONTEXT'
        ]
        if exact_junction_rows:
            max_exact_reads = max(to_float(r.get('normal_junction_max_reads'), 0.0) for r in exact_junction_rows)
            if max_exact_reads >= float(safety_cfg.get('normal_junction_reject_reads', 3)):
                worst = worst_safety_status(worst, 'FAIL')
                reasons.append('normal_junction_seen')
            else:
                worst = worst_safety_status(worst, 'CAUTION')
                reasons.append('low_level_normal_junction_seen')
        if gene_pair_junction_rows:
            if gene_pair_action == 'reject':
                worst = worst_safety_status(worst, 'FAIL')
                reasons.append('normal_junction_seen_gene_pair_only')
            else:
                worst = worst_safety_status(worst, 'CAUTION')
                reasons.append('gene_pair_only_normal_junction_context')
        if ambiguous_junction_rows:
            worst = worst_safety_status(worst, 'CAUTION')
            reasons.append('possible_normal_junction_strand_ambiguous')
        missing_layers = sorted({layer for r in rs for layer in r.get('safety_missing_layers','').split(';') if layer})
        event_rows.append({
            'event_id': eid, 'sample_id': e.get('sample_id',''), 'gene': e.get('gene',''), 'event_type': e.get('event_type',''),
            'mutation_source': e.get('mutation_source',''), 'normal_expression_status': rs[0].get('normal_expression_status','UNASSESSED'),
            'normal_junction_status': rs[0].get('normal_junction_assessment_status','UNASSESSED'),
            'normal_junction_match_level': (
                'EXACT_NORMAL_JUNCTION_MATCH' if exact_junction_rows
                else 'POSSIBLE_NORMAL_JUNCTION_STRAND_AMBIGUOUS' if ambiguous_junction_rows
                else 'GENE_PAIR_ONLY_CONTEXT' if junction_rows else ''
            ),
            'matched_normal_status': 'normal_alt_support' if matched_normal else rs[0].get('matched_normal_status',''),
            'normal_hspc_status': rs[0].get('normal_hspc_status','UNASSESSED'),
            'reference_proteome_status': rs[0].get('reference_proteome_status','UNASSESSED'),
            'normal_ligandome_status': rs[0].get('normal_ligandome_status','UNASSESSED'),
            'anchor_assessment_status': 'PARTIAL' if len({r.get('anchor_assessment_status') for r in rs}) > 1 else rs[0].get('anchor_assessment_status','UNASSESSED'),
            'normal_tissue_max_tpm': event_assessment.get('normal_tissue_max_tpm',''),
            'normal_tissue_max_tissue': event_assessment.get('normal_tissue_max_tissue',''),
            'critical_tissue_max_tpm': event_assessment.get('critical_tissue_max_tpm',''),
            'critical_tissue_name': event_assessment.get('critical_tissue_name',''),
            'normal_hspc_tpm': event_assessment.get('normal_hspc_tpm',''),
            'normal_hspc_unit': event_assessment.get('normal_hspc_unit',''),
            'safety_evidence_completeness': f"{min(to_float(r.get('safety_evidence_completeness'),0.0) for r in rs):.4f}",
            'safety_missing_layers': ';'.join(missing_layers),
            'event_safety_status': worst, 'event_safety_reason': combine_safety_reasons(*reasons, 'safety_evidence_incomplete' if missing_layers else '')
        })
    write_tsv(out_peptide_safety, rows, PEPTIDE_SAFETY_FIELDS)
    if out_event_safety:
        write_tsv(out_event_safety, event_rows, EVENT_SAFETY_FIELDS)
    # Review C6-02: record the safety-resource preflight next to the gate
    # outputs so the run manifest and the release-integrity audit can trace
    # it.  Recording only: the gate itself has already represented absent
    # resources as UNASSESSED layers with the priority caps in force.
    from .safety_preflight import run_safety_preflight, write_safety_preflight

    preflight_mode = str(safety_cfg.get('preflight_mode', 'exploratory')).lower()
    if preflight_mode not in {'exploratory', 'release'}:
        preflight_mode = 'exploratory'
    preflight_report = run_safety_preflight(
        paths={
            'reference_proteome': reference_proteome,
            'normal_hla_ligands': normal_hla_ligands,
            'normal_expression': normal_expression,
            'normal_junctions': normal_junctions,
        },
        event_types=[
            f"{e.get('event_type', '')} {e.get('mutation_source', '')}" for e in events.values()
        ],
        candidate_genes={g for e in events.values() for g in _split_genes(str(e.get('gene', '')))},
        mode=preflight_mode,
        waivers=safety_cfg.get('resource_waivers') or [],
    )
    write_safety_preflight(
        preflight_report,
        Path(out_peptide_safety).with_name('safety_preflight.json'),
        Path(out_peptide_safety).with_name('safety_preflight.tsv'),
    )
    return rows, event_rows


def load_peptide_safety(path: str | Path | None) -> dict[str, dict[str, str]]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    return {r.get('peptide_id',''): r for r in read_tsv(p) if r.get('peptide_id')}

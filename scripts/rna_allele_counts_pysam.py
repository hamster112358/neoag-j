#!/usr/bin/env python3
"""Count RNA ref/alt support at somatic VCF sites (pysam pileup).

Review C3-01/02/03 fixes:
  * a missing BAM index (BAI or CSI) aborts the run up front instead of
    surfacing later as swallowed pileup exceptions;
  * VCF-vs-BAM contig intersection is preflighted, so a VCF whose contigs
    are absent from the BAM aborts instead of producing an all-zero table
    that downstream reads as "no expression";
  * review C3-01: contig NAMES alone cannot catch a same-name cross-build
    input (hg19 VCF vs hg38 BAM both use chr1..chrY), so with
    --expected-build the BAM SQ lengths are fingerprinted against exact
    GRCh37/GRCh38 tables; a confirmed mismatch aborts unless
    --allow-build-mismatch records the override in the run summary.
    Without --expected-build the run reports the inferred build
    (BUILD_INFERRED) or BUILD_UNASSESSED;
  * per-site pileup failures are marked ERROR (blank numerics) and never
    silently become zeros;
  * multiallelic sites emit one row per ALT with per-ALT counts, and reads
    carrying neither REF nor a counted ALT stay in the depth as other_reads.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--somatic-vcf", type=Path, required=True)
    p.add_argument("--rna-bam", type=Path, required=True)
    p.add_argument("--output-tsv", type=Path, required=True)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--min-baseq", type=int, default=13)
    p.add_argument("--expected-build", default="",
                   help="expected BAM reference build (GRCh37/hg19 or GRCh38/hg38); "
                        "validated against BAM contig lengths (review C3-01)")
    p.add_argument("--allow-build-mismatch", action="store_true",
                   help="proceed despite a confirmed --expected-build mismatch; "
                        "the override is recorded in the run summary")
    return p.parse_args()


def iter_somatic_sites(vcf: Path):
    """Yield (chrom, pos, ref, alt_list) with every ALT of the site."""
    with open_text(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom, pos, _vid, ref, alt = parts[:5]
            if not alt or alt == ".":
                continue
            alts = [a for a in alt.split(",") if a and a != "."]
            if not alts:
                continue
            yield chrom, int(pos), ref, alts


def normalize_chrom(chrom: str, references: set[str]) -> str:
    if chrom in references:
        return chrom
    if chrom.startswith("chr") and chrom[3:] in references:
        return chrom[3:]
    if f"chr{chrom}" in references:
        return f"chr{chrom}"
    return chrom


# Review C3-01: same-name cross-build guard.  GRCh37 and GRCh38 share the
# chr1..chrY names, so contig overlap proves nothing; exact primary-assembly
# lengths do (UCSC and Ensembl agree on these lengths with or without the
# chr prefix, and hs37d5 shares the GRCh37 primary lengths).
GENOME_BUILD_FINGERPRINTS = {
    "GRCh37": {"1": 249250621, "2": 243199373, "10": 135534747, "17": 81195210, "X": 155270560},
    "GRCh38": {"1": 248956422, "2": 242193529, "10": 133797422, "17": 83257441, "X": 156040895},
}
GENOME_BUILD_ALIASES = {
    "grch37": "GRCh37", "hg19": "GRCh37", "b37": "GRCh37", "hs37d5": "GRCh37",
    "grch38": "GRCh38", "hg38": "GRCh38", "b38": "GRCh38",
}


def bare_chrom(chrom: str) -> str:
    return chrom[3:] if chrom.startswith("chr") else chrom


def fingerprint_build(references, lengths) -> tuple[str, str]:
    """Match BAM SQ contig lengths against the GRCh37/GRCh38 tables.

    Returns (build_or_empty, detail); an empty build means the header gave
    no usable evidence (no fingerprint contigs, or conflicting matches).
    """
    present: dict[str, int] = {}
    for name, length in zip(references, lengths):
        bare = bare_chrom(name)
        if bare in GENOME_BUILD_FINGERPRINTS["GRCh38"] and bare not in present:
            present[bare] = length
    if not present:
        return "", "no fingerprint contigs (1/2/10/17/X) in the BAM header"
    votes: dict[str, list[str]] = {}
    unknown: list[str] = []
    for contig, length in sorted(present.items()):
        hits = [b for b, table in GENOME_BUILD_FINGERPRINTS.items() if table.get(contig) == length]
        if len(hits) == 1:
            votes.setdefault(hits[0], []).append(contig)
        else:
            unknown.append(f"{contig}:{length}")
    if len(votes) > 1:
        return "", f"contigs fingerprint to different builds {votes}"
    if not votes:
        return "", f"no fingerprint contig matches GRCh37 or GRCh38 lengths (e.g. {unknown[:3]})"
    (build, contigs), = votes.items()
    detail = f"confirmed by {len(contigs)}/{len(present)} fingerprint contigs"
    if unknown:
        detail += f"; unmatched {unknown[:3]}"
    return build, detail


def check_genome_build(bam, expected_build: str, allow_mismatch: bool) -> tuple[str, str]:
    """Resolve --expected-build against the BAM length fingerprint (review C3-01).

    Status is BUILD_CONFIRMED / BUILD_MISMATCH / BUILD_INFERRED /
    BUILD_UNASSESSED.  Only a confirmed mismatch against an explicit
    --expected-build aborts, and only without --allow-build-mismatch;
    everything else proceeds with the status recorded in the summary.
    """
    key = (expected_build or "").strip().casefold()
    expected = GENOME_BUILD_ALIASES.get(key, "")
    if key and not expected:
        raise SystemExit(
            f"unknown --expected-build {expected_build!r}; accepted values: "
            "GRCh37/hg19/hs37d5 or GRCh38/hg38 (case-insensitive)"
        )
    inferred, detail = fingerprint_build(bam.references, bam.lengths)
    if expected:
        if inferred == expected:
            return "BUILD_CONFIRMED", f"BAM lengths match {expected} ({detail})"
        if inferred and inferred != expected:
            message = (
                f"--expected-build {expected} but BAM contig lengths fingerprint to "
                f"{inferred} ({detail}); allele counts across builds would be silently "
                "wrong (review C3-01)"
            )
            if not allow_mismatch:
                raise SystemExit(message + "; rerun with --allow-build-mismatch to proceed anyway")
            return "BUILD_MISMATCH", message + " [proceeded via --allow-build-mismatch]"
        return "BUILD_UNASSESSED", f"expected {expected}, but {detail}"
    if inferred:
        return "BUILD_INFERRED", f"BAM contig lengths fingerprint to {inferred} ({detail})"
    return "BUILD_UNASSESSED", detail


def variant_type(ref: str, alt: str) -> str:
    if len(ref) == len(alt) == 1:
        return "SNV"
    if len(alt) > len(ref) and alt.startswith(ref):
        return "INS"
    if len(ref) > len(alt) and ref.startswith(alt):
        return "DEL"
    if len(ref) == len(alt):
        return "MNV"
    return "COMPLEX"


def classify_observation(*, query_sequence: str, query_position: int, indel: int, ref: str, alts: list[str]) -> int | None:
    """Classify one pileup observation at the VCF anchor.

    Returns -1 for reference support, the index i for alts[i] support, or
    None for any other allele (review C3-02: other-allele reads must be kept
    in the depth instead of dropped).
    """
    ref_u = ref.upper()
    for index, alt in enumerate(alts):
        kind = variant_type(ref, alt)
        alt_u = alt.upper()
        if kind == "SNV":
            if query_sequence[query_position].upper() == alt_u:
                return index
        elif kind == "INS":
            expected = len(alt_u) - len(ref_u)
            if indel == expected:
                inserted = query_sequence[query_position + 1 : query_position + 1 + expected].upper()
                if inserted == alt_u[len(ref_u):]:
                    return index
        elif kind == "DEL":
            if indel == -(len(ref_u) - len(alt_u)):
                return index
        elif kind == "MNV" and indel == 0:
            if query_sequence[query_position : query_position + len(ref_u)].upper() == alt_u:
                return index
    if indel == 0 and query_sequence[query_position].upper() == ref_u[0]:
        return -1
    return None


def preflight_checks(bam, vcf_contigs: set[str]) -> None:
    """Fail fast on an unindexed BAM or a VCF whose contigs the BAM lacks.

    Review C3-03: pysam's check_index() accepts either a BAI or a CSI, and a
    review B4/01-era missing index must never surface later as swallowed
    pileup exceptions.  Review C3-01: zero overlapping contigs between the
    VCF and the BAM is a build mismatch, not a biological zero; same-name
    cross-build inputs (GRCh37 vs GRCh38) are caught separately by
    check_genome_build(), because contig names alone cannot distinguish them.
    """
    try:
        bam.check_index()
    except ValueError as exc:  # pysam raises ValueError when neither BAI nor CSI exists
        raise SystemExit(
            f"BAM has no usable index (.bai/.csi): {exc}; refusing to run because "
            "pileup failures would otherwise masquerade as zero coverage (review C3-03)"
        )
    references = set(bam.references)
    mapped = {normalize_chrom(chrom, references) for chrom in vcf_contigs}
    present = mapped & references
    if vcf_contigs and not present:
        raise SystemExit(
            f"none of the {len(vcf_contigs)} VCF contigs "
            f"(e.g. {sorted(vcf_contigs)[:4]}) exist among the BAM references "
            f"(e.g. {sorted(references)[:4]}) - possible coordinate/build mismatch between "
            "VCF and BAM; refusing to emit an all-zero allele table (review C3-01)"
        )


def count_site_alleles(
    bam,
    chrom: str,
    pos: int,
    ref: str,
    alts: list[str],
    *,
    references: set[str],
    min_mapq: int,
    min_baseq: int,
) -> dict:
    """Count one locus. Never converts a failure into zeros (review C3-01):
    pileup errors return status="ERROR" with the message, unknown contigs
    return status="CONTIG_NOT_IN_BAM"."""
    norm = normalize_chrom(chrom, references)
    if norm not in references:
        return {"status": "CONTIG_NOT_IN_BAM", "error": f"contig {chrom} (as {norm}) is not a BAM reference"}
    ref_count = 0
    other_count = 0
    alt_counts = [0] * len(alts)
    try:
        for col in bam.pileup(
            norm,
            pos - 1,
            pos,
            truncate=True,
            max_depth=50_000,
            min_base_quality=min_baseq,
            min_mapping_quality=min_mapq,
        ):
            if col.reference_pos != pos - 1:
                continue
            for pr in col.pileups:
                if pr.is_refskip or pr.is_del or pr.query_position is None:
                    continue
                observation = classify_observation(
                    query_sequence=pr.alignment.query_sequence,
                    query_position=pr.query_position,
                    indel=pr.indel,
                    ref=ref,
                    alts=alts,
                )
                if observation is None:
                    other_count += 1
                elif observation == -1:
                    ref_count += 1
                else:
                    alt_counts[observation] += 1
    except (ValueError, OSError) as exc:
        return {"status": "ERROR", "error": f"{type(exc).__name__}: {exc}"}
    return {
        "status": "ASSESSED",
        "error": "",
        "ref": ref_count,
        "alts": alt_counts,
        "other": other_count,
    }


def build_rows(sites, counts_by_site) -> list[dict[str, str]]:
    """Expand site results to one row per ALT (review C3-02)."""
    rows: list[dict[str, str]] = []
    for (chrom, pos, ref, alts), result in zip(sites, counts_by_site):
        for index, alt in enumerate(alts):
            vtype = variant_type(ref, alt)
            row = {
                "chrom": chrom,
                "pos": str(pos),
                "ref": ref,
                "alt": alt,
                "rna_ref_reads": "",
                "rna_alt_reads": "",
                "rna_other_reads": "",
                "rna_depth": "",
                "rna_vaf": "",
                "variant_key": f"{chrom}:{pos}{ref}>{alt}",
                "variant_type": vtype,
                "count_status": "",
                "rna_count_error": "",
            }
            if result["status"] != "ASSESSED":
                row["count_status"] = result["status"]
                row["rna_count_error"] = result.get("error", "")
            elif vtype == "COMPLEX":
                # Review C3-04: unassessed stays unassessed - blank numerics,
                # not a zero that a consumer could misread as no expression.
                row["count_status"] = "UNASSESSED_COMPLEX"
            else:
                alt_n = result["alts"][index]
                depth = result["ref"] + sum(result["alts"]) + result["other"]
                row["rna_ref_reads"] = str(result["ref"])
                row["rna_alt_reads"] = str(alt_n)
                row["rna_other_reads"] = str(result["other"])
                row["rna_depth"] = str(depth)
                row["rna_vaf"] = f"{alt_n / depth:.4f}" if depth > 0 else ""
                row["count_status"] = "ASSESSED"
            rows.append(row)
    return rows


def main() -> None:
    args = parse_args()
    import pysam

    if not args.rna_bam.is_file():
        raise SystemExit(f"RNA BAM not found: {args.rna_bam}")
    bam = pysam.AlignmentFile(str(args.rna_bam), "rb")
    sites = list(iter_somatic_sites(args.somatic_vcf))
    preflight_checks(bam, {chrom for chrom, _, _, _ in sites})
    build_status, build_detail = check_genome_build(bam, args.expected_build, args.allow_build_mismatch)
    print(f"build_check={build_status} {build_detail}", flush=True)
    references = set(bam.references)
    print(
        f"cross-check vcf_contigs={len({chrom for chrom, _, _, _ in sites})} "
        f"bam_contigs={len(references)} "
        f"matched={len({normalize_chrom(c, references) for c, _, _, _ in sites} & references)}",
        flush=True,
    )
    counts = [
        count_site_alleles(
            bam, chrom, pos, ref, alts,
            references=references, min_mapq=args.min_mapq, min_baseq=args.min_baseq,
        )
        for chrom, pos, ref, alts in sites
    ]
    bam.close()
    rows = build_rows(sites, counts)
    error_rows = [r for r in rows if r["count_status"] in ("ERROR", "CONTIG_NOT_IN_BAM")]
    with_support = sum(1 for r in rows if r["rna_alt_reads"] not in ("", "0"))

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys()) if rows else [
        "chrom", "pos", "ref", "alt", "rna_ref_reads", "rna_alt_reads", "rna_other_reads",
        "rna_depth", "rna_vaf", "variant_key", "variant_type", "count_status", "rna_count_error",
    ]
    with args.output_tsv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(
        f"wrote {args.output_tsv} sites={len(sites)} rows={len(rows)} "
        f"with_alt_support={with_support} unassessed_error_rows={len(error_rows)} "
        f"build_check={build_status}"
    )
    if error_rows:
        print(
            f"WARNING: {len(error_rows)} site rows could not be counted and were left "
            "blank with count_status ERROR/CONTIG_NOT_IN_BAM (review C3-01); "
            "they are NOT measured zeros.",
        )


if __name__ == "__main__":
    main()

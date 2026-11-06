#!/usr/bin/env python3
"""Count RNA ref/alt support at somatic VCF sites (pysam pileup).

Review C3-01/02/03 fixes:
  * a missing BAM index (BAI or CSI) aborts the run up front instead of
    surfacing later as swallowed pileup exceptions;
  * VCF-vs-BAM contig intersection is preflighted, so a cross-build input
    (e.g. hg19 VCF vs hg38 BAM) fails loudly instead of producing an
    all-zero table that downstream reads as "no expression";
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
    VCF and the BAM is a build mismatch, not a biological zero.
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
        f"with_alt_support={with_support} unassessed_error_rows={len(error_rows)}"
    )
    if error_rows:
        print(
            f"WARNING: {len(error_rows)} site rows could not be counted and were left "
            "blank with count_status ERROR/CONTIG_NOT_IN_BAM (review C3-01); "
            "they are NOT measured zeros.",
        )


if __name__ == "__main__":
    main()

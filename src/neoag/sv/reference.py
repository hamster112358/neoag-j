from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from pathlib import Path


def open_text(path: str | Path):
    p = Path(path)
    if str(p).endswith(".gz"):
        return gzip.open(p, "rt", encoding="utf-8", errors="ignore")
    return p.open("r", encoding="utf-8", errors="ignore")


RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1].upper()


class FastaReference:
    """Random-access FASTA reader.

    Review B4-02: the previous implementation loaded the whole genome into RAM
    as Python strings (multi-GB RSS per process).  Regions are now fetched by
    byte-level random access:

    * a samtools ``.fai`` index next to the FASTA is used when present and
      usable (plain FASTA, uniform line layout);
    * otherwise the file is scanned once in streaming mode, recording only byte
      offsets per contig (no sequence is retained);
    * the whole-file load survives only as a fallback for compressed or
      otherwise seek-incompatible FASTAs, where byte offsets are meaningless.

    Normal operation never materialises contig sequences in memory.
    """

    def __init__(self, path: str | Path):
        self.path = Path(path)
        # name -> (length, data_byte_offset, chars_per_line, bytes_per_line)
        self._index: dict[str, tuple[int, int, int, int]] = {}
        self._seqs: dict[str, str] | None = None
        if self.path.suffix == ".gz":
            self._load_whole_file()
            return
        if not self._load_fai():
            self._scan_plain()

    # -- indexing ----------------------------------------------------------
    def _load_fai(self) -> bool:
        """Read a samtools faidx sidecar.  Returns False if absent/unusable."""
        fai = Path(str(self.path) + ".fai")
        if not fai.is_file():
            return False
        try:
            index: dict[str, tuple[int, int, int, int]] = {}
            for line in fai.read_text(encoding="utf-8").splitlines():
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                name = parts[0]
                length, offset, line_len = int(parts[1]), int(parts[2]), int(parts[3])
                line_len_b = int(parts[4]) if len(parts) >= 5 else line_len
                if line_len_b < 0 or line_len < 1 or line_len_b < line_len:
                    # BGZF (line_len_b == -1) or multi-byte line layout: byte
                    # offsets cannot be applied to a plain file read.
                    return False
                chars_per_line = line_len - 1  # line_len includes the newline
                if chars_per_line < 1:
                    return False
                index[name] = (length, offset, chars_per_line, line_len_b)
        except (OSError, ValueError):
            return False
        if not index:
            return False
        self._index = index
        return True

    def _scan_plain(self) -> None:
        """One streaming pass recording per-contig byte offsets; no sequences kept."""
        index: dict[str, tuple[int, int, int, int]] = {}
        try:
            with self.path.open("rb") as fh:
                name: str | None = None
                length = 0
                data_offset = 0
                chars_per_line = 0
                stride = 0

                def flush() -> None:
                    if name is not None and stride > 0 and chars_per_line > 0:
                        index[name] = (length, data_offset, chars_per_line, stride)

                while True:
                    raw = fh.readline()
                    if not raw:
                        break
                    if raw.startswith(b">"):
                        flush()
                        name = raw[1:].split(None, 1)[0].decode("utf-8", "replace").strip()
                        length = 0
                        chars_per_line = 0
                        stride = 0
                        data_offset = fh.tell()
                        continue
                    if not name:
                        continue
                    stripped = raw.rstrip(b"\r\n")
                    if not stripped:
                        continue
                    if chars_per_line == 0:
                        chars_per_line = len(stripped)
                        stride = len(raw)
                    length += len(stripped)
                flush()
        except OSError:
            self._load_whole_file()
            return
        if index:
            self._index = index
        else:
            self._load_whole_file()

    def _load_whole_file(self) -> None:
        seqs: dict[str, list[str]] = {}
        name: str | None = None
        with open_text(self.path) as fh:
            for line in fh:
                if line.startswith(">"):
                    name = line[1:].strip().split()[0]
                    seqs.setdefault(name, [])
                elif name:
                    seqs[name].append(line.strip())
        self._seqs = {k: "".join(v).upper() for k, v in seqs.items()}

    # -- access ------------------------------------------------------------
    def _resolve_contig(self, chrom: str) -> str:
        if self._seqs is not None:
            if chrom in self._seqs:
                return chrom
            alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
            if alt in self._seqs:
                return alt
            raise KeyError(f"Contig not found in FASTA: {chrom}")
        if chrom in self._index:
            return chrom
        alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
        if alt in self._index:
            return alt
        raise KeyError(f"Contig not found in FASTA: {chrom}")

    def contig_length(self, chrom: str) -> int:
        name = self._resolve_contig(chrom)
        if self._seqs is not None:
            return len(self._seqs[name])
        return self._index[name][0]

    def fetch(self, chrom: str, start: int, end: int, strand: str = "+") -> str:
        name = self._resolve_contig(chrom)
        if self._seqs is not None:
            contig = self._seqs[name]
            s = max(1, int(start))
            e = min(len(contig), int(end))
            if e < s:
                return ""
            seq = contig[s - 1 : e]
            return revcomp(seq) if strand == "-" else seq
        length, data_offset, chars_per_line, stride = self._index[name]
        s = max(1, int(start))
        e = min(length, int(end))
        if e < s:
            return ""
        # Byte-range arithmetic: skip only the newline bytes between lines.
        seq = self._read_bases(data_offset, chars_per_line, stride, s - 1, e - s + 1)
        return revcomp(seq) if strand == "-" else seq

    def _read_bases(self, data_offset: int, chars_per_line: int, stride: int, zero_start: int, count: int) -> str:
        line_idx, off_in_line = divmod(zero_start, chars_per_line)
        byte_pos = data_offset + line_idx * stride + off_in_line
        chunks: list[bytes] = []
        remaining = count
        with self.path.open("rb") as fh:
            while remaining > 0:
                take = min(chars_per_line - off_in_line, remaining)
                fh.seek(byte_pos)
                chunk = fh.read(take)
                if not chunk:
                    break
                chunks.append(chunk)
                remaining -= len(chunk)
                byte_pos += len(chunk)
                off_in_line += len(chunk)
                if off_in_line >= chars_per_line:
                    byte_pos += stride - chars_per_line  # skip the line separator
                    off_in_line = 0
        return b"".join(chunks).decode("ascii", "replace").upper()


@dataclass
class Feature:
    chrom: str
    start: int
    end: int
    strand: str
    feature_type: str
    gene_id: str
    gene_name: str
    transcript_id: str = ""
    exon_number: str = ""
    phase: str = "."


@dataclass
class TranscriptModel:
    transcript_id: str
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    cds: list[Feature] = field(default_factory=list)
    exons: list[Feature] = field(default_factory=list)

    def sorted_cds(self) -> list[Feature]:
        if self.strand == "-":
            return sorted(self.cds, key=lambda f: f.start, reverse=True)
        return sorted(self.cds, key=lambda f: f.start)

    def sorted_exons(self) -> list[Feature]:
        if self.strand == "-":
            return sorted(self.exons, key=lambda f: f.start, reverse=True)
        return sorted(self.exons, key=lambda f: f.start)

    def cds_sequence_and_map(self, ref: FastaReference) -> tuple[str, list[int]]:
        seq_parts: list[str] = []
        pos_map: list[int] = []
        for f in self.sorted_cds():
            if self.strand == "-":
                seq_parts.append(ref.fetch(f.chrom, f.start, f.end, "-"))
                pos_map.extend(list(range(f.end, f.start - 1, -1)))
            else:
                seq_parts.append(ref.fetch(f.chrom, f.start, f.end, "+"))
                pos_map.extend(list(range(f.start, f.end + 1)))
        return "".join(seq_parts).upper(), pos_map


def parse_gtf_attributes(attr: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for raw in attr.strip().rstrip(";").split(";"):
        raw = raw.strip()
        if not raw:
            continue
        if " " in raw:
            k, v = raw.split(" ", 1)
            out[k] = v.strip().strip('"')
        elif "=" in raw:
            k, v = raw.split("=", 1)
            out[k] = v.strip().strip('"')
    return out


class GtfAnnotation:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.features: list[Feature] = []
        self.transcripts: dict[str, TranscriptModel] = {}
        self.gene_features: list[Feature] = []
        self._parse()

    def _parse(self) -> None:
        with open_text(self.path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                chrom, _source, ftype, start, end, _score, strand, phase, attrs = parts
                ad = parse_gtf_attributes(attrs)
                gene_id = ad.get("gene_id", "")
                gene_name = ad.get("gene_name", gene_id)
                transcript_id = ad.get("transcript_id", "")
                feat = Feature(
                    chrom=chrom,
                    start=int(start),
                    end=int(end),
                    strand=strand,
                    feature_type=ftype,
                    gene_id=gene_id,
                    gene_name=gene_name,
                    transcript_id=transcript_id,
                    exon_number=ad.get("exon_number", ""),
                    phase=phase,
                )
                if ftype == "gene":
                    self.gene_features.append(feat)
                if ftype in {"CDS", "exon"} and transcript_id:
                    tx = self.transcripts.get(transcript_id)
                    if not tx:
                        tx = TranscriptModel(transcript_id, gene_id, gene_name, chrom, strand)
                        self.transcripts[transcript_id] = tx
                    if ftype == "CDS":
                        tx.cds.append(feat)
                    elif ftype == "exon":
                        tx.exons.append(feat)
                    self.features.append(feat)

    def features_at(self, chrom: str, pos: int, feature_types: set[str] | None = None) -> list[Feature]:
        chroms = {chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom}
        out = []
        for f in self.features:
            if f.chrom in chroms and f.start <= pos <= f.end and (feature_types is None or f.feature_type in feature_types):
                out.append(f)
        return out

    def gene_at(self, chrom: str, pos: int) -> str:
        feats = self.features_at(chrom, pos, {"CDS", "exon"})
        if feats:
            cds = [f for f in feats if f.feature_type == "CDS"]
            chosen = (cds or feats)[0]
            return chosen.gene_name or chosen.gene_id
        chroms = {chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom}
        for f in self.gene_features:
            if f.chrom in chroms and f.start <= pos <= f.end:
                return f.gene_name or f.gene_id
        return ""

    def best_transcript_at(self, chrom: str, pos: int) -> TranscriptModel | None:
        feats = self.features_at(chrom, pos, {"CDS"}) or self.features_at(chrom, pos, {"exon"})
        if not feats:
            return None
        tids = {f.transcript_id for f in feats if f.transcript_id in self.transcripts}
        if not tids:
            return None
        # Prefer the transcript with the longest CDS.
        return sorted((self.transcripts[t] for t in tids), key=lambda tx: sum(f.end - f.start + 1 for f in tx.cds), reverse=True)[0]

    @staticmethod
    def _tx_cds_len(tx: TranscriptModel) -> int:
        return sum(f.end - f.start + 1 for f in tx.cds)

    @staticmethod
    def _tx_distance(tx: TranscriptModel, pos: int) -> int:
        spans = tx.exons or tx.cds
        if not spans:
            return 10**12
        if any(f.start <= pos <= f.end for f in spans):
            return 0
        return min(abs(pos - f.start) if pos < f.start else abs(pos - f.end) for f in spans)

    def best_coding_transcript_near(self, chrom: str, pos: int, max_distance_bp: int = 200_000) -> tuple[TranscriptModel | None, int, str]:
        """Return a coding transcript at, within the same gene as, or near a breakpoint.

        Long-read SV breakpoints often fall in introns rather than CDS/exons. For
        first-pass neoantigen triage, mapping to the closest coding transcript of
        the same gene preserves candidate fusions while keeping the reconstruction
        explicitly heuristic.
        """
        exact = self.best_transcript_at(chrom, pos)
        if exact:
            return exact, 0, "overlapping_cds_or_exon"
        chroms = {chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom}
        gene = self.gene_at(chrom, pos)
        candidates = [
            tx for tx in self.transcripts.values()
            if tx.chrom in chroms and tx.cds and (not gene or tx.gene_name == gene or tx.gene_id == gene)
        ]
        if candidates:
            best = sorted(candidates, key=lambda tx: (self._tx_distance(tx, pos), -self._tx_cds_len(tx)))[0]
            return best, self._tx_distance(best, pos), "same_gene_nearest_coding_transcript"
        nearby = [
            tx for tx in self.transcripts.values()
            if tx.chrom in chroms and tx.cds and self._tx_distance(tx, pos) <= max_distance_bp
        ]
        if not nearby:
            return None, 0, "no_nearby_coding_transcript"
        best = sorted(nearby, key=lambda tx: (self._tx_distance(tx, pos), -self._tx_cds_len(tx)))[0]
        return best, self._tx_distance(best, pos), "nearby_coding_transcript"

    def transcript_by_id(self, transcript_id: str) -> TranscriptModel | None:
        return self.transcripts.get(transcript_id)

    def are_potential_readthrough_partners(
        self,
        chrom1: str,
        gene1: str,
        chrom2: str,
        gene2: str,
        max_gap_bp: int = 100_000,
    ) -> bool:
        """True when two distinct genes sit head-to-tail next to each other.

        Review B4-13: such a pair can produce a NATIVE read-through splice
        junction whose coordinates a projected DNA fusion junction can match
        exactly, without any tumour genomic rearrangement.  The check is a
        conservative candidate flag (same contig, same strand, small gap), not
        a classification.
        """
        from .identity import normalize_contig

        if not gene1 or not gene2 or gene1 == gene2:
            return False
        if normalize_contig(chrom1) != normalize_contig(chrom2):
            return False
        contigs = {normalize_contig(chrom1)}

        def _spans(gene: str) -> list[Feature]:
            return [
                f for f in self.gene_features
                if (f.gene_name or f.gene_id) == gene and normalize_contig(f.chrom) in contigs
            ]

        for a in _spans(gene1):
            for b in _spans(gene2):
                if a.strand != b.strand:
                    continue
                if a.end < b.start and b.start - a.end - 1 <= max_gap_bp:
                    return True
                if b.end < a.start and a.start - b.end - 1 <= max_gap_bp:
                    return True
        return False

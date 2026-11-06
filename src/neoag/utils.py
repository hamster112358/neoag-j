from __future__ import annotations
import csv, gzip, json, math, re, shutil
from pathlib import Path
from typing import Any, Iterable, Mapping, TextIO

MISSING = {"", "NA", "NaN", "nan", "None", ".", "null", "NULL"}


def is_gzip_file(path: str | Path) -> bool:
    """True when file begins with gzip magic (handles misnamed plain-text .gz)."""
    p = Path(path)
    if not p.is_file():
        return False
    with p.open("rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def open_text_maybe_gz(path: str | Path, *, encoding: str = "utf-8") -> TextIO:
    """Open text file, using gzip only when content is actually gzipped."""
    p = Path(path)
    if is_gzip_file(p):
        return gzip.open(p, "rt", encoding=encoding)
    return p.open(encoding=encoding)


def copy_if_different(src: str | Path, dst: str | Path) -> Path:
    """Copy src to dst, skipping the copy when both paths already name the same file."""
    src_path = Path(src)
    dst_path = Path(dst)
    dst_path.parent.mkdir(parents=True, exist_ok=True)
    if dst_path.exists() and src_path.resolve() == dst_path.resolve():
        return dst_path
    shutil.copy2(src_path, dst_path)
    return dst_path


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Missing TSV: {p}")
    with p.open("r", encoding="utf-8", newline="") as fh:
        return [dict(row) for row in csv.DictReader(fh, delimiter="\t")]

def read_csv(path: str | Path) -> list[dict[str, str]]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Missing CSV: {p}")
    with p.open("r", encoding="utf-8", newline="") as fh:
        return [dict(row) for row in csv.DictReader(fh)]

def write_tsv(path: str | Path, rows: Iterable[Mapping[str, Any]], fieldnames: list[str] | None = None) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if fieldnames is None:
        fieldnames = []
        for row in rows:
            for k in row.keys():
                if k not in fieldnames:
                    fieldnames.append(k)
    with p.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: "" if row.get(k) is None else row.get(k) for k in fieldnames})

def write_json(path: str | Path, obj: Any) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")

def to_float(value: Any, default: float = 0.0) -> float:
    if value is None:
        return default
    if isinstance(value, str) and value.strip() in MISSING:
        return default
    try:
        x = float(value)
        if math.isnan(x) or math.isinf(x):
            return default
        return x
    except Exception:
        return default

def clamp(x: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, x))

def norm_rank(rank: Any, cap: float = 10.0) -> float:
    r = to_float(rank, cap)
    return clamp(1.0 - min(r, cap) / cap)

def norm_tpm(tpm: Any, high: float = 20.0) -> float:
    return clamp(to_float(tpm, 0.0) / high)

def truthy(v: Any) -> bool:
    if isinstance(v, bool):
        return v
    return str(v).strip().lower() in {"1", "true", "yes", "y", "pass", "present", "positive"}

def safe_id(text: Any) -> str:
    s = str(text if text is not None else "")
    s = "".join(ch if ch.isalnum() else "_" for ch in s)
    return "_".join(x for x in s.split("_") if x)[:160]

def first(row: Mapping[str, Any], keys: list[str], default: str = "") -> str:
    lower = {k.lower(): k for k in row.keys() if isinstance(k, str)}
    for k in keys:
        if k in row and str(row[k]).strip() not in MISSING:
            return str(row[k])
        lk = k.lower()
        if lk in lower and str(row[lower[lk]]).strip() not in MISSING:
            return str(row[lower[lk]])
    return default

def split_ws_or_tab(line: str) -> list[str]:
    if "\t" in line:
        return [x.strip() for x in line.strip().split("\t")]
    return [x.strip() for x in re.split(r"\s+", line.strip()) if x.strip()]


def toml_literal(value: Any) -> str:
    """Serialize one TOML scalar/array with full escaping (review B1-02).

    JSON string escaping is a safe subset of TOML basic strings, so quotes,
    backslashes and spaces in paths survive a tomllib round-trip instead of
    silently breaking the document structure.
    """
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return repr(value)
    if isinstance(value, (list, tuple)):
        return "[" + ", ".join(toml_literal(v) for v in value) + "]"
    return json.dumps(str(value), ensure_ascii=False)


def write_toml(path: str | Path, doc: Mapping[str, Any]) -> Path:
    """Write a flat-table TOML document and verify it parses before saving.

    Review B1-02: runtime configs were built by interpolating paths into a
    hand-written template, where one quote or backslash in a path silently
    changed the document's meaning.  Values are emitted from typed data and
    the rendered text must round-trip through tomllib, so a broken document
    fails here instead of mid-run.
    """
    import tomllib

    lines: list[str] = []
    scalars: list[str] = []
    tables: list[tuple[str, Mapping[str, Any]]] = []
    for key, value in doc.items():
        if isinstance(value, Mapping):
            tables.append((str(key), value))
        else:
            if value is None:
                continue  # absent is honest; TOML has no null
            scalars.append(f"{key} = {toml_literal(value)}")
    lines.extend(scalars)
    for name, table in tables:
        lines.append("")
        lines.append(f"[{name}]")
        for key, value in table.items():
            if isinstance(value, Mapping) or value is None:
                continue
            lines.append(f"{key} = {toml_literal(value)}")
    text = "\n".join(lines).strip() + "\n"
    parsed = tomllib.loads(text)  # raises on any malformed emission
    for name, table in tables:
        for key in table:
            if table[key] is not None and key not in parsed.get(name, {}):
                raise ValueError(f"write_toml: key {name}.{key} lost in round-trip")
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")
    return p

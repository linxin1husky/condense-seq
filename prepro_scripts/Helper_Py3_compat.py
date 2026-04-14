from __future__ import annotations

import bz2
from dataclasses import dataclass
import gzip
import lzma
import re
import subprocess
import sys
from typing import Iterable, Iterator, Literal, Optional, TextIO, Tuple, overload

"""Helper utilities for Condense-seq preprocessing scripts.

Python 3.11+ compatible. Keeps legacy helpers while adding reusable utilities.
"""

_COMPRESSED_SUFFIXES = {
    ".gz": gzip.open,
    ".bz2": bz2.open,
    ".xz": lzma.open,
}


def open_any(path: str, mode: str = "rt", encoding: str = "utf-8", newline: Optional[str] = None) -> TextIO:
    """Open a text/binary file, handling .gz/.bz2/.xz files automatically."""
    for suffix, opener in _COMPRESSED_SUFFIXES.items():
        if path.endswith(suffix):
            if "b" in mode:
                return opener(path, mode)
            return opener(path, mode, encoding=encoding, newline=newline)
    if "b" in mode:
        raise ValueError(
            "Binary mode requires a supported compressed suffix: "
            + ", ".join(sorted(_COMPRESSED_SUFFIXES.keys()))
            + f". Got: {path}. The current implementation only supports"
            + ".gz/.bz2/.xz for binary files."
        )
    return open(path, mode, encoding=encoding, newline=newline)


def gzopen(fname: str) -> TextIO:
    """Legacy wrapper for backward compatibility. Prefer open_any."""
    return open_any(fname, "rt")


def iter_lines(path: str, skip_empty: bool = True, comment_prefix: Optional[str] = None) -> Iterator[str]:
    """Iterate over lines, optionally skipping empty/comment lines."""
    with open_any(path, "rt") as f:
        for line in f:
            line = line.strip()
            if skip_empty and not line:
                continue
            if comment_prefix and line.startswith(comment_prefix):
                continue
            yield line


def chr_key(name: str) -> Tuple[int, object]:
    """Sorting key for chromosome names like chr1/1/chrX/MT."""
    if name is None:
        return (2, "")
    s = name.strip()
    if s.lower().startswith("chr"):
        s = s[3:]
    s_upper = s.upper()
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if s_upper.isdigit():
        return (0, int(s_upper))
    if s_upper in special:
        return (0, special[s_upper])
    return (1, s_upper)


def rev_cmp(seq: str) -> str:
    """Reverse complement of a DNA sequence with IUPAC handling."""
    if not seq:
        return ""
    trans = str.maketrans(
        "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
        "TGCAYRSWMKVHDBNtgacayrswmkvhdbn",
    )
    return seq.translate(trans)[::-1]


def parse_fasta(path: str, upper: bool = True) -> Iterator[Tuple[str, str]]:
    """Yield (name, seq) records from a FASTA file."""
    name = None
    seq_chunks = []
    for line in iter_lines(path, skip_empty=True, comment_prefix=None):
        if line.startswith(">"):
            if name is not None:
                seq = "".join(seq_chunks)
                yield (name, seq.upper() if upper else seq)
            name = line[1:].split()[0]
            seq_chunks = []
        else:
            seq_chunks.append(line)
    if name is not None:
        seq = "".join(seq_chunks)
        yield (name, seq.upper() if upper else seq)


def genome_sizes(path: str) -> dict:
    """Return a dict of contig sizes from a FASTA file."""
    sizes = {}
    for name, seq in parse_fasta(path, upper=False):
        sizes[name] = len(seq)
    return sizes


def gc_content(seq: str) -> float:
    """Compute GC content of a sequence."""
    if not seq:
        return 0.0
    seq = seq.upper()
    gc = sum(1 for nt in seq if nt in ("G", "C"))
    return float(gc) / float(len(seq))


def bin_gc_from_fasta(path: str, win_size: int) -> dict:
    """Compute GC content for each genome bin of size win_size."""
    if win_size <= 0:
        raise ValueError("win_size must be positive")
    out = {}
    for name, seq in parse_fasta(path, upper=True):
        num = 0
        for i in range(0, len(seq), win_size):
            chunk = seq[i : i + win_size]
            if not chunk:
                continue
            out[f"{name}_{num}"] = gc_content(chunk)
            num += 1
    return out


def read_gtab_header(cols: Iterable[str]) -> Tuple[str, int, int, bool]:
    """Parse GTAB-like header and return (data_type, col_st, col_ed, has_gc)."""
    cols = list(cols)
    if "Chromosome" not in cols:
        raise ValueError("Header missing 'Chromosome' column")

    has_gc = "GCcontent" in cols
    col_ed = cols.index("GCcontent") if has_gc else len(cols)

    if "Position" in cols or "PhysicalPosition" in cols:
        data_type = "point"
        pos_key = "Position" if "Position" in cols else "PhysicalPosition"
        col_st = cols.index(pos_key) + 1
    elif "Start" in cols and "End" in cols:
        data_type = "binned"
        col_st = cols.index("End") + 1
    else:
        raise ValueError("Unrecognized GTAB header format")

    return data_type, col_st, col_ed, has_gc


def parse_gtab_row(cols: Iterable[str], col_st: int) -> Tuple[str, int, int, list[float]]:
    """Parse a GTAB-like row into chromosome, interval, and numeric values."""
    cols = list(cols)
    if col_st == 2:
        chrom = cols[0]
        start = int(cols[1])
        end = start + 1
    else:
        chrom = cols[0]
        start = int(cols[1])
        end = int(cols[2])
    values = [float(value) for value in cols[col_st:]]
    return chrom, start, end, values


def iter_gtab_records(path: str) -> Iterator[Tuple[list[str], int, str, int, int, list[float]]]:
    """Iterate GTAB-like records, yielding header labels plus parsed row data."""
    labels: Optional[list[str]] = None
    col_st: Optional[int] = None
    for line in open_any(path, "rt"):
        cols = line.strip().split()
        if not cols:
            continue
        if labels is None or col_st is None:
            _, col_st, _, _ = read_gtab_header(cols)
            labels = cols[col_st:]
            continue
        chrom, start, end, values = parse_gtab_row(cols, col_st)
        yield labels, col_st, chrom, start, end, values


class BinLayout:
    """Reusable fixed-window bin layout for interval overlap iteration."""

    def __init__(self, bin_size: int, bin_step: int):
        if bin_size <= 0:
            raise ValueError("bin_size must be positive")
        if bin_step <= 0:
            raise ValueError("bin_step must be positive")
        self.bin_size = int(bin_size)
        self.bin_step = int(bin_step)

    def iter_overlaps(self, start: int, end: int) -> Iterator[Tuple[int, int]]:
        if end <= start:
            return

        idx = start // self.bin_step
        bin_start = self.bin_step * idx
        bin_end = bin_start + self.bin_size

        idx_start = idx
        while bin_start <= start < bin_end:
            if idx < idx_start:
                idx_start = idx
            idx -= 1
            if idx < 0:
                break
            bin_start = self.bin_step * idx
            bin_end = bin_start + self.bin_size

        idx_end = (end - 1) // self.bin_step
        if idx_start > idx_end:
            raise ValueError("invalid overlap interval")

        for idx in range(idx_start, idx_end + 1):
            bin_start = self.bin_step * idx
            bin_end = bin_start + self.bin_size
            left = max(start, bin_start)
            right = min(end, bin_end)
            overlap_len = right - left
            if overlap_len > 0:
                yield idx, overlap_len


@dataclass(frozen=True)
class SamRecord:
    qname: str
    flag: int
    rname: str
    pos: int
    mapq: int
    cigar: str
    rnext: str
    pnext: int
    tlen: int
    seq: str
    qual: str
    tags: dict[str, str | int]


def iter_samtools_view(
    path: str,
    samtools: str = "samtools",
    extra_args: Optional[list[str]] = None,
    include_header: bool = False,
) -> Iterator[list[str]]:
    """Yield SAM columns from `samtools view` without using shell=True."""
    cmd = [samtools, "view"]
    if include_header:
        cmd.append("-h")
    if extra_args:
        cmd.extend(extra_args)
    cmd.append(path)

    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(f"Required executable not found: {samtools}") from exc

    assert proc.stdout is not None
    try:
        for line in proc.stdout:
            if not include_header and line.startswith("@"):
                continue
            yield line.rstrip("\n").split("\t")
    finally:
        proc.stdout.close()
        stderr = proc.stderr.read() if proc.stderr is not None else ""
        returncode = proc.wait()
        if proc.stderr is not None:
            proc.stderr.close()
        if returncode != 0:
            raise RuntimeError(
                f"{samtools} view failed for {path} with exit code {returncode}: "
                + stderr.strip()
            )


def parse_sam_tags(cols: list[str]) -> dict[str, str | int]:
    """Parse optional SAM tags into a dictionary."""
    tags: dict[str, str | int] = {}
    for col in cols[11:]:
        parts = col.split(":", 2)
        if len(parts) != 3:
            continue
        tag, tag_type, value = parts
        if tag_type == "i":
            try:
                tags[tag] = int(value)
                continue
            except ValueError:
                pass
        tags[tag] = value
    return tags


def parse_sam_record(cols: list[str]) -> SamRecord:
    """Convert a SAM row into a structured SamRecord."""
    if len(cols) < 11:
        raise ValueError("SAM record must contain at least 11 columns")
    return SamRecord(
        qname=cols[0],
        flag=int(cols[1]),
        rname=cols[2].strip(),
        pos=int(cols[3]),
        mapq=int(cols[4]),
        cigar=cols[5],
        rnext=cols[6],
        pnext=int(cols[7]),
        tlen=int(cols[8]),
        seq=cols[9],
        qual=cols[10],
        tags=parse_sam_tags(cols),
    )


def cigar_reference_span(cigar: str) -> int:
    """Return reference-consuming span implied by a CIGAR string."""
    if cigar == "*":
        return 0
    span = 0
    parts = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    if not parts or "".join(num + op for num, op in parts) != cigar:
        raise ValueError(f"Malformed CIGAR string: {cigar}")
    for count_str, op in parts:
        count = int(count_str)
        if op in {"M", "D", "N", "=", "X"}:
            span += count
    return span


def fragment_bounds_from_record(record: SamRecord) -> tuple[int, int]:
    """Return half-open fragment bounds for a paired-end SAM record."""
    pos0 = record.pos - 1
    if record.tlen == 0:
        raise ValueError("Cannot derive fragment bounds for zero template length")
    if record.tlen > 0:
        start = pos0
        end = pos0 + record.tlen
        return start, end

    end = pos0 + cigar_reference_span(record.cigar)
    start = end + record.tlen
    return start, end


def passes_common_paired_filters(
    record: SamRecord,
    chr_choices: Optional[set[str]] = None,
    min_len: Optional[int] = None,
    max_len: Optional[int] = None,
    mm_cutoff: Optional[int] = None,
    require_proper_pair: bool = True,
) -> bool:
    """Return True when the read passes shared paired-end filtering rules."""
    pos0 = record.pos - 1
    if chr_choices and record.rname not in chr_choices:
        return False
    if pos0 < 0 or record.flag & 0x4 != 0 or record.rname == "*":
        return False
    if record.flag & 0x8 != 0:
        return False
    if require_proper_pair and record.flag & 0x2 == 0:
        return False

    frag_len = abs(record.tlen)
    if min_len is not None and frag_len < min_len:
        return False
    if max_len is not None and frag_len > max_len:
        return False

    if mm_cutoff is not None:
        nm = record.tags.get("NM")
        if isinstance(nm, int) and nm > mm_cutoff:
            return False
    return True


def fragment_center_positions(
    start: int,
    end: int,
    even_policy: str = "split",
) -> list[tuple[int, float]]:
    """Return center positions and weights for a half-open fragment interval."""
    if end <= start:
        raise ValueError("Fragment end must be greater than start")
    length = end - start
    if length % 2 != 0:
        return [(start + length // 2, 1.0)]
    if even_policy == "split":
        mid_right = start + length // 2
        return [(mid_right - 1, 0.5), (mid_right, 0.5)]
    if even_policy == "left":
        return [(start + length // 2 - 1, 1.0)]
    raise ValueError(f"Unsupported even_policy: {even_policy}")


def record_center_positions(
    record: SamRecord,
    even_policy: str = "split",
) -> list[tuple[int, float]]:
    """Return legacy center positions directly from a paired-end record."""
    pos0 = record.pos - 1
    tlen = record.tlen
    if tlen == 0:
        raise ValueError("Cannot derive center positions for zero template length")

    if tlen > 0:
        if tlen % 2 != 0:
            return [(pos0 + tlen // 2, 1.0)]
        if even_policy == "split":
            mid_right = pos0 + tlen // 2
            return [(mid_right - 1, 0.5), (mid_right, 0.5)]
        if even_policy == "left":
            return [(pos0 + tlen // 2 - 1, 1.0)]
        raise ValueError(f"Unsupported even_policy: {even_policy}")

    end_pos = pos0 + cigar_reference_span(record.cigar)
    if tlen % 2 != 0:
        return [(end_pos + tlen // 2 - 1, 1.0)]
    if even_policy == "split":
        mid_right = end_pos + tlen // 2
        return [(mid_right - 1, 0.5), (mid_right, 0.5)]
    if even_policy == "left":
        return [(end_pos + tlen // 2 - 1, 1.0)]
    raise ValueError(f"Unsupported even_policy: {even_policy}")


def str2bool(v: str) -> bool:
    """Standard CLI bool parser."""
    v = v.lower()
    if v in ("yes", "true", "t", "y", "1"):
        return True
    if v in ("no", "false", "f", "n", "0"):
        return False
    raise ValueError("Boolean value expected.")


@overload
def read_titration(
    fname: str,
    return_conc: Literal[False] = False,
    conc_col: int = 0,
    frac_col: int = -3,
    tnum_col: int = -1,
    delim: Optional[str] = None,
) -> dict[int, float]:
    ...


@overload
def read_titration(
    fname: str,
    return_conc: Literal[True] = True,
    conc_col: int = 0,
    frac_col: int = -3,
    tnum_col: int = -1,
    delim: Optional[str] = None,
) -> Tuple[dict[int, float], dict[int, float]]:
    ...


def read_titration(
    fname: str,
    return_conc: bool = False,
    conc_col: int = 0,
    frac_col: int = -3,
    tnum_col: int = -1,
    delim: Optional[str] = None,
) -> dict[int, float] | Tuple[dict[int, float], dict[int, float]]:
    """Read titration file; optionally return concentration and fraction maps."""
    tnum_tfrac = {}
    tnum_conc = {}
    first = True
    with open_any(fname, "rt") as f:
        for line in f:
            if first:
                first = False
                continue
            line = line.strip()
            if not line:
                continue
            cols = line.split(delim) if delim is not None else line.split()
            try:
                tnum = int(cols[tnum_col])
            except Exception:
                continue
            total_frac = float(cols[frac_col])
            if tnum in tnum_tfrac:
                raise ValueError(f"Duplicate titration number: {tnum}")
            tnum_tfrac[tnum] = total_frac
            if return_conc:
                tnum_conc[tnum] = float(cols[conc_col])
    return (tnum_conc, tnum_tfrac) if return_conc else tnum_tfrac


if __name__ == "__main__":
    sys.stderr.write(
        "Warning: This script is a helper module and is not intended to be run directly.\n"
        "Please run the main pipeline script instead.\n"
    )
    sys.exit(1)

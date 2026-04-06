from __future__ import annotations

import bz2
import gzip
import lzma
import sys
from typing import Iterable, Iterator, Optional, TextIO, Tuple

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


def str2bool(v: str) -> bool:
    """Standard CLI bool parser."""
    v = v.lower()
    if v in ("yes", "true", "t", "y", "1"):
        return True
    if v in ("no", "false", "f", "n", "0"):
        return False
    raise ValueError("Boolean value expected.")


def read_titration(
    fname: str,
    return_conc: bool = False,
    conc_col: int = 0,
    frac_col: int = -3,
    tnum_col: int = -1,
    delim: Optional[str] = None,
) -> dict:
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

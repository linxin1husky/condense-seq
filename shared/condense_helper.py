from __future__ import annotations

import bz2
import gzip
import lzma
import re
from typing import Iterable, Iterator, TextIO


"""Generic helpers shared by preprocessing and postprocessing scripts."""

_COMPRESSED_SUFFIXES = {
    ".gz": gzip.open,
    ".bz2": bz2.open,
    ".xz": lzma.open,
}


def open_text(
    path: str,
    mode: str = "rt",
    encoding: str = "utf-8",
    newline: str | None = None,
) -> TextIO:
    """Open plain or compressed text files with consistent defaults."""
    if "b" in mode:
        raise ValueError("open_text only supports text modes; use a binary opener for binary data")
    for suffix, opener in _COMPRESSED_SUFFIXES.items():
        if path.endswith(suffix):
            return opener(path, mode, encoding=encoding, newline=newline)
    return open(path, mode, encoding=encoding, newline=newline)


def iter_nonempty_lines(path: str, comment_prefix: str | None = None) -> Iterator[str]:
    """Yield stripped non-empty lines, optionally skipping comments."""
    with open_text(path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if comment_prefix is not None and line.startswith(comment_prefix):
                continue
            yield line


def coerce_value(value: str, na_to_nan: bool = True):
    """Convert numeric strings to float and optionally NA-like tokens to np.nan."""
    if na_to_nan and value in {"", "NA", "NaN", "nan"}:
        return float("nan")
    try:
        return float(value)
    except ValueError:
        return value


def to_column_view(row_dict: dict) -> dict:
    """Convert {row_id: {field: value}} into {field: {row_id: value}}."""
    field_dict = {}
    for row_id, field_values in row_dict.items():
        for field, value in field_values.items():
            if field not in field_dict:
                field_dict[field] = {}
            field_dict[field][row_id] = value
    return field_dict


def chr_key(name: str):
    """Sorting key for chromosome names like chr1/1/chrX/MT."""
    if name is None:
        return (2, "")
    s = str(name).strip()
    if s.lower().startswith("chr"):
        s = s[3:]
    s_upper = s.upper()
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if s_upper.isdigit():
        return (0, int(s_upper))
    if s_upper in special:
        return (0, special[s_upper])
    return (1, s_upper)


def parse_gtab_header(cols: Iterable[str]) -> tuple[str, int, int, bool]:
    """Parse a GTAB-like header and return (data_type, col_st, col_ed, has_gc)."""
    cols = list(cols)
    if "Chromosome" not in cols:
        raise ValueError("Header missing 'Chromosome' column")

    has_gc = "GCcontent" in cols
    col_ed = cols.index("GCcontent") if has_gc else len(cols)

    if "Position" in cols or "PhysicalPosition" in cols:
        pos_key = "Position" if "Position" in cols else "PhysicalPosition"
        return "point", cols.index(pos_key) + 1, col_ed, has_gc
    if "Start" in cols and "End" in cols:
        return "binned", cols.index("End") + 1, col_ed, has_gc
    raise ValueError("Unrecognized GTAB header format")


def parse_gtf_attributes(attribute: str) -> dict[str, str]:
    """Parse the GTF attribute column into a tag-value dictionary."""
    tag_value = {}
    for item in attribute.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if " " in item:
            tag, value = item.split(None, 1)
        elif "=" in item:
            tag, value = item.split("=", 1)
        else:
            continue
        tag_value[tag] = value.strip().strip('"')
    return tag_value


def merge_intervals(intervals: list[list[int]]) -> list[list[int]]:
    """Merge overlapping intervals represented as [start, end]."""
    merged = []
    for start, end in sorted(intervals):
        if merged and merged[-1][1] >= start:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged

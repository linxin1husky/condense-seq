from __future__ import annotations

import math
import sys

__all__ = [
    "RegularGenomeBins",
    "GenomicIntervalIndex",
    "bin_hash",
    "double_hash",
    "binary_search",
]


def _binary_search_left(sorted_values, target):
    """Return the leftmost insertion index for target in sorted_values."""
    st, ed = 0, len(sorted_values)
    while st < ed:
        mid = (st + ed) // 2
        if sorted_values[mid] < target:
            st = mid + 1
        else:
            ed = mid
    return st


def binary_search(sortlist, target):
    """Legacy public name retained for compatibility."""
    return _binary_search_left(sortlist, target)


def _overlap_len(start1, end1, start2, end2):
    return max(0, min(end1, end2) - max(start1, start2))


def _infer_max_pos(intervals):
    if not intervals:
        raise ValueError("Cannot infer max_pos from an empty interval mapping")
    return max(end for _, end in intervals.values())


def _infer_domain_size(max_pos):
    if max_pos <= 0:
        return 1
    # Match the old intent: split large genomic coordinates into coarse domains.
    return max(1, int(10 ** (int(math.log10(max_pos)) // 2)))


def _is_interval_mapping(value):
    if value is None or not isinstance(value, dict):
        return False
    if not value:
        return True
    first_value = next(iter(value.values()))
    return isinstance(first_value, (tuple, list)) and len(first_value) == 2


def _normalize_bin_constructor_args(args, kwargs):
    kwargs = dict(kwargs)
    silent = kwargs.pop("silent", False)
    ID_interval_kw = kwargs.pop("ID_interval", None)
    max_pos_kw = kwargs.pop("max_pos", None)
    bin_size_kw = kwargs.pop("bin_size", None)
    bin_step_kw = kwargs.pop("bin_step", None)
    if kwargs:
        unknown = ", ".join(sorted(kwargs))
        raise TypeError(f"Unexpected keyword argument(s): {unknown}")

    ID_interval = ID_interval_kw
    bin_size = bin_size_kw
    bin_step = bin_step_kw
    max_pos = max_pos_kw

    if args and _is_interval_mapping(args[0]):
        # Python3/current style: bin_hash(ID_interval, bin_size, bin_step, max_pos)
        if ID_interval is not None:
            raise TypeError("ID_interval provided both positionally and by keyword")
        ID_interval = args[0]
        rest = args[1:]
        if len(rest) > 3:
            raise TypeError("bin_hash accepts at most 4 positional arguments in ID_interval-first form")
        if len(rest) >= 1:
            if bin_size is not None:
                raise TypeError("bin_size provided both positionally and by keyword")
            bin_size = rest[0]
        if len(rest) >= 2:
            if bin_step is not None:
                raise TypeError("bin_step provided both positionally and by keyword")
            bin_step = rest[1]
        if len(rest) >= 3:
            if max_pos is not None:
                raise TypeError("max_pos provided both positionally and by keyword")
            max_pos = rest[2]
    else:
        # Legacy style: bin_hash(bin_size, bin_step, max_pos=None, ID_interval=None)
        if len(args) > 3:
            raise TypeError("bin_hash accepts at most 3 positional arguments in legacy form")
        if len(args) >= 1:
            if bin_size is not None:
                raise TypeError("bin_size provided both positionally and by keyword")
            bin_size = args[0]
        if len(args) >= 2:
            if bin_step is not None:
                raise TypeError("bin_step provided both positionally and by keyword")
            bin_step = args[1]
        if len(args) >= 3:
            if max_pos is not None:
                raise TypeError("max_pos provided both positionally and by keyword")
            max_pos = args[2]

    if bin_size is None or bin_step is None:
        raise TypeError("bin_hash requires bin_size and bin_step")
    bin_size = int(bin_size)
    bin_step = int(bin_step)
    if bin_size <= 0 or bin_step <= 0:
        raise ValueError("bin_size and bin_step must be positive")

    if ID_interval is None:
        if max_pos is None:
            raise ValueError("max_pos must be provided when ID_interval is None")
        max_pos = int(max_pos)
        n_bins = int(max_pos // bin_step) + 1
        ID_interval = {i: (i * bin_step, i * bin_step + bin_size) for i in range(n_bins)}
    else:
        ID_interval = {ID: (int(st), int(ed)) for ID, (st, ed) in ID_interval.items()}
        if max_pos is None:
            max_pos = _infer_max_pos(ID_interval)
        max_pos = int(max_pos)

    return ID_interval, bin_size, bin_step, max_pos, silent


def _normalize_interval_constructor_args(args, kwargs):
    kwargs = dict(kwargs)
    silent = kwargs.pop("silent", False)
    domain_size = kwargs.pop("domain_size", None)
    max_pos = kwargs.pop("max_pos", None)
    if kwargs:
        unknown = ", ".join(sorted(kwargs))
        raise TypeError(f"Unexpected keyword argument(s): {unknown}")
    if not args:
        raise TypeError("double_hash requires ID_interval")
    if len(args) > 3:
        raise TypeError("double_hash accepts at most 3 positional arguments")

    ID_interval = args[0]
    if len(args) >= 2:
        if domain_size is not None:
            raise TypeError("domain_size provided both positionally and by keyword")
        domain_size = args[1]
    if len(args) >= 3:
        if max_pos is not None:
            raise TypeError("max_pos provided both positionally and by keyword")
        max_pos = args[2]

    ID_interval = {ID: (int(st), int(ed)) for ID, (st, ed) in ID_interval.items()}
    if max_pos is None:
        max_pos = _infer_max_pos(ID_interval)
    max_pos = int(max_pos)
    if domain_size is None:
        domain_size = _infer_domain_size(max_pos)
    domain_size = int(domain_size)
    if domain_size <= 0:
        raise ValueError("domain_size must be positive")
    return ID_interval, domain_size, max_pos, silent


class RegularGenomeBins:
    """Fast overlap index for regular fixed-size or sliding genomic bins."""

    def __init__(self, *args, **kwargs):
        ID_interval, bin_size, bin_step, max_pos, silent = _normalize_bin_constructor_args(args, kwargs)
        self.ID_value = {}
        self.ID_interval = ID_interval
        self.bin_size = bin_size
        self.bin_step = bin_step
        self.max_pos = max_pos
        self.idx_ID = {}
        self.ID_idx = {}

        for ID, (st, ed) in self.ID_interval.items():
            if st % self.bin_step != 0:
                raise ValueError("RegularGenomeBins requires starts aligned to bin_step")
            if ed != st + self.bin_size:
                raise ValueError("RegularGenomeBins requires fixed-length intervals")
            idx = st // self.bin_step
            if idx in self.idx_ID:
                raise ValueError(f"Duplicate bin index: {idx}")
            self.idx_ID[idx] = ID
            self.ID_idx[ID] = idx

        if not silent:
            print("hash function is built", file=sys.stderr)

    def find_overlapping_position(self, pos):
        pos = int(pos)
        find_IDs = []
        idx = pos // self.bin_step
        min_idx = max(0, (pos - self.bin_size + 1 + self.bin_step - 1) // self.bin_step)
        for idx in range(idx, min_idx - 1, -1):
            ID = self.idx_ID.get(idx)
            if ID is None:
                continue
            st, ed = self.ID_interval[ID]
            if st <= pos < ed:
                find_IDs.append(ID)
        return find_IDs

    def find_overlapping_interval(self, start, end):
        start, end = int(start), int(end)
        if end <= start:
            return []
        min_idx = max(0, (start - self.bin_size + 1 + self.bin_step - 1) // self.bin_step)
        max_idx = (min(end, self.max_pos + 1) - 1) // self.bin_step
        find_IDs = []
        for idx in range(min_idx, max_idx + 1):
            ID = self.idx_ID.get(idx)
            if ID is None:
                continue
            st, ed = self.ID_interval[ID]
            if _overlap_len(st, ed, start, end) > 0:
                find_IDs.append(ID)
        return find_IDs

    def add_position_signal(self, pos, value):
        find_IDs = self.find_overlapping_position(pos)
        for ID in find_IDs:
            self.ID_value.setdefault(ID, 0.0)
            self.ID_value[ID] += value
        return find_IDs

    def add_interval_signal(self, start, end, value):
        find_IDs = self.find_overlapping_interval(start, end)
        for ID in find_IDs:
            st, ed = self.ID_interval[ID]
            length = _overlap_len(st, ed, int(start), int(end))
            self.ID_value.setdefault(ID, 0.0)
            self.ID_value[ID] += value * length
        return find_IDs

    def find(self, pos):
        return self.find_overlapping_position(pos)

    def find_range(self, rst, red):
        return self.find_overlapping_interval(rst, red)

    def insert(self, *args):
        if len(args) == 2:
            return self.add_position_signal(args[0], args[1])
        if len(args) == 3:
            return self.add_interval_signal(args[0], args[1], args[2])
        raise TypeError("insert expects (pos, value) or (start, end, value)")

    def insert_range(self, rst, red, value):
        return self.add_interval_signal(rst, red, value)

    def keys(self):
        return self.ID_value.keys()

    def values(self):
        return self.ID_value.values()

    def __getitem__(self, id):
        return self.ID_value.get(id, 0.0)

    def __setitem__(self, id, value):
        self.ID_value[id] = value

    def ID(self, id):
        return self.ID_value[id]

    def get(self):
        return self.ID_value

    def clear(self):
        self.ID_value = {}


class GenomicIntervalIndex:
    """Fast overlap index for irregular genomic annotations."""

    def __init__(self, *args, **kwargs):
        ID_interval, domain_size, max_pos, silent = _normalize_interval_constructor_args(args, kwargs)
        self.ID_value = {}
        self.ID_interval = ID_interval
        self.domain_size = domain_size
        self.max_pos = max_pos
        self.domain_num = max_pos // domain_size + 1
        self.domain_IDs = {i: [] for i in range(self.domain_num)}
        self._edID = sorted((ed, ID) for ID, (st, ed) in self.ID_interval.items())
        self._edlist = [ed for ed, ID in self._edID]
        self._IDlist = [ID for ed, ID in self._edID]

        for domain in range(self.domain_num):
            dst = domain * self.domain_size
            ded = min(dst + self.domain_size, self.max_pos + 1)
            idx = _binary_search_left(self._edlist, dst + 1)
            for j in range(idx, len(self._edlist)):
                ID = self._IDlist[j]
                st, ed = self.ID_interval[ID]
                if st < ded and ed > dst:
                    self.domain_IDs[domain].append(ID)

        if not silent:
            print("hash function is built", file=sys.stdout)

    def __str__(self):
        rows = ["ID\tst\ted\tvalue"]
        for ID, value in self.ID_value.items():
            st, ed = self.ID_interval[ID]
            rows.append(f"{ID}\t{st}\t{ed}\t{value}")
        return "\n".join(rows)

    def _candidate_ids_for_interval(self, start, end):
        if end <= start:
            return []
        domain1 = max(0, start // self.domain_size)
        domain2 = min(self.domain_num - 1, (end - 1) // self.domain_size)
        IDs = set()
        for domain in range(domain1, domain2 + 1):
            IDs.update(self.domain_IDs.get(domain, []))
        return IDs

    def find_overlapping_position(self, pos):
        pos = int(pos)
        domain = pos // self.domain_size
        IDs = self.domain_IDs.get(domain, [])
        edID = sorted((self.ID_interval[ID][1], ID) for ID in IDs)
        find_IDs = []
        for ed, ID in edID:
            st, _ = self.ID_interval[ID]
            if st <= pos < ed:
                find_IDs.append(ID)
        return find_IDs

    def find_overlapping_interval(self, start, end):
        start, end = int(start), int(end)
        IDs = self._candidate_ids_for_interval(start, end)
        edID = sorted((self.ID_interval[ID][1], ID) for ID in IDs)
        find_IDs = []
        for ed, ID in edID:
            st, ed = self.ID_interval[ID]
            if _overlap_len(st, ed, start, end) > 0:
                find_IDs.append(ID)
        return find_IDs

    def add_position_signal(self, pos, value):
        find_IDs = self.find_overlapping_position(pos)
        for ID in find_IDs:
            self.ID_value.setdefault(ID, 0.0)
            self.ID_value[ID] += value
        return find_IDs

    def add_interval_signal(self, start, end, value):
        find_IDs = self.find_overlapping_interval(start, end)
        for ID in find_IDs:
            st, ed = self.ID_interval[ID]
            length = _overlap_len(st, ed, int(start), int(end))
            self.ID_value.setdefault(ID, 0.0)
            self.ID_value[ID] += value * length
        return find_IDs

    def find(self, pos):
        return self.find_overlapping_position(pos)

    def find_range(self, rst, red):
        return self.find_overlapping_interval(rst, red)

    def insert(self, *args):
        if len(args) == 2:
            return self.add_position_signal(args[0], args[1])
        if len(args) == 3:
            return self.add_interval_signal(args[0], args[1], args[2])
        raise TypeError("insert expects (pos, value) or (start, end, value)")

    def insert_range(self, rst, red, value):
        return self.add_interval_signal(rst, red, value)

    def keys(self):
        return self.ID_value.keys()

    def values(self):
        return self.ID_value.values()

    def __getitem__(self, id):
        return self.ID_value.get(id, 0.0)

    def __setitem__(self, id, value):
        self.ID_value[id] = value

    def ID(self, id):
        return self.ID_value[id]

    def get(self):
        return self.ID_value

    def clear(self):
        self.ID_value = {}


bin_hash = RegularGenomeBins
double_hash = GenomicIntervalIndex

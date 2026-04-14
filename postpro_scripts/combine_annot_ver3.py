#!/usr/bin/env python3
"""
combine_annot_ver3.py

Unified annotation combiner.

This script merges the capabilities of:
- postpro_scripts/combine_annot.py
- postpro_scripts/combine_annot_edit.py
- postpro_scripts/combine_annot_ver2.py
- postpro_scripts/combine_annot_ver2.1.py

into a single, Python3-compatible entrypoint.

Two mutually exclusive input modes:
- Ncov mode:   --Ncov + --Nlen (optionally --Bin)  -> compute condensability-like metric vs control
- score mode:  --score + --binsize/--binstep       -> treat input as a generic score file

Output:
- point inputs -> <out>_anot.cn (SNP/PhysicalPosition schema)
- binned inputs -> <out>_anot.txt (BinID/Start/End schema)
"""

from __future__ import annotations

import argparse
import copy
import math
import sys
from bisect import bisect_left
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple


def rev_cmp(seq: str) -> str:
    dic = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    out = []
    for nt in seq.upper():
        out.append(dic.get(nt, "N"))
    return "".join(out)[::-1]


def AT_content(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return float("nan")
    at = 0
    for nt in seq:
        if nt in "AT":
            at += 1
    return float(at) / float(len(seq))


def C_motif(seq: str, motif: str = "CG", both: bool = True) -> int:
    seq = seq.upper()
    motif = motif.upper()
    n = 0
    for i in range(0, max(0, len(seq) - len(motif) + 1)):
        if seq[i : i + len(motif)] == motif:
            n += 1
    if both:
        r = rev_cmp(seq)
        for i in range(0, max(0, len(r) - len(motif) + 1)):
            if r[i : i + len(motif)] == motif:
                n += 1
    return n


def _parse_name_file_pairs(raw: Optional[List[str]], flag: str) -> Dict[str, str]:
    if not raw:
        return {}
    if len(raw) % 2 != 0:
        raise SystemExit(f"{flag}: expected pairs: <name1> <file1> <name2> <file2> ...")
    out: Dict[str, str] = {}
    it = iter(raw)
    for name, fname in zip(it, it):
        out[name] = fname
    return out


def read_genome_size_from_fasta(ref_fname: str) -> Dict[str, int]:
    genome_size: Dict[str, int] = {}
    key: Optional[str] = None
    with open(ref_fname) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                key = line[1:].split()[0]
                if key in genome_size:
                    raise ValueError(f"duplicate contig in fasta: {key}")
                genome_size[key] = 0
                continue
            if key is None:
                raise ValueError("FASTA appears to be missing a header line starting with '>'")
            genome_size[key] += len(line)
    return genome_size


def read_score_file(
    score_fname: str,
    chr_list: Iterable[str],
    bin_size: int,
    col_st: int = 3,
    offset: Optional[int] = None,
) -> Tuple[Dict[int, Tuple[str, int]], Dict[str, Dict[int, int]], List[Dict[int, float]], List[str], str]:
    """
    Read a score file.

    Supported formats:
    - point cn-like:  ID  chr  PhysicalPosition  score1 score2 ...
    - binned:         ID  chr  Start End         score1 score2 ...
      (Start is used as the bin start)
    """
    chr_set = set(chr_list)
    ID_chrst: Dict[int, Tuple[str, int]] = {}
    chrst_ID: Dict[str, Dict[int, int]] = {}
    ID_score_list: List[Dict[int, float]] = []
    names: List[str] = []
    data_type: Optional[str] = None

    first = True
    with open(score_fname) as f:
        for line in f:
            cols = line.strip().split()
            if not cols:
                continue
            if first:
                if len(cols) < 4:
                    raise ValueError("score file header too short")
                # cn file: point data
                if cols[2] == "PhysicalPosition":
                    data_type = "point"
                    col_st = 3
                    if offset is None:
                        offset = -(bin_size // 2)  # point is middle of bin
                else:
                    # otherwise, assume binned data
                    if cols[2] != "Start" or cols[3] != "End":
                        raise ValueError("unrecognized score header (expected PhysicalPosition or Start/End)")
                    data_type = "binned"
                    col_st = 4
                    offset = 0
                names = [c.rsplit("/", 1)[-1].rsplit(".", 1)[0] for c in cols[col_st:]]
                ID_score_list = [{} for _ in range(len(names))]
                first = False
                continue

            if data_type is None or offset is None:
                raise AssertionError("internal state error while reading score file")

            ID = int(cols[0])
            chrname = cols[1]
            if chr_set and chrname not in chr_set:
                continue

            # For point: cols[2] is PhysicalPosition
            # For binned: cols[2] is Start
            pos = int(cols[2])
            scores = [float(s) for s in cols[col_st:]]
            st = pos + int(offset)

            chrst_ID.setdefault(chrname, {})[st] = ID
            ID_chrst[ID] = (chrname, st)
            for i, s in enumerate(scores):
                ID_score_list[i][ID] = s

    if data_type is None:
        raise ValueError("empty score file?")
    return ID_chrst, chrst_ID, ID_score_list, names, data_type


def read_Ncov_file_as_scores(
    Ncov_fname: str,
    chr_list: Iterable[str],
    NCP_len: int,
    metric_mode: str,
) -> Tuple[Dict[int, Tuple[str, int]], Dict[str, Dict[int, int]], List[Dict[int, float]], List[str], str]:
    """
    Read an Ncov-style file:
      header: ID chr dyad <track1> <track2> ... <control>
      rows:   ID chr dyad count1 count2 ... control_count

    Converts each track to a score using the legacy metric (track vs control).
    Returns point-style data (PhysicalPosition implied as st + NCP_len/2).
    """
    chr_set = set(chr_list)

    def metric_v1(test: float, control: float) -> float:
        if test <= 0:
            test = 1.0
        if control <= 0:
            control = 1.0
        rcount = float(test) / float(control)
        return -math.log(rcount)

    def metric_v1_edit(test: float, control: float) -> float:
        if control <= 0:
            return float("nan")
        test = test + 1.0
        control = control + 1.0
        rcount = float(test) / float(control)
        return -math.log(rcount)

    if metric_mode == "v1":
        metric_fn = metric_v1
    elif metric_mode == "v1_edit":
        metric_fn = metric_v1_edit
    else:
        raise ValueError(f"unknown metric mode: {metric_mode}")

    ID_chrst: Dict[int, Tuple[str, int]] = {}
    chrst_ID: Dict[str, Dict[int, int]] = {}
    ID_score_list: List[Dict[int, float]] = []
    names: List[str] = []

    first = True
    with open(Ncov_fname) as f:
        for line in f:
            cols = line.strip().split()
            if not cols:
                continue
            if first:
                names_all = cols[3:]
                if len(names_all) < 2:
                    raise ValueError("Ncov header must include at least one track and one control column")
                names = names_all[:-1]  # exclude control from outputs
                ID_score_list = [{} for _ in range(len(names))]
                first = False
                continue

            ID = int(cols[0])
            chrname = cols[1]
            if chr_set and chrname not in chr_set:
                continue
            dyad = int(cols[2])
            counts = [float(c) for c in cols[3:]]
            if sum(counts) <= 0:
                continue
            control = counts[-1]
            st = dyad - (NCP_len // 2)

            chrst_ID.setdefault(chrname, {})[st] = ID
            ID_chrst[ID] = (chrname, st)
            for i in range(len(names)):
                ID_score_list[i][ID] = metric_fn(counts[i], control)

    return ID_chrst, chrst_ID, ID_score_list, names, "point"


def get_seq(ref_fname: str, chrname: str, st_ID: Dict[int, int], win: int) -> Dict[int, str]:
    """
    Stream FASTA and extract sequences for each (start -> ID) in st_ID.
    This is memory-light (does not load full chromosome sequence into RAM).
    """
    if not st_ID:
        return {}

    stID = [(st, st_ID[st]) for st in sorted(st_ID.keys())]
    k = 0
    pos, ID = stID[k]

    ID_seq: Dict[int, str] = {}
    left: List[Tuple[int, str]] = []
    pt = -1  # 0-based coordinate of end-of-processed sequence (exclusive) minus 1 (legacy style)
    in_chr = False

    with open(ref_fname) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if in_chr:
                    break
                in_chr = (line[1:].split()[0] == chrname)
                continue
            if not in_chr:
                continue

            # Fast-forward to near first requested position.
            if not left and pt + len(line) < pos:
                pt += len(line)
                continue

            # Fill any previously-started sequences that spilled across lines.
            if left:
                new_left: List[Tuple[int, str]] = []
                for leftID, seq in left:
                    need = win - len(seq)
                    take = min(len(line), need)
                    seq2 = seq + line[:take]
                    if len(seq2) >= win:
                        ID_seq[leftID] = seq2[:win]
                    else:
                        new_left.append((leftID, seq2))
                left = new_left

            # Start sequences at any positions that fall within this line.
            while pt + len(line) >= pos and k < len(stID):
                loc = pos - pt - 1  # convert to local index within line
                frag = line[loc : min(loc + win, len(line))]
                if len(frag) >= win:
                    ID_seq[ID] = frag[:win]
                else:
                    left.append((ID, frag))
                k += 1
                if k >= len(stID):
                    break
                pos, ID = stID[k]

            if not left and len(ID_seq) == len(stID):
                break
            pt += len(line)

    # Whatever remains at end-of-contig becomes shorter sequence (keeps old behavior).
    for leftID, seq in left:
        ID_seq[leftID] = seq

    if len(ID_seq) != len(stID):
        missing = len(stID) - len(ID_seq)
        raise RuntimeError(f"failed to extract {missing} sequences for chr={chrname}")
    return ID_seq


class bin_hash:
    """Fast interval hash when bin_size/bin_step are constant."""

    def __init__(self, ID_interval: Dict[int, Tuple[int, int]], bin_size: int, bin_step: int, max_pos: int):
        self.ID_value: Dict[int, float] = {}
        self.ID_interval = ID_interval
        self.bin_size = int(bin_size)
        self.bin_step = int(bin_step)
        self.max_pos = int(max_pos)

        self.idx_ID: Dict[int, int] = {}
        self.ID_idx: Dict[int, int] = {}
        for ID, (st, ed) in ID_interval.items():
            if st % self.bin_step != 0:
                raise ValueError("bin_hash requires starts aligned to bin_step")
            if ed != st + self.bin_size:
                raise ValueError("bin_hash requires fixed-length bins")
            idx = st // self.bin_step
            if idx in self.idx_ID:
                raise ValueError("duplicate bin index")
            self.idx_ID[idx] = ID
            self.ID_idx[ID] = idx

    def find(self, pos: int) -> List[int]:
        find_IDs: List[int] = []
        idx = int(pos) // self.bin_step
        st = self.bin_step * idx
        ed = st + self.bin_size
        while st <= pos < ed:
            ID = self.idx_ID.get(idx)
            if ID is not None:
                find_IDs.append(ID)
            idx -= 1
            if idx < 0:
                break
            st = self.bin_step * idx
            ed = st + self.bin_size
        return find_IDs

    def insert(self, pos: int, value: float) -> List[int]:
        find_IDs = self.find(pos)
        for ID in find_IDs:
            self.ID_value[ID] = self.ID_value.get(ID, 0.0) + float(value)
        return find_IDs

    def find_range(self, rst: int, red: int) -> List[int]:
        find_IDs: List[int] = []
        if red <= rst:
            return find_IDs

        idx = int(rst) // self.bin_step
        min_idx = idx
        st = self.bin_step * idx
        ed = st + self.bin_size
        while st <= rst < ed:
            min_idx = min(min_idx, idx)
            idx -= 1
            if idx < 0:
                break
            st = self.bin_step * idx
            ed = st + self.bin_size

        red = min(int(red), self.max_pos + 1)
        max_idx = (red - 1) // self.bin_step

        for idx in range(min_idx, max_idx + 1):
            ID = self.idx_ID.get(idx)
            if ID is not None:
                find_IDs.append(ID)
        return find_IDs

    def insert_range(self, rst: int, red: int, value: float) -> List[int]:
        find_IDs = self.find_range(rst, red)
        for ID in find_IDs:
            idx = self.ID_idx[ID]
            st = self.bin_step * idx
            ed = st + self.bin_size
            a, b = max(st, rst), min(ed, red)
            length = max(0, b - a)
            self.ID_value[ID] = self.ID_value.get(ID, 0.0) + float(value) * float(length)
        return find_IDs

    def get(self) -> Dict[int, float]:
        return self.ID_value


class double_hash:
    """
    Interval hash for arbitrary intervals.
    Uses coarse domains (chunks) to reduce the search space.
    """

    def __init__(self, ID_interval: Dict[int, Tuple[int, int]], domain_size: int, max_pos: int):
        self.ID_value: Dict[int, float] = {}
        self.ID_interval = ID_interval
        self.domain_size = int(domain_size)
        self.max_pos = int(max_pos)

        # Precompute per-domain interval lists sorted by end position.
        self.domain_num = self.max_pos // self.domain_size + 1
        self.domain_edlist: List[List[int]] = [[] for _ in range(self.domain_num)]
        self.domain_IDlist: List[List[int]] = [[] for _ in range(self.domain_num)]

        # Global end-sorted list for initial filtering.
        edID = sorted([(ed, ID) for ID, (st, ed) in ID_interval.items()], key=lambda x: x[0])
        edlist = [ed for ed, _ in edID]
        IDlist = [ID for _, ID in edID]

        for i in range(self.domain_num):
            dst = i * self.domain_size
            ded = min(dst + self.domain_size, self.max_pos + 1)
            idx1 = bisect_left(edlist, dst)
            if idx1 >= len(edlist):
                continue
            pairs: List[Tuple[int, int]] = []
            for j in range(idx1, len(edlist)):
                ID = IDlist[j]
                st, ed = ID_interval[ID]
                if st < ded:
                    pairs.append((ed, ID))
            pairs.sort(key=lambda x: x[0])
            self.domain_edlist[i] = [ed for ed, _ in pairs]
            self.domain_IDlist[i] = [ID for _, ID in pairs]

    def find(self, pos: int) -> List[int]:
        find_IDs: List[int] = []
        domain = int(pos) // self.domain_size
        if domain < 0 or domain >= self.domain_num:
            return find_IDs
        edlist = self.domain_edlist[domain]
        IDlist = self.domain_IDlist[domain]
        idx = bisect_left(edlist, pos)
        for k in range(idx, len(IDlist)):
            ID = IDlist[k]
            st, ed = self.ID_interval[ID]
            if st <= pos < ed:
                find_IDs.append(ID)
        return find_IDs

    def insert(self, pos: int, value: float) -> List[int]:
        find_IDs = self.find(pos)
        for ID in find_IDs:
            self.ID_value[ID] = self.ID_value.get(ID, 0.0) + float(value)
        return find_IDs

    def find_range(self, rst: int, red: int) -> List[int]:
        if red <= rst:
            return []
        domain1 = int(rst) // self.domain_size
        domain2 = int(red) // self.domain_size
        domain1 = max(0, domain1)
        domain2 = min(self.domain_num - 1, domain2)

        out: set[int] = set()
        for domain in range(domain1, domain2 + 1):
            for ID in self.domain_IDlist[domain]:
                st, ed = self.ID_interval[ID]
                if st < red and ed > rst:
                    out.add(ID)
        return list(out)

    def insert_range(self, rst: int, red: int, value: float) -> List[int]:
        find_IDs = self.find_range(rst, red)
        for ID in find_IDs:
            st, ed = self.ID_interval[ID]
            a, b = max(st, rst), min(ed, red)
            length = max(0, b - a)
            self.ID_value[ID] = self.ID_value.get(ID, 0.0) + float(value) * float(length)
        return find_IDs

    def get(self) -> Dict[int, float]:
        return self.ID_value


def read_BS_file(fname: str, Int_dict, chrname: str) -> Tuple[Dict[int, float], Dict[int, float]]:
    ID_meC: Dict[int, float] = {}
    with open(fname) as f:
        for line in f:
            cols = line.strip().split()
            if len(cols) < 11:
                continue
            chr0, st, _ed, _name, _score, _strand, _a, _b, _c, reads, frac = cols[:11]
            if chr0 != chrname:
                continue
            reads_i = int(reads)
            if reads_i <= 0:
                continue
            pos = int(st)
            frac_f = 0.01 * float(frac)
            findIDs = Int_dict.insert(pos, 1.0)
            for ID in findIDs:
                ID_meC[ID] = ID_meC.get(ID, 0.0) + frac_f
    ID_C = Int_dict.get()
    return ID_C, ID_meC


def read_chip_file(fname: str, Int_dict, chrname: str, unit: str = "signal") -> Dict[int, float]:
    with open(fname) as f:
        for line in f:
            cols = line.strip().split()
            if len(cols) < 9:
                continue
            chr0, st, ed, _peakname, _score, _strand, signal, pvalue, qvalue = cols[:9]
            if chr0 != chrname:
                continue
            if unit == "signal":
                score = float(signal)
            elif unit == "pvalue":
                score = float(pvalue)
            elif unit == "qvalue":
                score = float(qvalue)
            else:
                raise ValueError(f"unknown chip unit: {unit}")
            Int_dict.insert_range(int(st), int(ed), score)
    return Int_dict.get()


def read_bedgraph_file(fname: str, Int_dict, chrname: str) -> Dict[int, float]:
    with open(fname) as f:
        for line in f:
            if not line.startswith("chr"):
                continue
            cols = line.strip().split()
            if len(cols) < 4:
                continue
            chr0, st, ed, value = cols[:4]
            if chr0 != chrname:
                continue
            Int_dict.insert_range(int(st), int(ed), float(value))
    return Int_dict.get()


def combine_all(
    *,
    input_mode: str,
    score_fname: Optional[str],
    Ncov_fname: Optional[str],
    metric_mode: str,
    ref_fname: str,
    bin_size: int,
    bin_step: int,
    bs_fname: Dict[str, str],
    chip_fname: Dict[str, str],
    bedgraph_fname: Dict[str, str],
    full_seq: bool,
    chr_list: List[str],
    genome_size: Dict[str, int],
    out_fname: str,
) -> None:
    # Read base scores
    if input_mode == "score":
        ID_chrst, chrst_ID, ID_score_list, names, data_type = read_score_file(score_fname, chr_list, bin_size)
    elif input_mode == "Ncov":
        ID_chrst, chrst_ID, ID_score_list, names, data_type = read_Ncov_file_as_scores(
            Ncov_fname, chr_list, bin_size, metric_mode
        )
    else:
        raise ValueError(f"unknown input_mode: {input_mode}")

    chr_list = sorted(chrst_ID.keys())
    print("base score reading is done", file=sys.stderr)

    # Extract sequences for each bin (needed at least for ATcontent).
    ID_seq: Dict[int, str] = {}
    for chrname in chr_list:
        ID_seq.update(get_seq(ref_fname, chrname, chrst_ID[chrname], bin_size))
    print("reference reading is done", file=sys.stderr)

    # Build per-chromosome interval dictionary templates for annotations.
    chr_ID_interval: Dict[str, Dict[int, Tuple[int, int]]] = {}
    for chrname in chr_list:
        st_ID = chrst_ID[chrname]
        ID_interval: Dict[int, Tuple[int, int]] = {}
        for st, ID in st_ID.items():
            ID_interval[ID] = (int(st), int(st) + int(bin_size))
        chr_ID_interval[chrname] = ID_interval

    def new_intdict(chrname: str):
        if bin_step and bin_step > 0:
            return bin_hash(chr_ID_interval[chrname], bin_size, bin_step, genome_size[chrname])
        return double_hash(chr_ID_interval[chrname], 100000, genome_size[chrname])

    # Read BS
    bs_ID_C: Dict[str, Dict[int, float]] = {}
    bs_ID_meC: Dict[str, Dict[int, float]] = {}
    if bs_fname:
        for chrname in chr_list:
            for bs, fname in bs_fname.items():
                intdict = new_intdict(chrname)
                ID_C, ID_meC = read_BS_file(fname, intdict, chrname)
                bs_ID_C.setdefault(bs, {}).update(ID_C)
                bs_ID_meC.setdefault(bs, {}).update(ID_meC)
        print("BS reading is done", file=sys.stderr)

    # Read chip
    chip_ID_signal: Dict[str, Dict[int, float]] = {}
    if chip_fname:
        for chrname in chr_list:
            for chip, fname in chip_fname.items():
                intdict = new_intdict(chrname)
                ID_signal = read_chip_file(fname, intdict, chrname)
                chip_ID_signal.setdefault(chip, {}).update(ID_signal)
        print("Chip reading is done", file=sys.stderr)

    # Read bedgraph
    bg_ID_value: Dict[str, Dict[int, float]] = {}
    if bedgraph_fname:
        for chrname in chr_list:
            for bg, fname in bedgraph_fname.items():
                intdict = new_intdict(chrname)
                ID_value = read_bedgraph_file(fname, intdict, chrname)
                bg_ID_value.setdefault(bg, {}).update(ID_value)
        print("bedgraph file reading is done", file=sys.stderr)

    print("writing annotation file", file=sys.stderr)

    if data_type == "point":
        out_path = out_fname + "_anot.cn"
        header = "SNP\tChromosome\tPhysicalPosition"
    else:
        out_path = out_fname + "_anot.txt"
        header = "BinID\tChromosome\tStart\tEnd"

    with open(out_path, "w") as out:
        s = header
        for name in names:
            s += "\t" + name
        if full_seq:
            s += "\tSequence"
        s += "\tATcontent"
        for bs in sorted(bs_ID_C.keys()):
            s += f"\tCNumber({bs})\tmeCNumber({bs})"
        for chip in sorted(chip_ID_signal.keys()):
            s += "\t" + chip
        for bg in sorted(bg_ID_value.keys()):
            s += "\t" + bg
        out.write(s + "\n")

        IDs = sorted(ID_score_list[0].keys()) if ID_score_list else sorted(ID_chrst.keys())
        for ID in IDs:
            chrname, st = ID_chrst[ID]
            if data_type == "point":
                pos = int(st) + (bin_size // 2)
                row = f"{ID}\t{chrname}\t{pos}"
            else:
                ed = int(st) + int(bin_size)
                row = f"{ID}\t{chrname}\t{int(st)}\t{ed}"

            for ID_score in ID_score_list:
                val = ID_score.get(ID)
                if val is None:
                    row += "\tNA"
                elif isinstance(val, float):
                    if math.isnan(val):
                        row += "\tNA"
                    else:
                        row += "\t" + str(round(val, 5))
                else:
                    row += "\t" + str(val)

            seq = ID_seq[ID]
            if full_seq:
                row += "\t" + seq
            row += "\t" + str(round(AT_content(seq), 5))

            for bs in sorted(bs_ID_C.keys()):
                row += "\t" + str(int(bs_ID_C.get(bs, {}).get(ID, 0.0)))
                row += "\t" + str(round(bs_ID_meC.get(bs, {}).get(ID, 0.0), 5))

            for chip in sorted(chip_ID_signal.keys()):
                row += "\t" + str(round(chip_ID_signal.get(chip, {}).get(ID, 0.0), 5))

            for bg in sorted(bg_ID_value.keys()):
                row += "\t" + str(round(bg_ID_value.get(bg, {}).get(ID, 0.0), 5))

            out.write(row + "\n")

    print(f"Done: {out_path}", file=sys.stderr)


def _str2bool(v: str) -> bool:
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    if v.lower() in ("no", "false", "f", "n", "0"):
        return False
    raise argparse.ArgumentTypeError("Boolean value expected.")


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Combine epigenetic annotations with score/Ncov inputs")

    input_grp = parser.add_mutually_exclusive_group(required=True)
    input_grp.add_argument("--score", dest="score_fname", type=str, help="score file (point or binned)")
    input_grp.add_argument("--Ncov", dest="Ncov_fname", type=str, help="Ncov-style file (dyad + tracks + control)")

    parser.add_argument("--ref", dest="ref_fname", type=str, required=True, help="reference genome fasta")
    parser.add_argument("--chr", dest="chr_list", type=str, nargs="+", help="target chromosome list (default: all)")

    # score-mode binning
    parser.add_argument("--binsize", dest="bin_size", type=int, default=None, help="bin size (score mode) / Nlen (Ncov mode)")
    parser.add_argument("--binstep", dest="bin_step", type=int, default=0, help="bin step (0 => use interval hash)")

    # Ncov-mode params
    parser.add_argument("--Nlen", dest="NCP_len", type=int, default=None, help="NCP length (required for --Ncov unless --Bin is used)")
    parser.add_argument("--Bin", dest="Bin", type=int, nargs=2, help="legacy: provide bin size and step (Ncov mode)")
    parser.add_argument("--metric", dest="metric_mode", choices=["v1", "v1_edit"], default="v1_edit", help="metric behavior for --Ncov")

    # optional annotations
    parser.add_argument("--bs", dest="bs_fname", type=str, nargs="+", help="pairs: <name> <bedMethyl file> ...")
    parser.add_argument("--chip", dest="chip_fname", type=str, nargs="+", help="pairs: <name> <narrowPeak file> ...")
    parser.add_argument("--bedgraph", dest="bedgraph_fname", type=str, nargs="+", help="pairs: <name> <bedgraph file> ...")

    parser.add_argument("--full-seq", dest="full_seq", type=_str2bool, nargs="?", const=True, default=False, help="write full sequence column")
    parser.add_argument("-o", dest="out_fname", type=str, default="output", help="output prefix")

    args = parser.parse_args(argv)

    genome_size = read_genome_size_from_fasta(args.ref_fname)
    chr_list = sorted(args.chr_list) if args.chr_list else sorted(genome_size.keys())

    bs_fname = _parse_name_file_pairs(args.bs_fname, "--bs")
    chip_fname = _parse_name_file_pairs(args.chip_fname, "--chip")
    bedgraph_fname = _parse_name_file_pairs(args.bedgraph_fname, "--bedgraph")

    if args.score_fname:
        input_mode = "score"
        if args.bin_size is None:
            raise SystemExit("--binsize is required with --score")
        bin_size = int(args.bin_size)
        bin_step = int(args.bin_step or 0)
        combine_all(
            input_mode=input_mode,
            score_fname=args.score_fname,
            Ncov_fname=None,
            metric_mode=args.metric_mode,
            ref_fname=args.ref_fname,
            bin_size=bin_size,
            bin_step=bin_step,
            bs_fname=bs_fname,
            chip_fname=chip_fname,
            bedgraph_fname=bedgraph_fname,
            full_seq=bool(args.full_seq),
            chr_list=chr_list,
            genome_size=genome_size,
            out_fname=args.out_fname,
        )
        return 0

    input_mode = "Ncov"
    if args.Bin:
        bin_size, bin_step = int(args.Bin[0]), int(args.Bin[1])
    else:
        if args.NCP_len is None:
            raise SystemExit("--Nlen is required with --Ncov (unless --Bin is provided)")
        bin_size = int(args.NCP_len)
        bin_step = 0

    combine_all(
        input_mode=input_mode,
        score_fname=None,
        Ncov_fname=args.Ncov_fname,
        metric_mode=args.metric_mode,
        ref_fname=args.ref_fname,
        bin_size=bin_size,
        bin_step=bin_step,
        bs_fname=bs_fname,
        chip_fname=chip_fname,
        bedgraph_fname=bedgraph_fname,
        full_seq=bool(args.full_seq),
        chr_list=chr_list,
        genome_size=genome_size,
        out_fname=args.out_fname,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

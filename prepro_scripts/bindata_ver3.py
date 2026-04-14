import sys, math, gzip
from argparse import ArgumentParser
import Helper_Py3_compat as Helper_Py3


def _validate_inputs(bin_size, bin_step, bin_value):
    if bin_size <= 0:
        raise ValueError("bin_size must be positive")
    if bin_step <= 0:
        raise ValueError("bin_step must be positive")
    if bin_value not in {"sum", "mean"}:
        raise ValueError("bin_value must be either 'sum' or 'mean'")


def _log_progress(chrom, start, order):
    current_order = int(math.log10(max(start, 1)))
    if order is None:
        print(chrom + " st " + str(start) + " is reading", file=sys.stderr)
        return current_order
    if current_order > order:
        print(chrom + " st " + str(start) + " is reading", file=sys.stderr)
        return order + 1
    return order


def _accumulate_bin(chr_Bdata, chr_Bcount, chrom, idx, labels, values, overlap_len, bin_value):
    name_data = chr_Bdata.setdefault(chrom, {}).setdefault(idx, {})
    name_count = None
    if bin_value == "mean":
        name_count = chr_Bcount.setdefault(chrom, {}).setdefault(idx, {})

    for name, value in zip(labels, values):
        name_data[name] = name_data.get(name, 0.0) + value * overlap_len
        if name_count is not None:
            name_count[name] = name_count.get(name, 0) + 1


def _build_output_row(chrom, idx, chr_size, name_Bdata, chr_Bcount, labels, bin_size, bin_step, bin_value):
    values = []
    for name in labels:
        bin_data = name_Bdata.get(name, 0.0)
        if bin_value == "mean":
            bin_count = chr_Bcount.get(name, 0)
            if bin_count <= 0:
                value = "NA"
            else:
                value = float(bin_data) / bin_count
        else:
            value = bin_data
        values.append(value)

    bin_start = bin_step * idx
    bin_end = min(bin_start + bin_size, chr_size[chrom])
    row = [chrom, str(bin_start), str(bin_end)]
    row += [str(value) for value in values]
    return "\t".join(row)


def _write_output(out_fname, chr_Bdata, chr_Bcount, chr_size, chr_list, labels, bin_size, bin_step, bin_value, skip_zero):
    print("writing bin data file", file=sys.stderr)

    with gzip.open(out_fname + "_Bdata.gtab.gz", "wt", encoding="utf-8", newline="\n") as f:
        header = "Chromosome\tStart\tEnd\t" + "\t".join(labels)
        print(header, file=f)

        for chrom in chr_list:
            chrom_bins = chr_Bdata.get(chrom)
            if not chrom_bins:
                continue

            last_idx = chr_size[chrom] // bin_step
            for idx in range(last_idx + 1):
                name_Bdata = chrom_bins.get(idx)
                if name_Bdata is None:
                    continue

                if skip_zero and sum(name_Bdata.values()) <= 0:
                    continue

                name_Bcount = chr_Bcount.get(chrom, {}).get(idx, {})
                row = _build_output_row(
                    chrom,
                    idx,
                    chr_size,
                    name_Bdata,
                    name_Bcount,
                    labels,
                    bin_size,
                    bin_step,
                    bin_value,
                )
                print(row, file=f)

    print("Done", file=sys.stderr)


def Bin_data(data_fname,
             chr_size,
             bin_size,
             bin_step,
             bin_value,
             skip_zero,
             chr_list,
             out_fname):
    _validate_inputs(bin_size, bin_step, bin_value)

    print("reading data file", file=sys.stderr)

    chr_Bdata = {}
    chr_Bcount = {}
    labels = None
    order = None
    chr_set = set(chr_list)
    layout = Helper_Py3.BinLayout(bin_size, bin_step)

    for row_labels, _col_st, chrom, start, end, values in Helper_Py3.iter_gtab_records(data_fname):
        if labels is None:
            labels = row_labels

        if chrom not in chr_set:
            continue

        order = _log_progress(chrom, start, order)

        end = min(end, chr_size[chrom])
        if end <= start:
            continue

        for idx, overlap_len in layout.iter_overlaps(start, end):
            _accumulate_bin(chr_Bdata, chr_Bcount, chrom, idx, labels, values, overlap_len, bin_value)

    if labels is None:
        raise ValueError("input file has no GTAB records")

    _write_output(out_fname, chr_Bdata, chr_Bcount, chr_size, chr_list, labels, bin_size, bin_step, bin_value, skip_zero)


if __name__ == '__main__':
    parser = ArgumentParser(description='Binning data gtab file')
    parser.add_argument(metavar = '-f',
                        dest='data_fname',
                        type=str,
                        help='data gt file')
    parser.add_argument(metavar='-x',
                        dest='ref_fname',
                        type=str,
                        help='reference sequence filename')
    parser.add_argument('--Bsize',
                        dest="bin_size",
                        type=int,
                        default=167,
                        help='bin window in bp')
    parser.add_argument('--Bstep',
                        dest="bin_step",
                        type=int,
                        default=25,
                        help='bin moving step in bp')
    parser.add_argument('--Bvalue',
                        dest="bin_value",
                        type=str,
                        default='sum',
                        help='binning value choice (sum/mean)')    
    parser.add_argument('--skip',
                        dest="skip_zero",
                        type=Helper_Py3.str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='skip the zero coverage positions')
    parser.add_argument('--chr',
                        dest="chr_list",
                        type=str,
                        nargs='+',
                        help='tagert chromosome list')
    parser.add_argument('-o',
                        dest='out_fname',
                        default='output',
                        type=str,
                        help='output prefix filename')
    
    args = parser.parse_args()

    # get length for each chromosome
    chr_size = Helper_Py3.genome_sizes(args.ref_fname)

    if not args.chr_list:
        chr_list = sorted(chr_size.keys(), key=Helper_Py3.chr_key)
    else:
        chr_list = sorted(args.chr_list, key=Helper_Py3.chr_key)

    Bin_data (args.data_fname,
              chr_size,
              args.bin_size,
              args.bin_step,
              args.bin_value,
              args.skip_zero,
              chr_list,
              args.out_fname)

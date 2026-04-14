import sys
import argparse
import copy
import gzip
import Helper_Py3_compat as Helper_Py3


def _collect_input_labels(fnames):
    data, label = [], []
    for fname in fnames:
        data.append(fname)
        label.append(fname.rsplit('/', 1)[-1].split('.')[0])
    return data, label


def bin_count (fnames,
               genome_size,
               win_size,
               skip_zero,
               bin_GC,
               tlen_option,
               mm_cutoff,
               min_len,
               max_len,
               chr_list,
               out_fname):

    # gather whole file names
    data, label = _collect_input_labels(fnames)
    chr_choices = set(chr_list)

    # partition the genome
    chr_bins = {}
    for chr in chr_list:
        chr_bins[chr] = [0] * int((genome_size[chr] + win_size - 1) // win_size)
        
    # output genome count dictionary
    out = [copy.deepcopy(chr_bins) for i in range(len(data))]

    # mean tlen per each bin for control reads
    if tlen_option:
        g_tlen = copy.deepcopy(chr_bins)
    
    # open the sam file
    for k in range(len(data)):
        g_count = out[k]      
        filename = data[k]
        print("reading %s" % (filename), file=sys.stderr)
        for cols in Helper_Py3.iter_samtools_view(filename, extra_args=["-F", "0x10"]):
            record = Helper_Py3.parse_sam_record(cols)
            if not Helper_Py3.passes_common_paired_filters(
                record,
                chr_choices=chr_choices,
                min_len=min_len,
                max_len=max_len,
                mm_cutoff=mm_cutoff,
                require_proper_pair=True,
            ):
                continue

            # Keep the historical single-position policy for even-length fragments.
            NCPscore = Helper_Py3.record_center_positions(record, even_policy="left")

            # collect valid data
            for NCPpos, score in NCPscore:
                index = int(NCPpos) // int(win_size)
                g_count[record.rname][index] += score

                # record tlen for control
                if tlen_option and k == len(data)-1:
                    g_tlen[record.rname][index] += abs(record.tlen)

    # summarize the output
    print("writing bin file", file=sys.stderr)

    f = gzip.open(out_fname + '_bin.gtab.gz', 'wt', encoding='utf-8', newline='\n')
    s = 'Chromosome\tStart\tEnd'
    for i in range(len(label)):
        s += '\t' + label[i]
    if bin_GC is not None:
        s += '\t' + 'GCcontent'
    if tlen_option:
        s += '\t' + 'Meantlen'
    print(s, file=f)
        

    for chr in sorted(g_count.keys(), key=Helper_Py3.chr_key):
        for i in range(len(g_count[chr])):
            st, ed = i*win_size, (i+1)*win_size 
            s = chr + "\t" + str(st) + '\t' + str(ed)
            total = 0
            for j in range(len(out)):
                total += out[j][chr][i]
                s += '\t%s' % (out[j][chr][i])
            if skip_zero and total <= 0:
                continue
            if bin_GC is not None:
                BinID = chr + '_' + str(i)
                s += '\t%f' % (bin_GC[BinID])
            if tlen_option:
                if out[-1][chr][i] > 0:
                    mean_tlen = float(g_tlen[chr][i])/out[-1][chr][i]
                else:
                    mean_tlen = 0
                s += '\t%f' % (mean_tlen)
            print(s, file=f)

    f.close()

    print("Done", file=sys.stderr)
    

if __name__ == '__main__':
    str2bool = Helper_Py3.str2bool

    parser = argparse.ArgumentParser(description='Divide the genome into bins and counts the read number',
                                     epilog="Notes from Xin: make sure to activate the 'condense-seq' environment before running.")
    parser.add_argument('-f',
                        dest="fnames",
                        type=str,
                        nargs='+',
                        help='SAM/Bam filenames')
    parser.add_argument('-x',
                        dest='ref_fname',
                        type=str,
                        help='reference sequence filename')
    parser.add_argument('-w',
                        dest="win_size",
                        type=int,
                        default=1001,
                        help='bin window size in bp')
    parser.add_argument('-m',
                        dest="mm_cutoff",
                        type=int,
                        default=10,
                        help='mismatch cut-off in bp')
    parser.add_argument('--min',
                        dest="min_len",
                        type=int,
                        default=0,
                        help='minimum length for selection in bp')
    parser.add_argument('--max',
                        dest="max_len",
                        type=int,
                        default=sys.maxsize, # sys.maxint is python 2 annotation, not python 3 annotation.
                        help='maximum length for selection in bp')
    parser.add_argument('--chr',
                        dest="chr_list",
                        type=str,
                        nargs='+',
                        help='target chromosome list')
    parser.add_argument('--skip',
                        dest="skip_zero",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='skip the zero count bins')
    parser.add_argument('--gc',
                        dest="GC_option",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='GC content option')
    parser.add_argument('--tlen',
                        dest="tlen_option",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='tlen option')    
    parser.add_argument('-o',
                        dest='out_fname',
                        default='bc_out',
                        type=str,
                        help='output prefix filename')
    
    args = parser.parse_args()

    if not args.ref_fname:
        print("Error: there is no reference file input.", file=sys.stderr)
        sys.exit(1)

    genome_size = Helper_Py3.genome_sizes(args.ref_fname)

    win_size = args.win_size
    bin_GC = None
    if args.GC_option:
        bin_GC = Helper_Py3.bin_gc_from_fasta(args.ref_fname, win_size)
    
    chr_list = []
    if not args.chr_list:
        chr_list = sorted(genome_size.keys(), key=Helper_Py3.chr_key)
    else:
        chr_list = sorted(args.chr_list, key=Helper_Py3.chr_key)

    bin_count (args.fnames,
               genome_size,
               win_size,
               args.skip_zero,
               bin_GC,
               args.tlen_option,
               args.mm_cutoff,
               args.min_len,
               args.max_len,
               chr_list,
               args.out_fname
               )

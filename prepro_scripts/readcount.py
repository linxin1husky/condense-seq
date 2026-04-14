import sys
from argparse import ArgumentParser
import gzip
import Helper_Py3_compat as Helper_Py3


def _collect_input_labels(fnames):
    data, label = [], []
    for fname in fnames:
        data.append(fname)
        label.append(fname.rsplit('/', 1)[-1].split('.')[0])
    return data, label


def _record_mid_count(chr_profile, ref_id, mid_pos, name):
    chrom_profile = chr_profile.setdefault(ref_id, {})
    pos_counts = chrom_profile.setdefault(mid_pos, {})
    pos_counts[name] = pos_counts.get(name, 0) + 1


def read_count (fnames,
                genome_size,
                min_len,
                max_len,
                mm_cutoff,
                chr_list,
                scale,
                out_fname):

    # gather whole file names
    data, label = _collect_input_labels(fnames)
    chr_choices = set(chr_list)

    # make genome start-end profile dictionary
    chr_profile = {}
    
    # open the sam file
    for i in range(len(data)):
        name = label[i]
        filename = data[i]
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

            start, end = Helper_Py3.fragment_bounds_from_record(record)
            mid_pos = (start + end) // 2
            _record_mid_count(chr_profile, record.rname, mid_pos, name)
            
    # summarize the output
    print("writing read count file", file=sys.stderr)
    
    f = gzip.open(out_fname + '_count.gtab.gz', 'wt', encoding='utf-8')
    s = 'Chromosome\tPosition'
    for i in range(len(label)):
        s += '\t' + label[i]
    print(s, file=f)

    for chr in sorted(chr_profile.keys(), key=Helper_Py3.chr_key):
        for i in sorted(chr_profile[chr]):
            counts = []
            for j in range(len(label)):
                name = label[j]
                count = chr_profile[chr][i].get(name, 0)
                counts.append(count)
            if sum(counts) == 0:
                continue
            s = chr + "\t" + str(i) + '\t'
            s += '\t'.join([str(round(count*scale, 5)) for count in counts])
            print(s, file=f)
        
    f.close()
    print("Done", file=sys.stderr)
        

if __name__ == '__main__':
    parser = ArgumentParser(description='Count reads along the genome')
    parser.add_argument(metavar='-f1',
                        dest="fnames",
                        type=str,
                        nargs='+',
                        help='SAM/Bam filenames')
    parser.add_argument(metavar='-x',
                        dest='ref_fname',
                        type=str,
                        help='reference sequence filename')
    parser.add_argument('-m',
                        dest="mm_cutoff",
                        type=int,
                        default=10,
                        help='mismatch cut-off in bp')
    parser.add_argument('--min',
                        dest="min_len",
                        type=int,
                        nargs='?',
                        default=0,
                        const=120,
                        help='minimum length for selection in bp')
    parser.add_argument('--max',
                        dest="max_len",
                        type=int,
                        nargs='?',
                        default=sys.maxsize,
                        const=170,
                        help='maximum length for selection in bp')
    parser.add_argument('--chr',
                        dest="chr_list",
                        type=str,
                        nargs='+',
                        help='target chromosome list')
    parser.add_argument('--scale',
                        dest="scale",
                        type=float,
                        default=1.0,
                        help='scaling factor')
    parser.add_argument('-o',
                        dest='out_fname',
                        default='output',
                        type=str,
                        help='output prefix filename')
    
    args = parser.parse_args()

    # get length for each chromosome
    if not args.ref_fname:
        print("Error: there is no reference file input.", file=sys.stderr)
        sys.exit(1)
        
    genome_size = Helper_Py3.genome_sizes(args.ref_fname)

    if not args.chr_list:
        chr_list = sorted(genome_size.keys(), key=Helper_Py3.chr_key)
    else:
        chr_list = sorted(args.chr_list, key=Helper_Py3.chr_key)

    #print chr_list

    read_count(
        args.fnames,
        genome_size,
        args.min_len,
        args.max_len,
        args.mm_cutoff,
        chr_list,
        args.scale,
        args.out_fname,
    )

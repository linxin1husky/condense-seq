import sys
from argparse import ArgumentParser
import math
import gzip
import Helper_Py3_compat as Helper_Py3


def _collect_input_labels(fnames):
    data, label = [], []
    for fname in fnames:
        data.append(fname)
        label.append(fname.rsplit('/', 1)[-1].split('.')[0])
    return data, label


def _record_mid_scores(chr_mid, ref_id, name, midscore_list):
    chrom_mid = chr_mid.setdefault(ref_id, {})
    for mid, score in midscore_list:
        pos_scores = chrom_mid.setdefault(mid, {})
        pos_scores[name] = pos_scores.get(name, 0) + score


def NCP_occ (fnames,
             genome_size,
             min_len,
             max_len,
             mm_cutoff,
             NCP_len,
             sigma,
             skip_zero,
             chr_list,
             out_fname):

    # gather whole file names
    data, label = _collect_input_labels(fnames)
    chr_choices = set(chr_list)

    # make genome start-end profile dictionary
    chr_mid = {}
    
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

            midscore_list = Helper_Py3.record_center_positions(record, even_policy="split")
            _record_mid_scores(chr_mid, record.rname, name, midscore_list)

    chr_profile = {}
    # mark start-end position of nucleosome
    # no smoothing
    if NCP_len:
        for chr in sorted(chr_mid.keys(), key=Helper_Py3.chr_key):
            for pos in chr_mid[chr]:
                st = pos - NCP_len // 2
                ed = pos + NCP_len // 2 + 1
                for name in chr_mid[chr][pos]:
                    score = chr_mid[chr][pos][name]
                    if chr not in chr_profile:
                        chr_profile[chr] = {}
                    if st not in chr_profile[chr]:
                        chr_profile[chr][st] = {}
                    if ed not in chr_profile[chr]:
                        chr_profile[chr][ed] = {}
                    if name not in chr_profile[chr][st]:
                        chr_profile[chr][st][name] = 0
                    if name not in chr_profile[chr][ed]:
                        chr_profile[chr][ed][name] = 0
                    chr_profile[chr][st][name] += score
                    chr_profile[chr][ed][name] -= score
                    
    # gaussain smoothing
    if sigma:
        for chr in sorted(chr_mid.keys(), key=Helper_Py3.chr_key):
            for pos in chr_mid[chr]:
                for name in chr_mid[chr][pos]:
                    score = chr_mid[chr][pos][name]
                    for i in range(int(math.floor(score))):
                        offset = int(round(sigma* math.sqrt(2*math.log(float(score)/(i+1)))))
                        st = max(0, pos - offset)
                        ed = min(genome_size[chr], pos + offset + 1)
                        if chr not in chr_profile:
                            chr_profile[chr] = {}
                        if st not in chr_profile[chr]:
                            chr_profile[chr][st] = {}
                        if ed not in chr_profile[chr]:
                            chr_profile[chr][ed] = {}
                        if name not in chr_profile[chr][st]:
                            chr_profile[chr][st][name] = 0
                        if name not in chr_profile[chr][ed]:
                            chr_profile[chr][ed][name] = 0
                        chr_profile[chr][st][name] += 1
                        chr_profile[chr][ed][name] -= 1
                
    # summarize the output
    print("writing occupancy file", file=sys.stderr)

    f = gzip.open(out_fname + '_occ.gtab.gz', 'wt', encoding='utf-8', newline='\n')
    s = 'Chromosome\tPosition'
    for i in range(len(label)):
        s += '\t' + label[i]
    print(s, file=f)

    for chr in sorted(chr_profile.keys(), key=Helper_Py3.chr_key):
        previous = [0]*len(label)
        if skip_zero:
            start = min(chr_profile[chr])
            end = max(chr_profile[chr]) + 1
        else:
            start = 0
            end = genome_size[chr]
        for i in range(start, end):
            #s = str(ID) + "\t" + chr + "\t" + str(i)
            s = chr + "\t" + str(i)
            for j in range(len(label)):
                name = label[j]
                past = previous[j]
                current = past + chr_profile[chr].get(i, {}).get(name, 0)
                s += "\t" + str(current)
                previous[j] = current
            if skip_zero and sum(previous) == 0:
                continue
            print(s, file=f)
            #ID += 1
        
    f.close()
    print("Done", file=sys.stderr)
        

if __name__ == '__main__':
    parser = ArgumentParser(description='Calculate occupancy along the genome')
    parser.add_argument(metavar='-f1',
                        dest="fnames",
                        type=str,
                        nargs='+',
                        help='SAM/Bam filenames')
    parser.add_argument(metavar='-x',
                        dest='ref_fname',
                        type=str,
                        help='reference sequence filename')
    parser.add_argument('--min',
                        dest="min_len",
                        type=int,
                        default=126,
                        help='minimum length for selection in bp, default:126bp')
    parser.add_argument('--max',
                        dest="max_len",
                        type=int,
                        default=184,
                        help='maximum length for selection in bp, default:184bp')
    parser.add_argument('-m',
                        dest="mm_cutoff",
                        type=int,
                        default=10,
                        help='mismatch cut-off in bp')
    parser.add_argument('--Nlen',
                        dest="NCP_len",
                        type=int,
                        default=73,
                        help='Mono-nucleosomal occupancy length in bp, default:73bp')
    parser.add_argument('--sigma',
                        dest="sigma",
                        type=int,
                        nargs='?',
                        const=20,
                        default=None,
                        help='Gaussian smoothing binwidth, default:20bp')
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
                        help='target chromosome list')
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

    if args.sigma:
        NCP_len = None
    else:
        NCP_len = args.NCP_len

    NCP_occ (args.fnames,
             genome_size,
             args.min_len,
             args.max_len,
             args.mm_cutoff,
             NCP_len,
             args.sigma, 
             args.skip_zero,
             chr_list,
             args.out_fname
             )

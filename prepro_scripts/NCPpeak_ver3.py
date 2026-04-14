import sys
import argparse
import gzip
import Helper_Py3_compat as Helper_Py3


def _collect_input_labels(fnames):
    data, label = [], []
    for fname in fnames:
        data.append(fname)
        label.append(fname.rsplit('/', 1)[-1].split('.')[0])
    return data, label


def _record_ncp_scores(chr_ncp, ref_id, name, ncp_scores):
    chrom_scores = chr_ncp.setdefault(ref_id, {})
    for ncp_pos, score in ncp_scores:
        pos_scores = chrom_scores.setdefault(ncp_pos, {})
        pos_scores[name] = pos_scores.get(name, 0) + score


def _write_gtab(path, chr_scores, labels):
    with gzip.open(path, 'wt', encoding='utf-8', newline='\n') as f:
        print('Chromosome\tPosition\t' + '\t'.join(labels), file=f)
        for chrom in sorted(chr_scores.keys(), key=Helper_Py3.chr_key):
            for pos in sorted(chr_scores[chrom].keys()):
                row = [chrom, str(pos)]
                score_map = chr_scores[chrom][pos]
                for name in labels:
                    row.append(str(score_map.get(name, 0)))
                print('\t'.join(row), file=f)


def NCP_peak (fnames,
              min_len,
              max_len,
              mm_cutoff,
              NCP_len,
              overlap,
              skip_zero,
              chr_list,
              out_fname):

    # gather whole file names
    data, label = _collect_input_labels(fnames)
    chr_choices = set(chr_list) if chr_list else None

    # make genome counts dictionary
    chr_NCP = {}
       
    # open the sam file
    for i in range(len(data)):
        name = label[i]
        #chr_cov = out_cov[i]
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

            ncp_scores = Helper_Py3.record_center_positions(record, even_policy="split")
            _record_ncp_scores(chr_NCP, record.rname, name, ncp_scores)

    # Select NCP positions with highest score and minimized overlapping.
    print("peak calling: filtering NCP positions", file=sys.stderr)
        
    # find the index where the target would be inserted in right order
    def binary_insert (sortlist, target):
        st, ed = 0, len(sortlist)-1
        while st <= ed:
            mid = (st+ed) // 2
            if sortlist[mid] == target:
                return mid
            elif sortlist[mid] > target:
                ed = mid - 1 
            elif sortlist[mid] < target:
                st = mid + 1
        return st

    chr_peak = {}
    for chr in sorted(chr_NCP.keys(), key=Helper_Py3.chr_key):
        name_temp = {}
        for NCPpos in sorted(chr_NCP[chr].keys()):
            for name in label:
                score = chr_NCP[chr][NCPpos].get(name)
                if score is None:
                    continue
                if name not in name_temp:
                    name_temp[name] = []
                name_temp[name].append([score, NCPpos])
        for name in name_temp:
            selected = []
            temp = sorted(name_temp[name], key=lambda x: x[0], reverse=True)
            for k in range(len(temp)):
                score, NCPpos = temp[k]
                idx = binary_insert (selected, NCPpos)
                neighbor = selected[max(0,idx-1):min(idx+1,len(selected))]
                check = True
                for pos in neighbor:
                    if abs(pos - NCPpos) < NCP_len - overlap:
                        check = False
                        break
                if check:
                    selected.insert(idx, NCPpos)
                    if chr not in chr_peak:
                        chr_peak[chr] = {}
                    if NCPpos not in chr_peak[chr]:
                        chr_peak[chr][NCPpos] = {}
                    assert name not in chr_peak[chr][NCPpos]
                    chr_peak[chr][NCPpos][name] = score


    # summarize the output
    print("writing NPS file", file=sys.stderr)
    _write_gtab(out_fname + '_NPS.gtab.gz', chr_NCP, label)

    
    print("writing peak file", file=sys.stderr)
    _write_gtab(out_fname + '_peak.gtab.gz', chr_peak, label)

    print("Done", file=sys.stderr)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calling NCP peaks and scores')
    parser.add_argument('-f',
                        dest="fnames",
                        type=str,
                        nargs='+',
                        help='SAM/Bam filenames (last file used as control)')
    parser.add_argument('-m',
                        dest="mm_cutoff",
                        type=int,
                        default=10,
                        help='mismatch cut-off in bp')
    parser.add_argument('--min',
                        dest="min_len",
                        type=int,
                        default=126,
                        help='minimum length for selection in bp (default:126bp)')
    parser.add_argument('--max',
                        dest="max_len",
                        type=int,
                        default=184,
                        help='maximum length for selection in bp (default:184bp)')
    parser.add_argument('--Nlen',
                        dest="NCP_len",
                        type=int,
                        default=147,
                        help='Mono-nucleosomal length in bp')
    parser.add_argument('--ovlap',
                        dest="overlap",
                        type=int,
                        default=40,
                        help='maximum allowed overlap between NCPS in bp')
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
    
    chr_list = []
    if not args.chr_list:
        chr_list = None
    else:
        chr_list = sorted(args.chr_list, key=Helper_Py3.chr_key)

    NCP_peak (args.fnames,
              args.min_len,
              args.max_len,
              args.mm_cutoff,
              args.NCP_len,
              args.overlap,
              args.skip_zero,
              chr_list,
              args.out_fname
              )

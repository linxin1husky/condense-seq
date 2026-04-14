import sys
import argparse
import gzip
import Helper_Py3_compat as Helper_Py3

"""
Step 1 of condense-seq pipeline after bowtie2 alignment.
"""


def _collect_input_labels(fnames):
    data, label = [], []
    for fname in fnames:
        data.append(fname)
        label.append(fname.rsplit('/', 1)[-1].split('.')[0])
    return data, label


def _record_profile_delta(chr_profile, ref_id, name, start, end):
    chrom_profile = chr_profile.setdefault(ref_id, {})
    start_scores = chrom_profile.setdefault(start, {})
    end_scores = chrom_profile.setdefault(end, {})
    start_scores[name] = start_scores.get(name, 0) + 1
    end_scores[name] = end_scores.get(name, 0) - 1

def NCP_count (fnames,
               genome_size,
               min_len,
               max_len,
               mm_cutoff,
               skip_zero,
               chr_list,
               out_fname):

    # gather whole file names
    data, label = _collect_input_labels(fnames)
    chr_choices = set(chr_list)

    # make genome start-end profile dictionary
    chr_profile = {}
    
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

            start, end = Helper_Py3.fragment_bounds_from_record(record)
            _record_profile_delta(chr_profile, record.rname, name, start, end)
            
    # summarize the output
    print("writing coverage file", file=sys.stderr)
    
    f = gzip.open(out_fname + '_cov.gtab.gz', 'wt', encoding='utf-8', newline='\n')
    s = 'Chromosome\tPosition'
    for i in range(len(label)):
        s += '\t' + label[i]
    print(s, file=f)

    for chr in sorted(chr_profile.keys(), key=Helper_Py3.chr_key):
        previous = [0]*len(label)
        if skip_zero:
            # start = min(chr_profile[chr])
            end = max(chr_profile[chr]) + 1
        else:
            # start = 0
            end = genome_size[chr]
        for i in range(end):
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
        
    f.close()
    print("Done", file=sys.stderr)
        

if __name__ == '__main__':
    str2bool = Helper_Py3.str2bool

    parser = argparse.ArgumentParser(
        description='Calculate coverage along the genome',
        epilog="""Example of usage: python3 condense-seq/prepro_scripts/coverage_ver3.py 
        -f <input BAM 1> <input BAM 2> <input BAM 3> -x <Reference genome in FASTA> --chr chr1 --skip -o <Output path prefix>""")
    parser.add_argument('-f',
                        dest="fnames",
                        type=str,
                        nargs='+',
                        help='SAM/Bam filenames')
    parser.add_argument('-x',
                        dest='ref_fname',
                        type=str,
                        help='reference FASTA sequence filename')
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
    parser.add_argument('--skip',
                        dest="skip_zero",
                        type=str2bool,
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

    NCP_count (args.fnames,
               genome_size,
               args.min_len,
               args.max_len,
               args.mm_cutoff,
               args.skip_zero,
               chr_list,
               args.out_fname
               )

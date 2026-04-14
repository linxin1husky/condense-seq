import sys
import argparse 
import Helper_Py3_compat as Helper_Py3


def _collect_input_labels(fnames):
    data, label = [], []
    for fname in fnames:
        data.append(fname)
        label.append(fname.split('/')[-1].split('.')[0])
    return data, label

def NCP_count (fnames,
               mm_cutoff,
               chr_list,
               out_fname):

    # gather whole file names
    data, label = _collect_input_labels(fnames)
    chr_choices = set(chr_list) if chr_list else None

    # make genome read length dictionary
    output_list = []
        
    # open the sam file
    for i in range(len(data)):
        name = label[i]
        filename = data[i]
        print("reading %s" % (filename), file=sys.stderr) 

        chr_rlen = {}
        for cols in Helper_Py3.iter_samtools_view(filename, extra_args=["-F", "0x10"]):
            record = Helper_Py3.parse_sam_record(cols)
            if not Helper_Py3.passes_common_paired_filters(
                record,
                chr_choices=chr_choices,
                mm_cutoff=mm_cutoff,
                require_proper_pair=True,
            ):
                continue

            # record read length
            rlen = abs(record.tlen)
            if record.rname not in chr_rlen:
                chr_rlen[record.rname] = {}
            if rlen not in chr_rlen[record.rname]:
                chr_rlen[record.rname][rlen] = 0
            chr_rlen[record.rname][rlen] += 1
            
        output_list.append(chr_rlen)

    # summarize the output
    print("writing rlen file", file=sys.stderr)

    for i in range(len(output_list)):
        chr_rlen = output_list[i]
        name = label[i]
        
        f = open(out_fname + '_' + name + '_rlen.txt', 'w')
        s = 'Chromosome\tReadLength\tCounts'
        print(s, file=f)

        for chr in sorted(chr_rlen.keys(), key=Helper_Py3.chr_key):
            for rlen in sorted(chr_rlen[chr].keys()):
                count = chr_rlen[chr][rlen]
                s = chr + "\t" + str(rlen) + "\t" + str(count)
                print(s, file=f)

        f.close()

    print("Done", file=sys.stderr)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get read length distribution')
    parser.add_argument('-f1',
                        dest="fnames",
                        type=str,
                        nargs='+',
                        help='SAM/Bam filenames (last file used as control)')
    parser.add_argument('-m',
                        dest="mm_cutoff",
                        type=int,
                        default=10,
                        help='mismatch cut-off in bp')
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

    NCP_count (args.fnames,
               args.mm_cutoff,
               chr_list,
               args.out_fname
               )

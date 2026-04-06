import sys, subprocess, re
from argparse import ArgumentParser
import gzip
import Helper_Py3_compat as Helper_Py3

def read_count (fnames,
                genome_size,
                min_len,
                max_len,
                mm_cutoff,
                chr_list,
                scale,
                out_fname):

    # gather whole file names
    data, label = [], []
    for fname in fnames:
        data.append(fname)
        label.append(fname.rsplit('/', 1)[-1].split('.')[0])

    # make genome start-end profile dictionary
    chr_profile = {}
    
    # open the sam file
    for i in range(len(data)):
        name = label[i]
        filename = data[i]
        print("reading %s" % (filename), file=sys.stderr)
        samtools_proc = subprocess.Popen(
            ["samtools", "view", "-F 0x10", filename],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        )

        for line in samtools_proc.stdout:
            # skip the header
            if line.startswith('@'):
                continue
            
            cols = line.strip().split()
            _, flag, ref_id, pos, _, cigar_str = cols[:6]

            tlen = int(cols[8])
            flag, pos = int(flag), int(pos)
            ref_id = ref_id.strip()
            pos-=1

            # non target chromosome
            if ref_id not in chr_list:
                continue

            # invalid: mapping failure
            if pos < 0:
                #type = 'invalid:mutant'
                continue
            if flag & 0x4 != 0:
                #type = 'invalid:mutant'
                continue
            if ref_id == '*':
                continue

            # invalid: mate pairing failture
            if flag & 0x8 != 0:
                #type = 'invalid:mutant'
                continue
            if flag & 0x2 == 0:
                continue

            # invalid: ambiguous mapping
            #mapQ = float(mapQ)
            #if mapQ < 10:
            #    type = 'invalid:multimap'

            AS,NM,MD = None, None, None
            for i in range(11, len(cols)):
                col = cols[i]
                if col.startswith('AS'):
                    AS = int(col[5:])
                elif col.startswith('NM'):
                    NM = int(col[5:])
                elif col.startswith('MD'):
                    MD = col[5:]

            # invalid: too large edit distance 
            if NM > mm_cutoff:
                #type = 'invalid:mutant'
                continue

            # invalid: not right size of read 
            if abs(tlen) < min_len or abs(tlen) > max_len:
                #type = 'invalid:mutant'
                continue

            # get NCP position and coverage range
            if tlen > 0: # left read
                end_pos = pos + tlen
            else: # right read
                end_pos = pos
                cigar_str=re.split('(\d+)',cigar_str)[1:]
                for i in range(len(cigar_str)//2):
                    s = cigar_str[2*i+1]
                    num = int(cigar_str[2*i])
                    if s == 'M' or s == 'D':
                        end_pos += num
                pos = end_pos + tlen

            mid_pos = (pos + end_pos) // 2

            # record mid positions of reads
            if ref_id not in chr_profile:
                chr_profile[ref_id] = {}
            if mid_pos not in chr_profile[ref_id]:
                chr_profile[ref_id][mid_pos] = {}
            if name not in chr_profile[ref_id][mid_pos]:
                chr_profile[ref_id][mid_pos][name] = 0
            chr_profile[ref_id][mid_pos][name] += 1
            
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
                try:
                    count = chr_profile[chr][i][name]
                except:
                    count = 0
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
                        help='tagert chromosome list')
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

import sys, subprocess, re
import argparse
import gzip
import Helper_Py3

"""
Step 1 of condense-seq pipeline after bowtie2 alignment.
"""

def NCP_count (fnames,
               genome_size,
               min_len,
               max_len,
               mm_cutoff,
               skip_zero,
               chr_list,
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
        #chr_cov = out_cov[i]
        filename = data[i]
        print("reading %s" % (filename), file=sys.stderr)
        samtools_proc = subprocess.Popen(["samtools",  "view", "-F 0x10", filename], 
                                         stdout=subprocess.PIPE, 
                                         stderr=subprocess.DEVNULL,
                                         text=True)

        for line in samtools_proc.stdout:
            # skip the header
            if line.startswith('@'):
                continue
            
            cols = line.strip().split()
            read_id, flag, ref_id, pos, mapQ, cigar_str = cols[:6]
            read_id=":".join(read_id.split(':')[3:7])

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
                cigar_str=re.split(r'(\d+)',cigar_str)[1:]
                for i in range(len(cigar_str)/2):
                    s = cigar_str[2*i+1]
                    num = int(cigar_str[2*i])
                    if s == 'M' or s == 'D':
                        end_pos += num
                pos = end_pos + tlen

            # record start-end positions of reads
            if ref_id not in chr_profile:
                chr_profile[ref_id] = {}
            if pos not in chr_profile[ref_id]:
                chr_profile[ref_id][pos] = {}
            if name not in chr_profile[ref_id][pos]:
                chr_profile[ref_id][pos][name] = 0
            if end_pos not in chr_profile[ref_id]:
                chr_profile[ref_id][end_pos] = {}
            if name not in chr_profile[ref_id][end_pos]:
                chr_profile[ref_id][end_pos][name] = 0
            chr_profile[ref_id][pos][name] += 1
            chr_profile[ref_id][end_pos][name] -= 1
            
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
                try:
                    current = past + chr_profile[chr][i][name]
                except:
                    current = past + 0
                s += "\t" + str(current)
                previous[j] = current
            if skip_zero and sum(previous) == 0:
                continue
            print(s, file=f)
        
    f.close()
    print("Done", file=sys.stderr)
        

if __name__ == '__main__':
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

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
                        help='tagert chromosome list')
    parser.add_argument('-o',
                        dest='out_fname',
                        default='output',
                        type=str,
                        help='output prefix filename')
    
    args = parser.parse_args()

    # get length for each chromosome
    genome_size = {}
    if not args.ref_fname:
        print("Error: there is no reference file input.", file=sys.stderr)
        sys.exit(1)
        
    for line in Helper_Py3.gzopen(args.ref_fname):
        line = line.strip()
        if line.startswith('>'):
            key = line.split()[0][1:]
            assert key not in genome_size
            genome_size[key] = 0
            continue
        genome_size[key] += len(line)

    chr_list = []
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

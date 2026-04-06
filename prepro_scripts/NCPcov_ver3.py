import sys, argparse, math, gzip
import Helper_Py3_compat as Helper_Py3

def NCP_count (peak_fname,
               cov_fname,
               peak_choice,
               NCP_len,
               chr_list,
               out_fname):

    def _store_cov(store, chrom, dyad, covs, labels):
        if chrom not in store:
            store[chrom] = {}
        store[chrom][dyad] = {labels[i]: covs[i] for i in range(len(labels))}

    # read peak file
    print("reading peak file", file=sys.stderr)
    chr_peak = {}
    First = True
    
    for line in Helper_Py3.open_any(peak_fname, "rt"):
        cols = line.strip().split()
        if First:
            label = cols[2:]
            First = False
            continue
        chr, pos = cols[0], int(cols[1])
        if chr_list and chr not in chr_list:
            continue
        if chr not in chr_peak:
            chr_peak[chr] = {}
        if pos not in chr_peak[chr]:
            chr_peak[chr][pos] = {}
        counts = cols[2:]
        for i in range(len(label)):
            name = label[i]
            count = float(counts[i])
            if count <=0:
                continue
            if name not in chr_peak[chr][pos]:
                chr_peak[chr][pos][name] = 0
            chr_peak[chr][pos][name] += count

    # make a list of target NCP peaks
    if peak_choice == 'input':
        target = [label[-1]] # target only control peaks
    elif peak_choice == 'all':
        target = label # target all peaks
        
    chr_st, chr_ed = {}, {}
    for chr in sorted(chr_peak.keys(), key=Helper_Py3.chr_key):
        for NCPpos in sorted(chr_peak[chr].keys()):
            for name in target:
                try:
                    score = chr_peak[chr][NCPpos][name]
                    if chr not in chr_st:
                        chr_st[chr] = []
                        chr_ed[chr] = []
                    st = NCPpos - NCP_len // 2
                    ed = NCPpos + NCP_len // 2
                    if st not in chr_st:
                        chr_st[chr].append(st)
                        chr_ed[chr].append(ed)
                except:
                    continue           

    # reading coverage file and get NCP coverage
    print("reading coverage file", file=sys.stderr)
    chr_Ncov = {}
    NCPcovs = []
    
    First = True
    order = None
    prev_chr = None

    for line in Helper_Py3.open_any(cov_fname, "rt"):
        cols = line.strip().split()
        if First:
            label = cols[2:]
            First = False
            continue
        chr, pos = cols[0], int(cols[1])
        if chr not in chr_st:
            continue
        if prev_chr is not None and chr != prev_chr:
            while chr_st[prev_chr]:
                chr_st[prev_chr].pop(0)
                NCPcovs.insert(0, [0]*len(label))
            while chr_ed[prev_chr]:
                ed = chr_ed[prev_chr].pop(0)
                dyad = ed - NCP_len // 2
                covs = NCPcovs.pop()
                _store_cov(chr_Ncov, prev_chr, dyad, covs, label)
            prev_chr = chr
        if prev_chr is None:
            prev_chr = chr
        if order == None:
            print(chr + " pos " + str(pos) +" is reading", file=sys.stderr)
            order = int(math.log10(max(pos, 1)))
        elif int(math.log10(max(pos, 1))) > order:
            print(chr + " pos " + str(pos) +" is reading", file=sys.stderr)
            order += 1
        while chr_st[chr] and pos >= chr_st[chr][0]:
            chr_st[chr].pop(0)
            NCPcovs.insert(0, [0]*len(label))
        while chr_ed[chr] and pos > chr_ed[chr][0]:
            ed = chr_ed[chr].pop(0)
            dyad = ed - NCP_len // 2
            covs = NCPcovs.pop()
            _store_cov(chr_Ncov, chr, dyad, covs, label)
        if len(NCPcovs) <= 0:
            continue
        counts = cols[2:]
        for k in range(len(NCPcovs)):
            for u in range(len(label)):
                count = int(counts[u])
                NCPcovs[k][u] += count
    while chr_st[chr]:
        chr_st[chr].pop(0)
        NCPcovs.insert(0, [0]*len(label))
    while chr_ed[chr]:
        ed = chr_ed[chr].pop(0)
        dyad = ed - NCP_len // 2
        covs = NCPcovs.pop()
        _store_cov(chr_Ncov, chr, dyad, covs, label)
    
    # write NCP coverage file
    print("writing NCP coverage file", file=sys.stderr)
    
    f = gzip.open(out_fname + '_Ncov.gtab.gz', 'wt', encoding='utf-8', newline='\n')
    s = 'Chromosome\tPosition'
    for i in range(len(label)):
        s += '\t' + label[i]
    print(s, file=f)

    for chr in sorted(chr_Ncov.keys(), key=Helper_Py3.chr_key):
        for NCPpos in sorted(chr_Ncov[chr].keys()):
            s = chr + "\t" + str(NCPpos)
            for name in label:
                s += "\t" + str(chr_Ncov[chr][NCPpos].get(name, 0))
            print(s, file=f)
        
    f.close()

    print("Done", file=sys.stderr)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate NCP coverage')
    parser.add_argument('--peak',
                        dest="peak_fname",
                        type=str,
                        help='peak gtab file')
    parser.add_argument('--cov',
                        dest='cov_fname',
                        type=str,
                        help='coverage gtab file')
    parser.add_argument('--peak-choice',
                        dest="peak_choice",
                        type=str,
                        default='input',
                        help='NCP peak data choice (default:input control only, "all":all data)')
    parser.add_argument('--Nlen',
                        dest="NCP_len",
                        type=int,
                        default=171,
                        help='Mono-nucleosomal window in bp')
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

    if not args.chr_list:
        chr_list = None
    else:
        chr_list = sorted(args.chr_list, key=Helper_Py3.chr_key)

    NCP_count (args.peak_fname,
               args.cov_fname,
               args.peak_choice,
               args.NCP_len,
               chr_list,
               args.out_fname
               )

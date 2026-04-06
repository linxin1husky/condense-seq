import sys, math, gzip
from argparse import ArgumentParser
import Helper_Py3_compat as Helper_Py3

def Bin_data (data_fname,
              chr_size,
              bin_size,
              bin_step,
              bin_value,
              skip_zero,
              chr_list,
              out_fname):

    # reading data file and binning the data
    print >> sys.stderr, "reading data file"

    chr_Bdata = {}
    chr_Bcount = {}  
    
    First = True
    order = None
    for line in Helper_Py3.open_any(data_fname, "rt"):
        cols = line.strip().split()

        # find data type and range
        if First:
            _, col_st, _, _ = Helper_Py3.read_gtab_header(cols)
            labels = cols[col_st:]
            First = False
            continue

        if col_st == 2:
            chr, st, ed = cols[0], int(cols[1]), int(cols[1]) + 1
        else:
            chr, st, ed = cols[0], int(cols[1]), int(cols[2])
        
        if chr not in chr_list:
            continue

        if order == None:
            print >> sys.stderr, chr + " st " + str(st) +" is reading"
            order = int(math.log10(max(st, 1)))
        elif int(math.log10(max(st, 1))) > order:
            print >> sys.stderr, chr + " st " + str(st) +" is reading"
            order += 1

        # limit the input end by chromosome size
        ed = min(ed, chr_size[chr])

        # find bins overlap with data range
        idx = st // bin_step
        bst = bin_step*idx
        bed = bst + bin_size

        idx_st = idx
        while st >= bst and st < bed:
            if idx < idx_st:
                idx_st = idx
            idx -= 1
            if idx < 0:
                break
            bst = bin_step*idx
            bed = bst + bin_size

        idx_ed = (ed - 1) // bin_step
        assert idx_st <= idx_ed
        
        # binning the data
        values = [float(value) for value in cols[col_st:]]
        for idx in range(idx_st, idx_ed+1):
            bst = bin_step*idx
            bed = bst + bin_size
            a, b = max(st, bst), min(ed, bed)
            length = b - a

            if chr not in chr_Bdata:
                chr_Bdata[chr] = {}
            if idx not in chr_Bdata[chr]:
                chr_Bdata[chr][idx] = {}

            if bin_value == 'mean':
                if chr not in chr_Bcount:
                    chr_Bcount[chr] = {}
                if idx not in chr_Bcount[chr]:
                    chr_Bcount[chr][idx] = {}

            for name, value in zip(labels, values):
                if name not in chr_Bdata[chr][idx]:
                    chr_Bdata[chr][idx][name] = 0
                chr_Bdata[chr][idx][name] +=value*length

                if bin_value == 'mean':
                    if name not in chr_Bcount[chr][idx]:
                        chr_Bcount[chr][idx][name] = 0
                    chr_Bcount[chr][idx][name] +=1
                    
    # write bin data file
    print("writing bin data file", file=sys.stderr)
    
    f = gzip.open(out_fname + '_Bdata.gtab.gz', 'wt', encoding='utf-8', newline='\n')
    s = 'Chromosome\tStart\tEnd\t'
    s += '\t'.join(labels)
    print(s, file=f)

    for chr in chr_list:
        try:
            chr_Bdata[chr]
        except:
            continue
        
        idx = 0
        last_idx = chr_size[chr] // bin_step
        while idx <= last_idx:
            try:
                name_Bdata = chr_Bdata[chr][idx]
            except:
                idx +=1
                continue
            
            total = sum(name_Bdata.values())
            
            if skip_zero and total <=0:
                idx +=1
                continue

            Bvalues = []
            for name in labels:
                Bdata = name_Bdata[name]
                if bin_value == 'mean':
                    Bcount = chr_Bcount[chr][idx][name]
                    if Bcount <= 0:
                        Bvalue = 'NA'
                    else:
                        Bvalue = float(Bdata)/Bcount
                elif bin_value == 'sum':
                    Bvalue = Bdata
                Bvalues.append(Bvalue)

            bst = bin_step*idx
            bed = min(bst + bin_size, chr_size[chr])
            
            s = [chr, str(bst), str(bed)]
            s += [str(Bvalue) for Bvalue in Bvalues]

            s = '\t'.join(s)
            print(s, file=f)

            #ID +=1
            idx +=1
                        
    f.close()

    print("Done", file=sys.stderr)

    
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

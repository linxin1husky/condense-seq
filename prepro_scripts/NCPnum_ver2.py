import sys, argparse
import gzip
import Helper_Py3_compat as Helper_Py3

def _header_info(cols):
    data_type, col_st, col_ed, _ = Helper_Py3.read_gtab_header(cols)
    if data_type == "point":
        fields = ["Chromosome", "Position"]
    else:
        fields = ["Chromosome", "Start", "End"]
    return col_st, col_ed, fields

def NCP_number (tnum_fnames,
                tfname,
                input_mnum,
                chr_list,
                out_fname):
    """ Estimate physical number of molecules from coverage/count data and titration information"""
    # read titration file
    tnum_tfrac = Helper_Py3.read_titration(tfname)

    # sort by titration number
    tnums = sorted(tnum_fnames.keys())

    for i in range(len(tnums)):
        tnum = tnums[i]
        print("processing tnum %d" % (tnum), file=sys.stderr)

        mnum = input_mnum * tnum_tfrac[tnum] # molecule number of titraion

        # get total of all coverage/counts
        total_covs = []
        fnames = tnum_fnames[tnum]
        for fname in fnames:
            print("\t reading %s" % (fname.rsplit('/', 1)[-1]), file=sys.stderr)
                            
            First = True
            for line in Helper_Py3.gzopen(fname):
                line = line.strip()
                if not line:
                    continue
                cols = line.split()
                if First:
                    col_st, col_ed, _ = _header_info(cols)
                    if not total_covs:
                        total_covs = [0] * (col_ed - col_st)
                    else:
                        assert len(total_covs) == col_ed - col_st
                    First = False
                    continue

                # non target chromosome
                chr_name = cols[0]
                if chr_list and chr_name not in chr_list:
                    continue

                for i in range(col_st, col_ed):
                    total_covs[i - col_st] += float(cols[i])

        # convert to molecule number
        for fname in fnames:
            print("\t converting %s" % (fname.rsplit('/', 1)[-1]), file=sys.stderr)

            out_fname = fname.rsplit('_', 1)[0] + '_num.gtab.gz'
            f = gzip.open(out_fname, 'wt', encoding='utf-8')

            First = True
            for line in Helper_Py3.gzopen(fname):
                line = line.strip()
                if not line:
                    continue
                cols = line.split()
                if First:
                    col_st, col_ed, fields = _header_info(cols)
                    fields += cols[col_st:col_ed]
                    print('\t'.join(fields), file=f)
                    
                    First = False
                    continue

                # non target chromosome
                chr_name = cols[0]
                if chr_list and chr_name not in chr_list:
                    continue

                nums = []
                for i in range(col_st, col_ed):
                    cov = float(cols[i])
                    frac = cov / total_covs[i - col_st]
                    nums.append(str(int(round(mnum * frac))))

                print('\t'.join(cols[:col_st] + nums), file=f)

            f.close()
                    
    print("Done", file=sys.stderr)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Estimate physical number of molecules')
    parser.add_argument('-f',
                        dest="fnames_list",
                        type=str,
                        nargs='+',
                        action='append',
                        help='Ncov/bin/Bdata.gtab file list')
    parser.add_argument('-t',
                        dest='tfname',
                        type=str,
                        help='titration filename')
    parser.add_argument('--tnum',
                        dest="tnums",
                        type=int,
                        nargs='+',
                        help='titration number of each data')
    parser.add_argument('--mscale',
                        dest="mscale",
                        default = 1,
                        type=int,
                        help='total molecule number scale of input')
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

    # map tnum to filenames
    if not args.tnums:
        tnums = range(len(args.fnames_list))
    else:
        if len(args.fnames_list) != len(args.tnums):
            print("Error: mismatch of input file and titration number.", file=sys.stderr)
            sys.exit(1)
        tnums = args.tnums
        
    tnum_fnames = {}
    for tnum, fnames in zip(tnums, args.fnames_list):
        tnum_fnames[tnum] = fnames

    # list target chromosomes
    if not args.chr_list:
        chr_list = None
    else:
        chr_list = sorted(args.chr_list, key=Helper_Py3.chr_key)

    input_mnum = args.mscale * 1.6*(10**12)
    
    NCP_number (tnum_fnames,
                args.tfname,
                input_mnum,
                chr_list,
                args.out_fname
                )

import sys
from argparse import ArgumentParser
import math
import Helper_Py3_compat as Helper_Py3

def standardize (fnames,
                 chr_list,
                 colnums,
                 out_fname):

    count_list = None
    mean_list = None
    std_list = None
    fieldnames = None

    ## read all files 2 times to get statistics
    ## 1st: get mean, 2nd: get std
    for i in range(2):
        if i == 0:
            print("summing up data to compute means", file=sys.stderr)
        elif i == 1:
            print("summing up data to compute stds", file=sys.stderr)
            
        for fname in fnames:            
            print("\t reading %s" % (fname.rsplit('/', 1)[-1]), file=sys.stderr)

            First = True
            for line in Helper_Py3.open_any(fname, "rt"):
                line = line.strip()
                if not line:
                    continue
                cols = line.split()

                if First:
                    # check field names
                    if fieldnames is None:
                        fieldnames = cols
                    else:
                        # all files has same field names
                        assert cols == fieldnames

                    # set data range
                    _, col_st, col_ed, _ = Helper_Py3.read_gtab_header(cols)

                    if colnums is None:
                        colnums = list(range(col_st, col_ed))

                    if i == 0 and mean_list is None:
                        count_list = [0] * len(colnums)
                        mean_list = [0] * len(colnums)

                    if i == 1 and std_list is None:
                        std_list = [0] * len(colnums)

                    First = False
                    continue

                # non target chromosome
                chr_name = cols[0]
                if chr_list and chr_name not in chr_list:
                    continue

                # sum up the data
                for k, colnum in enumerate(colnums):
                    try:
                        value = float(cols[colnum])
                    except ValueError:
                        continue

                    if i == 0:
                        count_list[k] += 1
                        mean_list[k] += value
                    elif i == 1:
                        mean = mean_list[k]
                        std = (value - mean)**2
                        std_list[k] += std

        # take average of sums at last
        for k, _ in enumerate(colnums):
            count = count_list[k]
            if count <= 0:
                continue
            if i == 0:            
                mean_list[k] = float(mean_list[k]) / count
            elif i == 1:
                std_list[k] = math.sqrt(float(std_list[k]) / count)

    ##compute z-score and writing down output files
    print("data standardization and writing output", file=sys.stderr)
    for fname in fnames:
        out_fname = fname.rsplit('_', 1)[0] + '_zscore.gtab.gz'
        print("\t writing %s" % (out_fname.rsplit('/', 1)[-1]), file=sys.stderr)
        
        f = Helper_Py3.open_any(out_fname, "wt")
        
        First = True
        for line in Helper_Py3.open_any(fname, "rt"):
            line = line.strip()
            if not line:
                continue
            cols = line.split()

            if First:
                # set data range
                _, col_st, _, _ = Helper_Py3.read_gtab_header(cols)
                    
                row = cols[:col_st]
                for colnum in colnums:
                    row.append(cols[colnum])
                print('\t'.join(row), file=f)
                First = False
                continue

            # non target chromosome
            chr_name = cols[0]
            if chr_list and chr_name not in chr_list:
                continue

            # compute z-score and writing to output
            row = cols[:col_st]
            for k, colnum in enumerate(colnums):
                try:
                    value = float(cols[colnum])
                except ValueError:
                    zscore = 'NA'
                else:
                    mean = mean_list[k]
                    std = std_list[k]
                    if std <= 0:
                        zscore = 'NA'
                    else:
                        zscore = round(float(value-mean)/std, 5)
                row.append(str(zscore))

            print('\t'.join(row), file=f)

        f.close()
                    
    print("Done", file=sys.stderr)


if __name__ == '__main__':
    parser = ArgumentParser(description='standardize values in gtab file')
    parser.add_argument(metavar='-f',
                        dest="fnames",
                        type=str,
                        nargs='+',
                        help='gtab file list')
    parser.add_argument('--chr',
                        dest="chr_list",
                        type=str,
                        nargs='+',
                        help='target chromosome list')
    parser.add_argument('-c',
                        dest="colnums",
                        type=str,
                        nargs='+',
                        help='picked column numbers of files for standardization')
    parser.add_argument('-o',
                        dest='out_fname',
                        default='output',
                        type=str,
                        help='output prefix filename')

    args = parser.parse_args()

    # list target chromosomes
    if not args.chr_list:
        chr_list = None
    else:
        chr_list = sorted(args.chr_list, key=Helper_Py3.chr_key)

    # list target column numbers
    if not args.colnums:
        colnums = None
    else:
        colnums = [int(colnum) for colnum in args.colnums]

    standardize (args.fnames,
                 chr_list,
                 colnums,
                 args.out_fname
                 )

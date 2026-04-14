import sys
import argparse
import Helper_Py3_compat as Helper_Py3


def _coord_key(chrom):
    token = chrom[3:] if chrom.lower().startswith('chr') else chrom
    try:
        return int(token)
    except ValueError:
        return token

def concatenate (fnames,
                 colnums_list,
                 out_fname):

    # initialize the parameters
    files = []
    file_pts = []
    file_cols = []
    read_checks = []
    EOF_checks = []

    for i in range(len(fnames)):
        fname = fnames[i]
        file = Helper_Py3.open_any(fname, "rt")
        files.append(file)
        file_pts.append(None)
        file_cols.append(None)
        read_checks.append(True)
        EOF_checks.append(False)

    # start concatenate the files
    print("Concatenating files", file=sys.stderr)
    f = Helper_Py3.open_any(out_fname +'.gtab.gz', "wt")

    First = True
    while True:
        # read files line by line but skipping empty line
        for i in range(len(fnames)):
            file = files[i]
            read_check = read_checks[i]
            EOF_check = EOF_checks[i]

            while read_check and not EOF_check:
                line = file.readline()

                if not line:
                    EOF_checks[i] = True
                    break

                line = line.strip()
                if line:
                    cols = line.split()
                    file_cols[i] = cols

                    if not First:
                        chr = cols[0]
                        pos = cols[1:col_st]

                        chr = _coord_key(chr)

                        pos = [int(value) for value in pos]
                        pt = [chr] + pos
                        pt = tuple(pt)

                        file_pts[i] = pt
                    break

        # breakout the loop if all files at EOF
        if all(EOF_check == True for EOF_check in EOF_checks):
            break    

        # check the data fields in the first line of files
        if First:
            fields_list = []
            for i in range(len(fnames)):
                file_col = file_cols[i]
                colnums = colnums_list[i]

                # check data type and ranges
                data_type, col_st, col_ed, _ = Helper_Py3.read_gtab_header(file_col)

                # check the genomic fields
                if i == 0:
                    genomic_fields = file_col[:col_st]
                    fields_list += genomic_fields
                else:
                    assert file_col[:col_st] == genomic_fields # same genomic field

                # select data fields
                if len(colnums) <=0:
                    colnums = list(range(col_st, col_ed))
                    
                fields = [file_col[colnum] for colnum in colnums]
                fields_list += fields

            print('\t'.join(fields_list), file=f)
            First = False
            continue

        # get file pointer at minimum genomic position
        pts = []
        for i in range(len(file_pts)):
            file_pt = file_pts[i]
            EOF_check = EOF_checks[i]
            if EOF_check:
                continue
            pts.append(file_pt)
        data_pt = sorted(pts)[0]
        
        # read data of files and write to output
        data_list = []
        for i in range(len(fnames)):
            file_pt = file_pts[i]
            file_col = file_cols[i]
            read_check = read_checks[i]
            EOF_check = EOF_checks[i]
            colnums = colnums_list[i]

            if EOF_check:
                data = ['NA'] * len(colnums)
                data_list += data
                read_checks[i] = False

            else:
                if data_pt < file_pt:
                    data = ['NA'] * len(colnums)
                    data_list += data
                    read_checks[i] = False

                else:
                    assert data_pt == file_pt
                    data = [file_col[colnum] for colnum in colnums]
                    data_list += data
                    read_checks[i] = True

        #row = [str(ID)]
        row = []
        row += ['chr' + str(data_pt[0])]
        row += [str(value) for value in data_pt[1:]]
        row += data_list
        
        print('\t'.join(row), file=f)
        #ID +=1

    f.close()
    print("Done", file=sys.stderr)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='concatenate gtab files')
    parser.add_argument('-f',
                        dest="fnames",
                        type=str,
                        nargs='+',
                        help='gtab file list')
    parser.add_argument('-c',
                        dest="colnums_list",
                        type=str,
                        nargs='+',
                        help='colnum numbers of files to append (format: a1,a2 b1,b2 ...)')
    parser.add_argument('-o',
                        dest='out_fname',
                        default='output',
                        type=str,
                        help='output prefix filename')
    
    args = parser.parse_args()

    # parse the column numbers for each file
    if args.colnums_list:
        assert len(args.colnums_list) == len(args.fnames)
        colnums_list = []
        for colnums in args.colnums_list:
            colnums = [int(colnum) for colnum in colnums.split(',')]
            colnums_list.append(colnums)
    else:
        colnums_list = [[] for i in range(len(args.fnames))]
        
    concatenate (args.fnames,
                 colnums_list,
                 args.out_fname
                 )

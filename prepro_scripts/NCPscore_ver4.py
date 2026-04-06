import os, sys, subprocess, re
import argparse
import math
import gzip
import Helper_Py3_compat as Helper_Py3

# metric for condensability
def get_score (test, control, metric):
    if control <= 0:
        return "NA"
    ratio = float(test)/float(control)
    if metric == None:
        return ratio
    elif type(metric) == int:
        return -math.log(ratio, base=metric)
    else:
        assert metric == 'e'
        return -math.log(ratio)


def NCP_score (fnames_list,
               control_fnames,
               numc_choice,
               tfname,
               tnums,
               metric,
               chr_list,
               out_fname):

    if numc_choice:
        # read titration file
        tnum_tfrac = Helper_Py3.read_titration(tfname)

    # first-loop: get total coverage/counts for normalization
    # second-loop: compute condensabilty scores
    for u in range(2):

        if u == 0:
            print("summing data for normalization", file=sys.stderr)
            total_covs_list = []
        else:
            print("computing scores", file=sys.stderr)

        for k in range(len(fnames_list)):
            fnames = fnames_list[k]

            if u == 0:
                total_covs = []
                if k == 0:
                    control_covs = []
            else:
                total_covs = total_covs_list[k]
                if numc_choice:
                    tfrac = tnum_tfrac[tnums[k]]
                    offset = get_score(tfrac, 1, metric=metric)
                else:
                    if metric != None:
                        offset = 0
                    else:
                        offset = 1

            for fname, control_fname in zip(fnames, control_fnames):
                print(f"processing {(fname.rsplit('/', 1)[-1])}", file=sys.stderr)

                test_f = Helper_Py3.gzopen(fname)
                control_f = Helper_Py3.gzopen(control_fname)

                if u == 1:
                    out_fname = fname.rsplit('_', 1)[0] + '_score.gtab.gz'
                    f = gzip.open(out_fname, 'wt', encoding='utf-8')

                test_read, control_read = True, True
                test_EOF, control_EOF = False, False
                First = True
                #ID = 0
                while True:
                    # read test file line by line but skipping empty line
                    while test_read and not test_EOF:
                        test_line = test_f.readline()

                        if not test_line:
                            test_EOF = True
                            break

                        test_line = test_line.strip()
                        if test_line:
                            test_cols = test_line.split()

                            if not First:
                                #test_chr = test_cols[1]
                                #test_pos = test_cols[2:col_st]
                                test_chr = test_cols[0]
                                test_pos = test_cols[1:col_st]

                                # skip non-target chromosome
                                if chr_list and test_chr not in chr_list:
                                    continue

                                try:
                                    test_chr = int(test_chr[3:])
                                except:
                                    test_chr = test_chr[3:]

                                test_pos = [int(value) for value in test_pos]
                                test_pt = [test_chr] + test_pos
                                test_pt = tuple(test_pt)    

                            break

                    # read control file line by line but skipping empty line
                    while control_read and not control_EOF:
                        control_line = control_f.readline()

                        if not control_line:
                            control_EOF = True
                            break

                        control_line = control_line.strip()
                        if control_line:
                            control_cols = control_line.split()

                            if not First:
                                #control_chr = control_cols[1]
                                #control_pos = control_cols[2:col_st]
                                control_chr = control_cols[0]
                                control_pos = control_cols[1:col_st]

                                # skip non-target chromosome
                                if chr_list and control_chr not in chr_list:
                                    continue

                                try:
                                    control_chr = int(control_chr[3:])
                                except:
                                    control_chr = control_chr[3:]

                                control_pos = [int(value) for value in control_pos]
                                control_pt = [control_chr] + control_pos
                                control_pt = tuple(control_pt)

                            break

                    ## break out the loop if control files reach EOF
                    #if control_EOF:
                    #    break

                    # break out the loop if both files reach EOF
                    if test_EOF and control_EOF:
                        break

                    if First:
                        data_type, col_st, col_ed, _ = Helper_Py3.read_gtab_header(test_cols)
                        _, control_col_st, control_col_ed, _ = Helper_Py3.read_gtab_header(control_cols)
                        assert test_cols[:col_st] == control_cols[:control_col_st]  # same field
                        assert col_st == control_col_st and col_ed == control_col_ed

                        if u == 0:
                            if len(total_covs) == 0:
                                total_covs = [0] * (col_ed - col_st)
                                
                                if k == 0 and len(control_covs) == 0:
                                    control_covs = [0] * (col_ed - col_st)

                            else:
                                # all files has same data range
                                assert len(total_covs) == col_ed - col_st
                        else:
                            #print >> f, test_line
                            print('\t'.join(test_cols[:col_ed]), file=f)

                        First = False
                        continue

                    if test_EOF: # test data is empty
                        assert not control_EOF
                        data_pt = control_pt
                        test_data = [0] * (col_ed - col_st)
                        control_data = control_cols[col_st:col_ed]
                        test_read = False
                        control_read = True

                    elif control_EOF: # skip when control data is empty
                        assert not test_EOF
                        #data_pt = test_pt
                        #txest_data = test_cols[col_st:col_ed]
                        #control_data = [0] * (col_ed - col_st)
                        test_read = True
                        control_read = False
                        continue

                    else:
                        assert not test_EOF
                        assert not control_EOF

                        # skip when control data is empty
                        if test_pt < control_pt:
                            #data_pt = test_pt
                            #test_data = test_cols[col_st:col_ed]
                            #control_data = [0] * (col_ed - col_st)
                            test_read = True
                            control_read = False
                            continue

                        # test data is empty
                        if test_pt > control_pt:
                            data_pt = control_pt
                            test_data = [0] * (col_ed - col_st)
                            control_data = control_cols[col_st:col_ed]
                            test_read = False
                            control_read = True

                        # both test and control data are available
                        if test_pt == control_pt:
                            data_pt = test_pt
                            test_data = test_cols[col_st:col_ed]
                            control_data = control_cols[col_st:col_ed]
                            test_read = True
                            control_read = True

                    # skip when whole row of control is empty
                    if sum([float(count) for count in control_data]) <=0:
                        continue

                    # compare the test and control data to compute the score
                    s = []
                    for i in range(col_ed - col_st):
                        count = float(test_data[i])
                        control = float(control_data[i])

                        if control <= 0: # score is not defined when control data is empty
                            score = 'NA'
                            s.append(score)
                            continue

                        # add pseudo count
                        count+=1
                        control+=1

                        if u == 0:
                            total_covs[i] += count

                            if k == 0:
                                control_covs[i] += control

                        else:
                            # normalize by total
                            ncount = float(count)/total_covs[i]
                            ncontrol = float(control)/control_covs[i]

                            # get score
                            score = get_score(ncount, ncontrol, metric=metric)
                            if metric != None:
                                score += offset # offset correction
                            else:
                                score *= offset

                            s.append(str(round(score, 5)))

                    if u == 1:
                        #row = [str(ID)]
                        row = []
                        row += ['chr' + str(data_pt[0])]
                        row += [str(value) for value in data_pt[1:]]
                        row += s
                        print('\t'.join(row), file=f)
                        #ID +=1

                if u == 1:
                    f.close()

        if u == 0:
            total_covs_list.append(total_covs)

        print("", file=sys.stderr)   

    print("Done", file=sys.stderr)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get condensability scores')
    parser.add_argument('-f',
                        dest="fnames_list",
                        type=str,
                        nargs='+',
                        action='append',
                        help='Ncov/bin/Bdata.gtab file list')
    parser.add_argument('-i',
                        dest="control_fnames",
                        type=str,
                        nargs='+',
                        help='input control files (in same order of data files)')
    parser.add_argument('-t',
                        dest='tfname',
                        type=str,
                        help='titration filename')
    parser.add_argument('--tnum',
                        dest="tnums",
                        type=int,
                        nargs='+',
                        help='titration number for each data (in same order of data files)')
    parser.add_argument('--numc',
                        dest="numc_choice",
                        type=bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='correct scores using titration file')
    parser.add_argument('--metric',
                        dest='metric',
                        type=str,
                        default='-log[e]',
                        help='score metric (-log[base], default:-log[e], None:No-log)')    
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

    
    # check the control filenames
    for fnames in args.fnames_list:
        if len(fnames) != len(args.control_fnames):
            print("Error: mismatch of data and control files", file=sys.stderr)
            sys.exit(1)

            
    # check number correction option
    if args.numc_choice:
        numc_choice = True
        if not args.tnums or not args.tfname:
            print("No titration information", file=sys.stderr)
            print("starting score computation without correction", file=sys.stderr)
        # check tnums for filenames
        elif len(args.fnames_list) != len(args.tnums):
            print("Error: mismatch of file and titration number.", file=sys.stderr)
            sys.exit(1)
        else:
            print("starting score computation with correction", file=sys.stderr)
            
    else:
        print("starting score computation without correction", file=sys.stderr)
        numc_choice = False

    # check condensability score metric
    metric = args.metric
    if metric.startswith('-log'):
        metric = metric.split('[')[-1][:-1]
        try:
            metric = int(metric)
        except:
            assert metric == 'e'
    elif metric.startswith('None'):
        metric = None
    else:
        print("require proper format of metric argument", file=sys.stderr)
            
    # list target chromosomes
    chr_list = []
    if not args.chr_list:
        chr_list = None
    else:
        chr_list = sorted(args.chr_list, key=Helper_Py3.chr_key)
        
    NCP_score (args.fnames_list,
               args.control_fnames,
               numc_choice,
               args.tfname,
               args.tnums,
               metric,
               chr_list,
               args.out_fname
               )

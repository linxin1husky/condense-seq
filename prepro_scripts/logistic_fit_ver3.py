import sys
from argparse import ArgumentParser
import warnings
import Helper_Py3_compat as Helper_Py3
import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution

# 4-parameter logistic function (sigmoid type)
def sigmoid_func (x, top, rate, chalf, bottom):
    y = bottom + float(top-bottom)/(1+np.exp(rate*(x-chalf)))
    return y

# 4-parameter logistic function (Hill type)
def hill_func (x, top, rate, chalf, bottom):
    y = bottom + float(top-bottom)/(1.0 + (x/float(chalf))**rate)
    return y

# objective function to optimize parameters
def obj_func (parms, func, X, Y):
    warnings.filterwarnings("ignore")
    Y_pred = func(X, *parms)
    return np.sum((Y_pred - Y)**2)

# compute CP value from 4PL data
def get_CP (top, hill, chalf, bottom, percent):
    surv_frac = 1 - percent/100.0
    CP = chalf*(((float(top-bottom)/(surv_frac-bottom)) - 1)**(1.0/float(hill)))
    return CP

def logistic_fit (fnames,
                  tfname,
                  tnums,
                  model,
                  method,
                  min_rsq,
                  min_top,
                  max_top,
                  min_bottom,
                  max_bottom,
                  min_rate,
                  max_rate,
                  min_chalf,
                  max_chalf,
                  chr_list,
                  graph_option,
                  out_fname):
    
    # read titration file
    tnum_conc, tnum_tfrac = Helper_Py3.read_titration(
        tfname,
        return_conc=True,
        conc_col=0,
        frac_col=7,
        tnum_col=-1,
        delim="\t",
    )

    # covert titration number to physical concentration
    concs = [tnum_conc[tnum] for tnum in tnums]

    # set indices of data by the order of concentration
    conc_idx = sorted([(concs[i], i) for i in range(len(concs))])
    idx_list = [idx for conc, idx in conc_idx]
    
    # check input titration point    
    input_index = idx_list[0] # input data column index
    assert concs[input_index] == 0 # conc is zero at input

    # read files and start logistic regression
    rsq_list = []
    total_count, fail_count = 0, 0
    
    f = Helper_Py3.open_any(out_fname + '_4PL.gtab.gz', "wt")
    
    for fname in fnames:
        print("Processing %s" % (fname.rsplit('/', 1)[-1]),
              file=sys.stderr)
        First = True
        for line in Helper_Py3.open_any(fname, "rt"):
            line = line.strip()
            if not line:
                continue
            cols = line.split()
            
            if First:
                _, col_st, col_ed, _ = Helper_Py3.read_gtab_header(cols)

                assert col_ed - col_st == len(tnums)  # data col len == titration num

                s = cols[:col_st]
                s += ['Model',
                      'Method',
                      'Top',
                      'Rate',
                      'C-half',
                      'Bottom',
                      'R-squared']
                
                print('\t'.join(s), end='\n', file=f) # write header
                First = False
                continue

            # non target chromosome
            chr_name = cols[0]
            if chr_list and chr_name not in chr_list:
                continue

            # get molecule numbers
            nums = []
            for i in range(col_st, col_ed):
                num = int(cols[i])
                nums.append(num)

            input_num = nums[input_index] # input molecule number

            if input_num <=0: # skip if input has no molecules
                continue

            fracs = [float(num/input_num) for num in nums] # convert num to fraction
            
            # organize titration data for fitting
            X = np.asarray([concs[idx] for idx in idx_list])
            Y = np.asarray([fracs[idx] for idx in idx_list])

            # set guess and boundary of parameters
            top_guess = max(Y)
            bottom_guess = min(Y)
            chalf_guess = np.mean(X)
            rate_guess = 1.0

            if min_top is None:
                top_min = top_guess * (1-0.01)
            else:
                top_min = min_top

            if max_top is None:
                top_max = top_guess * (1+0.01)
            else:
                top_max = max_top

            if min_bottom is None:
                bottom_min = bottom_guess * (1-0.01)
            else:
                bottom_min = min_bottom

            if max_bottom is None:
                bottom_max = bottom_guess * (1+0.01)
            else:
                bottom_max = max_bottom

            if min_chalf is None:
                chalf_min = min(X)
            else:
                chalf_min = min_chalf

            if max_chalf is None:
                chalf_max = max(X)
            else:
                chalf_max = max_chalf

            if min_rate is None:
                rate_min = 0.0
            else:
                rate_min = min_rate

            if max_rate is None:
                rate_max = 100.0
            else:
                rate_max = max_rate
                
            # fitting the data with a logistic function
            if method == 'curve_fit':
                # set initial guess of parameters
                p0 = [top_guess, rate_guess, chalf_guess, bottom_guess]
                # set boundary values of parameters
                bounds = ([top_min, rate_min, chalf_min, bottom_min],
                          [top_max, rate_max, chalf_max, bottom_max])

                try:
                    if model == 'sigmoid':                        
                        p_opt, p_cov = curve_fit(sigmoid_func,
                                                 X,
                                                 Y,
                                                 p0,
                                                 bounds = bounds,
                                                 method='dogbox')
                    elif model == 'hill':
                        p_opt, p_cov = curve_fit(hill_func,
                                                 X,
                                                 Y,
                                                 p0,
                                                 bounds = bounds,
                                                 method='dogbox')
                    success = True
                except (RuntimeError, ValueError):
                    success = False

            elif method == 'evolution':
                # set boundary values of parameters
                bounds = [(top_min, top_max),
                          (rate_min, rate_max),
                          (chalf_min, chalf_max),
                          (bottom_min, bottom_max)]

                if model == 'sigmoid':
                    result = differential_evolution(obj_func,
                                                    args=(sigmoid_func, X, Y),
                                                    bounds=bounds,
                                                    seed=3)
                elif model == 'hill':
                    result = differential_evolution(obj_func,
                                                    args=(hill_func, X, Y),
                                                    bounds=bounds,
                                                    seed=3)

                p_opt = result.x
                success = result.success

            # check fitting quality
            if success:
                if model == 'sigmoid':
                    residuals = np.asarray(Y)- sigmoid_func(X, *p_opt)
                elif model == 'hill':
                    residuals = np.asarray(Y)- hill_func(X, *p_opt)

                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((np.asarray(Y)-np.mean(Y))**2)
                r_squared = 1 - (ss_res / ss_tot)

                # too poor fitting
                if r_squared < min_rsq:
                    success = False
                    pass

            if success:
                top, rate, chalf, bottom = p_opt
                rsq_list.append(r_squared)

            else:
                top, rate, chalf, bottom = 'NA', 'NA', 'NA', 'NA'
                r_squared = 'NA'
                fail_count +=1

            if graph_option:
            #if graph_option and not success:
            #if graph_option and success and rate>100:
            #if graph_option and success and r_squared < 0.7:

                print (top, rate, chalf, bottom)
                print (r_squared)
                #print (CP80)
                
                fig = plt.figure()
                plt.plot(X, Y, '.', markersize=10, alpha=0.2)

                if success:
                    X_pred = np.linspace(min(X), max(X), 1000)
                    if model == 'sigmoid':
                        Y_pred = sigmoid_func(X_pred, *p_opt)
                    elif model == 'hill':
                        Y_pred = hill_func(X_pred, *p_opt)
                    plt.plot(X_pred, Y_pred, 'k-', alpha=0.2)
                    plt.axhline(y=bottom + (top-bottom)*0.5, linestyle='--', color='b')
                    plt.axvline(x=chalf, linestyle='--', color='r')

                plt.xlabel("Concentration")
                plt.ylabel("Soluble fractin")
                plt.title("%s" % (r_squared))
                plt.show()
                plt.close()

            s = cols[:col_st]
            s += [model,
                  method,
                  str(top),
                  str(rate),
                  str(chalf),
                  str(bottom),
                  str(r_squared)]
        
            print('\t'.join(s), end='\n', file=f)
            total_count +=1

    f.close()

    # summarize output
    fail_pct = 0.0 if total_count <= 0 else 100 * float(fail_count) / total_count
    median_rsq = float("nan") if not rsq_list else np.median(rsq_list)
    print("fitting failure %d/%d (%.2f %%)"
          % (fail_count,
             total_count,
             fail_pct),
          file=sys.stderr)
    print("Median R-squared of fitting %.2f" % (median_rsq), file=sys.stderr)
    print("Done", file=sys.stderr)
    

if __name__ == '__main__':
    str2bool = Helper_Py3.str2bool

    parser = ArgumentParser(description='fitting condense-seq data with four parameter logistic function')
    parser.add_argument(metavar='-f',
                        dest="fnames",
                        type=str,
                        nargs='+',
                        help='concatenated num.gtab file list')
    parser.add_argument('-t',
                        dest='tfname',
                        type=str,
                        help='titration filename')
    parser.add_argument('--tnum',
                        dest="tnums",
                        type=int,
                        nargs='+',
                        help='titration number of each columns in num data')
    parser.add_argument('--min_tnum',
                        dest="min_tnum",
                        type=int,
                        default=3,
                        help='minimum titration data number for fitting')
    parser.add_argument('--model',
                        dest="model",
                        type=str,
                        default='sigmoid',
                        help='logistic model for fitting data (sigmoid or hill)')
    parser.add_argument('--method',
                        dest="method",
                        type=str,
                        default='evolution',
                        help='logistic regression method (curve_fit or evolution)')    
    parser.add_argument('--min_rsq',
                        dest="min_rsq",
                        type=float,
                        default=0.5,
                        help='minimum R-squared value for fitting quality')
    parser.add_argument('--min_top',
                        dest="min_top",
                        type=float,
                        help='lower bound of Top parameter in 4PL model')
    parser.add_argument('--max_top',
                        dest="max_top",
                        type=float,
                        help='upper bound of Top parameter in 4PL model')
    parser.add_argument('--min_bottom',
                        dest="min_bottom",
                        type=float,
                        help='lower bound of Bottom parameter in 4PL model')
    parser.add_argument('--max_bottom',
                        dest="max_bottom",
                        type=float,
                        help='upper bound of Bottom parameter in 4PL model')
    parser.add_argument('--min_rate',
                        dest="min_rate",
                        type=float,
                        help='lower bound of Rate parameter in 4PL model')
    parser.add_argument('--max_rate',
                        dest="max_rate",
                        type=float,
                        help='upper bound of Rate parameter in 4PL model')
    parser.add_argument('--min_chalf',
                        dest="min_chalf",
                        type=float,
                        help='lower bound of C-half parameter in 4PL model')
    parser.add_argument('--max_chalf',
                        dest="max_chalf",
                        type=float,
                        help='upper bound of C-half parameter in 4PL model')
    parser.add_argument('--chr',
                        dest="chr_list",
                        type=str,
                        nargs='+',
                        help='target chromosome list')
    parser.add_argument('--graph',
                        dest="graph_option",
                        type=str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='plot the graphs fitting the data')
    parser.add_argument('-o',
                        dest='out_fname',
                        default='output',
                        type=str,
                        help='output prefix filename')
    
    args = parser.parse_args()

    # check titration point of input control is included
    if 0 not in args.tnums:
        print("input control data should be included (titration number is zero)",
              file=sys.stderr)
        sys.exit(1)
        
    # check at least 3 titration points for fitting
    if len(args.tnums) < args.min_tnum:
        print("Need more data points for fitting",
              file=sys.stderr)
        sys.exit(1)
    
    # list target chromosomes
    if not args.chr_list:
        chr_list = None
    else:
        chr_list = sorted(args.chr_list, key=Helper_Py3.chr_key)

    if args.graph_option:
        import matplotlib.pyplot as plt

    # set logistic model
    model = args.model.lower()

    # set regression method
    method = args.method.lower()

    logistic_fit (args.fnames,
                  args.tfname,
                  args.tnums,
                  model,
                  method,
                  args.min_rsq,
                  args.min_top,
                  args.max_top,
                  args.min_bottom,
                  args.max_bottom,
                  args.min_rate,
                  args.max_rate,
                  args.min_chalf,
                  args.max_chalf,
                  chr_list,
                  args.graph_option,
                  args.out_fname
                  )

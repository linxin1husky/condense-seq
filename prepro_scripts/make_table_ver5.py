import sys, math, copy
from argparse import ArgumentParser
import Helper_Py3_compat as Helper_Py3

# get AT content of sequence
def AT_content(seq):
    seq = seq.upper()
    output = 0.0
    for nt in seq:
        if nt in "AT":
            output += 1.0
    return output/len(seq)

# get the number of C in the target motif in both strand
def C_motif (seq, motif='CG', both=True):
    seq = seq.upper()
    num = 0
    for i in range(len(seq)-1):
        if seq[i:i+2] == motif:
            num += 1
    if both:
        rev_seq = Helper_Py3.rev_cmp(seq)
        for i in range(len(rev_seq)-1):
            if rev_seq[i:i+2] == motif:
                num +=1
    return num

def binary_search (sortlist, target):
    st, ed = 0, len(sortlist)-1
    while st <= ed:
        mid = (st+ed) // 2
        if sortlist[mid] == target:
            return mid
        elif sortlist[mid] > target:
            ed = mid - 1 
        elif sortlist[mid] < target:
            st = mid + 1
    return st

# simple hash function when binning size/step is constant
class bin_hash:

    def __init__(self,
                 ID_interval,
                 bin_size,
                 bin_step,
                 max_pos):

        self.ID_value = {}    
        self.bin_size = bin_size
        self.bin_step = bin_step
        self.max_pos = max_pos

        # map bin idx to bin ID
        self.idx_ID = {}
        self.ID_idx = {}
        for ID in ID_interval:
            st, ed = ID_interval[ID]
            assert st % self.bin_step == 0
            assert ed == st + self.bin_size
            idx = st // self.bin_step
            assert idx not in self.idx_ID
            self.idx_ID[idx] = ID
            self.ID_idx[ID] = idx
            
        print("hash function is built", file=sys.stderr)
        
    def find(self, pos):
        find_IDs = []
        idx = pos // self.bin_step
        st = self.bin_step*idx
        ed = st + self.bin_size
        while pos >= st and pos < ed:
            ID = self.idx_ID.get(idx)
            if ID is not None:
                find_IDs.append(ID)
            idx -= 1
            if idx < 0:
                break
            st = self.bin_step*idx
            ed = st + self.bin_size
        return find_IDs

    def insert (self, pos, value):
        find_IDs = self.find(pos)
        for ID in find_IDs:
            if ID not in self.ID_value:
                self.ID_value[ID] = 0.0
            self.ID_value[ID] += value
        return find_IDs
        
    def find_range(self, rst, red):
        find_IDs = []

        idx = rst // self.bin_step
        min_idx = idx
        st = self.bin_step*idx
        ed = st + self.bin_size
        while rst >= st and rst < ed:
            if idx < min_idx:
                min_idx = idx
            idx -= 1
            if idx < 0:
                break
            st = self.bin_step*idx
            ed = st + self.bin_size

        red = min(red, self.max_pos + 1)
        max_idx = (red - 1) // self.bin_step

        for idx in range(min_idx, max_idx+1):
            ID = self.idx_ID.get(idx)
            if ID is not None:
                find_IDs.append(ID)
        return find_IDs

    def insert_range (self, rst, red, value):
        find_IDs = self.find_range(rst, red)
        for ID in find_IDs:
            idx = self.ID_idx[ID]
            st = self.bin_step*idx
            ed = st + self.bin_size
            a, b = max(st, rst), min(ed, red)
            length = b - a
            if ID not in self.ID_value:
                self.ID_value[ID] = 0.0
            self.ID_value[ID] += value*length
        return find_IDs

    def keys (self):
        return self.ID_value.keys()

    def values (self):
        return self.ID_value.values()

    def ID (self, id):
        return self.ID_value[id]
        
    def get (self):
        return self.ID_value

                
# build interval dictionary by using double hashing
class double_hash:
    def __init__(self,
                 ID_interval,
                 domain_size,
                 max_pos):

        self.ID_value = {}
        self.ID_interval = ID_interval

        edID = []
        for ID, interval in ID_interval.items():
            st, ed = interval
            edID.append([ed, ID])
        edID = sorted(edID, key=lambda x: (x[0], x[1]))

        edlist, IDlist = [], []
        for ed, ID in edID:
            edlist.append(ed)
            IDlist.append(ID)
        
        self.domain_size = domain_size
        self.max_pos = max_pos
        self.domain_IDs = {}
        self.domain_num = max_pos // domain_size + 1
        
        for i in range(self.domain_num):
            self.domain_IDs[i] = []
            dst = i * self.domain_size
            ded = min(dst + self.domain_size, max_pos+1)
            idx1 = binary_search(edlist, dst)
            if idx1 == len(edlist):
                continue
            for j in range(idx1, len(edlist)):
                ID = IDlist[j]
                st, ed = self.ID_interval[ID]
                if st < ded:
                    self.domain_IDs[i].append(ID)
                
        print("hash function is built", file=sys.stderr)

    def __str__ (self):
        print("%s\t%s\t%s\t%s" % ("ID", "st", "ed", "value"))
        for ID, value in self.ID_value.items():
            st, ed = self.ID_interval[ID]
            print("%d\t%d\t%d\t%f" % (ID, st, ed, value))
        return

    def find (self, pos):
        find_IDs = []
        
        domain = pos // self.domain_size
        IDs = self.domain_IDs[domain]

        edID = []
        for ID in IDs:
            st, ed  = self.ID_interval[ID]
            edID.append([ed,ID])
        edID = sorted(edID, key=lambda x: (x[0], x[1]))

        edlist, IDlist = [], []
        for ed, ID in edID:
            edlist.append(ed)
            IDlist.append(ID)
            
        idx = binary_search(edlist, pos)
        if idx == len(edlist):
            return find_IDs

        for i in range(idx, len(edlist)):
            ID = IDlist[i]
            st, ed = self.ID_interval[ID]
            if pos >= st and pos < ed:
                find_IDs.append(ID)
            
        return find_IDs

    def insert (self, pos, value):
        find_IDs = self.find(pos)
        for ID in find_IDs:
            if ID not in self.ID_value:
                self.ID_value[ID] = 0.0
            self.ID_value[ID] += value
        return find_IDs

    def find_range (self, rst, red):
        find_IDs = []
        domain1 = rst // self.domain_size
        domain2 = red // self.domain_size

        IDs = set([])
        for i in range(domain1, domain2 + 1):
            IDs |= set(self.domain_IDs[i])
        IDs = list(IDs)

        edID = []
        for ID in IDs:
            st, ed = self.ID_interval[ID]
            edID.append([ed, ID])
        edID = sorted(edID, key=lambda x: (x[0], x[1]))

        edlist, IDlist = [], []
        for ed, ID in edID:
            edlist.append(ed)
            IDlist.append(ID)

        idx1 = binary_search(edlist, rst)
        if idx1 == len(edlist):
            return find_IDs
        
        for j in range(idx1, len(edlist)):
            ID = IDlist[j]
            st, ed = self.ID_interval[ID]
            if st < red:
                find_IDs.append(ID)

        return find_IDs
    
    def insert_range (self, rst, red, value):
        find_IDs = self.find_range(rst, red)
        for ID in find_IDs:
            st, ed = self.ID_interval[ID]
            a, b = max(st, rst), min(ed, red)
            length = b - a
            if ID not in self.ID_value:
                self.ID_value[ID] = 0.0
            self.ID_value[ID] += value*length
        return find_IDs

    def keys (self):
        return self.ID_value.keys()

    def values (self):
        return self.ID_value.values()

    def ID (self, id):
        return self.ID_value[id]
        
    def get (self):
        return self.ID_value

# read gtab file
def read_gtab (fname,
               mode='row',
               chr_choices=None):

    if mode == 'row':
        ID_field_value = {}
    elif mode == 'col':
        field_ID_value = {}
    elif mode == 'both':
        ID_field_value, field_ID_value = {}, {}

    First = True
    data_type = None
    for line in Helper_Py3.open_any(fname, "rt"):
        line = line.strip()

        if not line:
            continue

        cols = line.split('\t')
        if First:
            data_type, col_st, col_ed, _ = Helper_Py3.read_gtab_header(cols)

            field_names = cols[col_st:col_ed]
            field_idxs = range(col_st, col_ed)
                    
            First = False
            continue
        
        if data_type == 'point':
            chr, pos = cols[:col_st]
            ID = (chr, int(pos))
        elif data_type == 'binned':
            chr, start, end = cols[:col_st]
            ID = (chr, int(start), int(end))

        if chr_choices is not None and chr not in chr_choices:
            continue
                    
        for field, k in zip(field_names, field_idxs):
            value = cols[k]
            try:
                value = float(value)
            except ValueError:
                pass

            if mode in ['row', 'both']:
                if ID not in ID_field_value:
                    ID_field_value[ID] = {}
                ID_field_value[ID][field] = value

            if mode in ['col', 'both']:
                if field not in field_ID_value:
                    field_ID_value[field] = {}
                field_ID_value[field][ID] = value
 
    if mode == 'row':
        return ID_field_value, field_names, data_type
    if mode == 'col':
        return field_ID_value, field_names, data_type
    if mode == 'both':
        return ID_field_value, field_ID_value, field_names, data_type


# read fasta and get sequence for each windows
def get_seq_from_FASTA (fasta_fname,
                        chr_ID_win,
                        mode='fullseq'):
    
    chrs = sorted(chr_ID_win.keys(), key=Helper_Py3.chr_key)
    #print chrs
    Find = False
    chr_ID_seq = {}
    for line in Helper_Py3.open_any(fasta_fname, "rt"):
        line = line.strip()
        if line.startswith(">"):
            if Find:
                while len(left) > 0:
                    leftID, win_size, seq = left.pop(0)
                    if mode == 'fullseq':
                        ID_seq[leftID] = seq
                    elif mode == 'ATcontent':
                        ID_seq[leftID] = AT_content(seq)
                assert len(ID_seq) == len(winID)
                assert chr not in chr_ID_seq
                chr_ID_seq[chr] = ID_seq
                #print chr
                #print len(ID_seq), len(winID)

            if len(chr_ID_seq) >= len(chrs):
                Find = False
                continue
                
            chr_name = line.split()[0][1:]
            if chr_name in chrs:
                chr = chr_name
                ID_win = chr_ID_win[chr]
                winID = []
                for ID, win in ID_win.items():
                    wst, wed = win
                    winID.append((wst, wed, ID))
                winID = sorted(winID)
                
                k = 0
                wst, wed, ID = winID[k]
                pt = -1
                ID_seq = {}
                left = []

                Find = True
            else:
                Find = False
            continue

        if Find:                    
            if len(left) == 0 and len(ID_seq) == len(winID):
                continue
            if len(left) == 0 and pt + len(line) < wst:
                pt += len(line)
                continue
            
            for i in range(len(left)):
                leftID, win_size, seq = left.pop(0)
                ed = min(len(line), win_size-len(seq))
                seq += line[:ed]
                if len(seq) == win_size:
                    if mode == 'fullseq':
                        ID_seq[leftID] = seq
                    elif mode == 'ATcontent':
                        ID_seq[leftID] = AT_content(seq)
                else:
                    left.append([leftID, win_size, seq])

            while pt + len(line) >= wst and k < len(winID):
                loc = wst - pt - 1
                win_size = wed - wst
                seq = line[loc:min(loc+win_size, len(line))]
                if len(seq) == win_size:
                    if mode == 'fullseq':
                        ID_seq[ID] = seq
                    elif mode == 'ATcontent':
                        ID_seq[ID] = AT_content(seq)
                else:
                    left.append([ID, win_size, seq])
                k += 1
                if k >= len(winID):
                    break
                wst, wed, ID =  winID[k]

            pt += len(line)

    if Find:
        while len(left) > 0:
            leftID, win_size, seq = left.pop(0)
            if mode == 'fullseq':
                ID_seq[leftID] = seq
            elif mode == 'ATcontent':
                ID_seq[leftID] = AT_content(seq)
        assert len(ID_seq) == len(winID)
        assert chr not in chr_ID_seq
        chr_ID_seq[chr] = ID_seq

    #print sorted(chr_ID_seq.keys(), key=Helper_Py3.chr_key)
    assert len(chrs) == len(chr_ID_seq)
    return chr_ID_seq


# read bisulfite-seq data (ENCODE bethylbed file) and count methylated C for each bin
def bin_BS_file (fnames, chr_intdict):
    chr_ID_C = {} # num of "detected" C in the target motif

    for fname in fnames:
        for line in Helper_Py3.gzopen(fname):
            cols = line.strip().split()
            chrname, st, ed, _, _, strand, _, _, _, reads, frac = cols[:11]

            if chrname not in chr_intdict:
                continue
            
            pos = int(st)
            reads, frac = int(reads), 0.01*float(frac)
            if reads <= 0: # skip "undetected" C in the target motif
                continue

            findIDs = chr_intdict[chrname].insert(pos, frac)
            for ID in findIDs:
                if chrname not in chr_ID_C:
                    chr_ID_C[chrname] = {}
                if ID not in chr_ID_C[chrname]:
                    chr_ID_C[chrname][ID] = 0
                chr_ID_C[chrname][ID] += 1

    chr_ID_meC = {}  # num of "detected" methylated C in the target motif
    for chr in chr_intdict:
        chr_ID_meC[chr] = chr_intdict[chr].get()
        
    return chr_ID_C, chr_ID_meC


# read chip-seq data (ENCODE peak-bed file) and get signal for each bin
def bin_chip_file (fnames, chr_intdict, unit='signal'):
    for fname in fnames:        
        for line in Helper_Py3.gzopen(fname):
            cols = line.strip().split()
            chrname, st, ed, peakname, _, strand, signal, pvalue, qvalue = cols[:9]

            if chrname not in chr_intdict:
                continue
            
            if unit == 'signal':
                score = float(signal)
            elif unit == 'pvalue':
                score = float(pvalue)
            elif unit == 'qvalue':
                score = float(qvalue)
                
            st, ed = int(st), int(ed)
            chr_intdict[chrname].insert_range(st, ed, score)

    chr_ID_value = {}
    for chr in chr_intdict:
        chr_ID_value[chr] = chr_intdict[chr].get()
        
    return chr_ID_value


# read bedgraph file and get value for each bin
def bin_bedgraph_file (fnames, chr_intdict):
    for fname in fnames:
        for line in Helper_Py3.gzopen(fname):
            if not line.startswith('chr'):
                continue
            cols = line.strip().split()
            chrname, st, ed, value = cols

            if chrname not in chr_intdict:
                continue
            
            st, ed = int(st), int(ed)
            value = float(value)
            chr_intdict[chrname].insert_range(st, ed, value)

    chr_ID_value = {}
    for chr in chr_intdict:
        chr_ID_value[chr] = chr_intdict[chr].get()
    
    return chr_ID_value

# read gtab file and get value for each bin
def bin_gtab (fnames, chr_intdict):
    for fname in fnames:
        First = True
        for line in Helper_Py3.gzopen(fname):
            cols = line.split('\t')
            if First:
                if cols[1] == 'Position':
                    data_type = 'point'
                    col_st = 2
                    col_ed = len(cols)
                else:
                    assert cols[1] == 'Start'
                    assert cols[2] == 'End'
                    data_type = 'binned'
                    col_st = 3
                    try:
                        col_ed = cols.index('GCcontent')
                    except ValueError:
                        col_ed = len(cols)

                First = False
                continue

            if data_type == 'point':
                chrname, st, ed = cols[0], int(cols[1]), int(cols[1])
                ed +=1
                value = float(cols[2])
            else:
                chrname, st, ed = cols[0], int(cols[1]), int(cols[2])
                value = float(cols[3])
                
            if chrname not in chr_intdict:
                continue
                
            chr_intdict[chrname].insert_range(st, ed, value)

    chr_ID_value = {}
    for chr in chr_intdict:
        chr_ID_value[chr] = chr_intdict[chr].get()
    
    return chr_ID_value

def make_table(data_fname,
               ref_fname,
               bin_size,
               bin_step,
               bs_fnames,
               chip_fnames,
               bedgraph_fnames,
               gtab_fnames,
               full_seq,
               chr_list,
               genome_size,
               out_fname):    

    # read data reference gtab file
    field_ID_data, field_names, data_type = read_gtab(data_fname,
                                                      mode='col',
                                                      chr_choices=chr_list)

    chr_ID_win = {}
    for ID in field_ID_data[field_names[0]]:
        if len(ID) == 3:
            chr, st, ed = ID
        else:
            chr, pos = ID
            st = pos - bin_size // 2
            ed = pos + bin_size // 2

            if bin_size % 2 != 0:
                ed +=1
        
        if chr not in chr_ID_win:
            chr_ID_win[chr] = {}
        chr_ID_win[chr][ID] = (st, ed)

    chr_list = sorted(chr_ID_win.keys(), key=Helper_Py3.chr_key)
    print("data gtab file reading is done", file=sys.stderr)

    # extract sequence information for each bin
    if full_seq:
        chr_ID_seq = get_seq_from_FASTA(ref_fname,
                                        chr_ID_win,
                                        mode='fullseq')
    else:
        chr_ID_AT = get_seq_from_FASTA(ref_fname,
                                       chr_ID_win,
                                       mode='ATcontent')
        
    print("reference reading is done", file=sys.stderr)
    
    # build interval dictionary for each chromosome
    chr_intdict = {}
    if len(bs_fnames) + len(chip_fnames) + len(bedgraph_fnames) + len(gtab_fnames) > 0:
        print("building interval dictionary", file=sys.stderr)
        for chr in chr_list:
            print("%s" % (chr), file=sys.stderr)
            if bin_step:
                Int_dict = bin_hash(chr_ID_win[chr],
                                    bin_size,
                                    bin_step,
                                    genome_size[chr])
            else:
                domain_size = 10**(int(math.log10(genome_size[chr])) // 2)
                Int_dict = double_hash(chr_ID_win[chr],
                                       domain_size,
                                       genome_size[chr])
            chr_intdict[chr] = Int_dict

    # read bisulfite-seq data (ENCODE bethylbed file) and count methylated C for each bin
    bs_chr_ID_C = {}
    bs_chr_ID_meC = {}
    if bs_fnames:
        for bs in bs_fnames:
            chr_intdict_copy = copy.deepcopy(chr_intdict)
            fnames = bs_fnames[bs]
            chr_ID_C, chr_ID_meC = bin_BS_file(fnames, chr_intdict_copy)
            bs_chr_ID_C[bs] = chr_ID_C
            bs_chr_ID_meC[bs] = chr_ID_meC
            del chr_intdict_copy, chr_ID_C, chr_ID_meC
            print("%s BS reading is done" % (bs), file=sys.stderr)

    # read chip-seq data (ENCODE peak-bed file) and get value for each bin
    chip_chr_ID_value = {}
    if chip_fnames:
        for chip in chip_fnames:
            chr_intdict_copy = copy.deepcopy(chr_intdict)
            fnames = chip_fnames[chip]
            chr_ID_value = bin_chip_file(fnames, chr_intdict_copy)
            chip_chr_ID_value[chip] = chr_ID_value
            del chr_intdict_copy, chr_ID_value
            print("%s Chip reading is done" % (chip), file=sys.stderr)

    # read bedgraph file and get value for each bin
    bg_chr_ID_value = {}
    if bedgraph_fnames:
        for bg in bedgraph_fnames:
            chr_intdict_copy = copy.deepcopy(chr_intdict)
            fnames = bedgraph_fnames[bg]
            chr_ID_value = bin_bedgraph_file(fnames, chr_intdict_copy)
            bg_chr_ID_value[bg] = chr_ID_value
            del chr_intdict_copy, chr_ID_value
            print("%s Bedgraph reading is done" % (bg), file=sys.stderr)

    # read gtab file and get value for each bin
    gtab_chr_ID_value = {}
    if gtab_fnames:
        for gtab in gtab_fnames:
            chr_intdict_copy = copy.deepcopy(chr_intdict)
            fnames = gtab_fnames[gtab]
            chr_ID_value = bin_gtab(fnames, chr_intdict_copy)
            gtab_chr_ID_value[gtab] = chr_ID_value
            del chr_intdict_copy, chr_ID_value
            print("%s gtab reading is done" % (gtab), file=sys.stderr)

    del chr_intdict

    # write annotation file
    print("writing annotation table file", file=sys.stderr)

    f = Helper_Py3.open_any(out_fname + '_table.gtab.gz', "wt")
    if data_type == 'point':
        s = 'Chromosome\tPosition'
    else:
        s = 'Chromosome\tStart\tEnd'
        
    for field_name in field_names:
        s += '\t' + field_name

    if full_seq:
        s += '\t' + 'Sequence'
    s += '\t' + 'ATcontent'

    bs_names = sorted(bs_chr_ID_C.keys())
    for bs in bs_names:
        s += '\t' + 'CNumber(%s)' % (bs)
        s += '\t' + 'meCNumber(%s)' % (bs)

    chip_names = sorted(chip_chr_ID_value.keys())
    for chip in chip_names:
        s += '\t' + chip

    bg_names = sorted(bg_chr_ID_value.keys())
    for bg in bg_names:
        s += '\t' + bg

    gtab_names = sorted(gtab_chr_ID_value.keys())
    for gtab in gtab_names:
        s += '\t' + gtab
    print(s, file=f)


    for chr in chr_list:
        ID_win = chr_ID_win[chr]
        for ID in sorted(ID_win.keys()):
            if len(ID) == 2:
                chr, pos = ID
                s = chr + "\t" + str(pos)
            else:
                chr, st, ed = ID
                s = chr + "\t" + str(st) + "\t" + str(ed)            

            for field_name in field_names:
                data = field_ID_data[field_name][ID]
                try:
                    data = round(data, 5)
                except TypeError:
                    pass
                s += '\t' + str(data)
            
            if full_seq:
                seq = chr_ID_seq[chr][ID]
                s += '\t' + seq
                s += '\t' + str(round(AT_content(seq), 5))
            else:
                s += '\t' + str(round(chr_ID_AT[chr][ID], 5))

            for bs in bs_names:
                C = int(bs_chr_ID_C.get(bs, {}).get(chr, {}).get(ID, 0))
                s += '\t' + str(C)
                meC = round(bs_chr_ID_meC.get(bs, {}).get(chr, {}).get(ID, 0.0), 5)
                s += '\t' + str(meC)

            for chip in chip_names:
                ID_value = chip_chr_ID_value[chip][chr]
                value = round(ID_value.get(ID, 0.0), 5)
                s += '\t' + str(value)

            for bg in bg_names:
                ID_value = bg_chr_ID_value[bg][chr]
                value = round(ID_value.get(ID, 0.0), 5)
                s += '\t' + str(value)

            for gtab in gtab_names:
                ID_value = gtab_chr_ID_value[gtab][chr]
                value = round(ID_value.get(ID, 0.0), 5)
                s += '\t' + str(value)

            print(s, file=f)

    f.close()

    print("Done", file=sys.stderr)
    
if __name__ == '__main__':
    parser = ArgumentParser(description='Make annotation table by putting together gtab file with other epigenetic data')
    parser.add_argument(metavar='-f',
                        dest="data_fname",
                        type=str,
                        help='reference data gtab files')
    parser.add_argument(metavar = '-x',
                        dest='ref_fname',
                        type=str,
                        help='reference genome file')
    parser.add_argument('--binsize',
                        dest="bin_size",
                        type=int,
                        default=171,
                        help='bin window size of cn file in bp (default: 171)')
    parser.add_argument('--binstep',
                        dest="bin_step",
                        type=int,
                        default=0,
                        help='bin step size for cn file in bp (regular binning case)')
    parser.add_argument('--bs',
                        dest="bs_fnames",
                        type=str,
                        nargs='+',
                        action='append',
                        help='bisulfite sequencing files (order: name, ENCODE bedMethyl files)')
    parser.add_argument('--chip',
                        dest="chip_fnames",
                        type=str,
                        nargs='+',
                        action='append',
                        help='chip sequencing files (order: name, ENCODE peak-bed files)')
    parser.add_argument('--bedgraph',
                        dest="bg_fnames",
                        type=str,
                        nargs='+',
                        action='append',
                        help='bedgraph files (order: name, bedgraph files)')
    parser.add_argument('--gtab',
                        dest="gtab_fnames",
                        type=str,
                        nargs='+',
                        action='append',
                        help='gtab files (order: name, gtab files)')
    parser.add_argument('--full-seq',
                        dest="full_seq",
                        type=Helper_Py3.str2bool,
                        nargs='?',
                        const=True,
                        default=False,
                        help='write full sequence information')
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

    # get length for each chromosome
    genome_size = Helper_Py3.genome_sizes(args.ref_fname)

    if not args.chr_list:
        chr_list = sorted(genome_size.keys(), key=Helper_Py3.chr_key)
    else:
        chr_list = sorted(args.chr_list, key=Helper_Py3.chr_key)

    bs_fnames = {}
    if args.bs_fnames:
        for bs_fname in args.bs_fnames:
            if len(bs_fname) <= 1:
                print("need name and file information", file=sys.stderr)
                sys.exit(1)
            bs, fnames = bs_fname[0], bs_fname[1:]
            assert bs not in bs_fnames
            bs_fnames[bs] = fnames

    chip_fnames = {}
    if args.chip_fnames:
        for chip_fname in args.chip_fnames:
            if len(chip_fname) <= 1:
                print("need name and file information", file=sys.stderr)
                sys.exit(1)
            chip, fnames = chip_fname[0], chip_fname[1:]
            assert chip not in chip_fnames
            chip_fnames[chip] = fnames

    bg_fnames = {}
    if args.bg_fnames:
        for bg_fname in args.bg_fnames:
            if len(bg_fname) <= 1:
                print("need name and file information", file=sys.stderr)
                sys.exit(1)
            bg, fnames = bg_fname[0], bg_fname[1:]
            assert bg not in bg_fnames
            bg_fnames[bg] = fnames

    gtab_fnames = {}
    if args.gtab_fnames:
        for gtab_fname in args.gtab_fnames:
            if len(gtab_fname) <= 1:
                print("need name and file information", file=sys.stderr)
                sys.exit(1)
            gtab, fnames = gtab_fname[0], gtab_fname[1:]
            assert gtab not in gtab_fnames
            gtab_fnames[gtab] = fnames

    if args.bin_step:
        bin_step = args.bin_step
    else:
        bin_step = None
                    
    make_table(args.data_fname,
               args.ref_fname,
               args.bin_size,
               bin_step,
               bs_fnames,
               chip_fnames,
               bg_fnames,
               gtab_fnames,
               args.full_seq,
               chr_list,
               genome_size,
               args.out_fname)

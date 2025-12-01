import gzip, sys, argparse

""" This the python3 helper script for any related scripts.
"""

# Replace chr_cmp() with python3_helper.chr_key()
def chr_key(name: str):
    """
    Generate a sorting key for chromosome names.

    Converts chromosome identifiers (e.g., 'chr1', 'chrX', 'chrMT')
    into sortable tuples such that numeric chromosomes are ordered
    naturally before non-numeric ones, with X/Y/MT converted to numeric
    equivalents (X=23, Y=24, M/MT=25).

    Args:
        name (str): Chromosome name string starting with 'chr'.

    Returns:
        tuple: A sorting key of the form (category, value), where
               numeric chromosomes return (0, int), and others return (1, str).

    Example:
        >>> sorted(['chrX', 'chr2', 'chr10'], key=chr_key)
        ['chr2', 'chr10', 'chrX']
    """
    assert name.startswith('chr')
    s = name[3:]
    special = {'X': 23, 'Y': 24, 'M': 25, 'MT': 25}
    if s.isdigit():
        return (0, int(s))
    if s in special:
        return (0, special[s])
    return (1, s)


def gzopen(fname):
    """
    Open a text or gzip-compressed file for reading.

    Automatically detects whether the input filename ends with '.gz'
    and opens it using gzip in binary mode ('rb'), otherwise opens it
    as a regular text file ('r').

    Args:
        fname (str): Path to the input file (plain text or .gz).

    Returns:
        file object: A readable file handle suitable for iteration.

    Example:
        >>> with gzopen('reads.fasta.gz') as f:
        ...     for line in f:
        ...         process(line)
    """
    if fname.endswith('.gz'):
        reading_file = gzip.open(fname, 'rt')
    else:
        reading_file = open(fname, 'r')
    return reading_file


def rev_cmp(seq):
    """
    Return the reverse complement of a DNA sequence.

    Translates each nucleotide (A↔T, C↔G) and reverses the resulting string.
    Ambiguous bases (e.g., 'N') are preserved as is.

    Args:
        seq (str): Input DNA sequence (uppercase letters only).

    Returns:
        str: Reverse-complemented DNA sequence.

    Example:
        >>> rev_cmp("ATCGN")
        'NCGAT'
    """
    dic = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    output = ''
    for nt in seq:
        output += dic[nt]
    return output[::-1]

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

if __name__ == "__main__":
    sys.stderr.write(
        "⚠️  Warning: This script is a helper module and is not intended to be run directly.\n"
        "Please run the main pipeline script instead.\n"
    )
    sys.exit(1)

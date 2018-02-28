import os
import logging
import gzip
try:
    maketrans = ''.maketrans
except AttributeError:
    # fallback for Python 2
    from string import maketrans
from functools import partial
from itertools import *
import numpy as np
from collections import Counter

## Variables ##

# Codon lookup table
codon_dict = {"ATT":"I","ATC":"I","ATA":"I","CTT":"L","CTC":"L","CTA":"L","CTG":"L","TTA":"L","TTG":"L","GTT":"V","GTC":"V","GTA":"V","GTG":"V","TTT":"F","TTC":"F","ATG":"M","TGT":"C","TGC":"C","GCT":"A","GCC":"A","GCA":"A","GCG":"A","GGT":"G","GGC":"G","GGA":"G","GGG":"G","CCT":"P","CCC":"P","CCA":"P","CCG":"P","ACT":"T","ACC":"T","ACA":"T","ACG":"T","TCT":"S","TCC":"S","TCA":"S","TCG":"S","AGT":"S","AGC":"S","TAT":"Y","TAC":"Y","TGG":"W","CAA":"Q","CAG":"Q","AAT":"N","AAC":"N","CAT":"H","CAC":"H","GAA":"E","GAG":"E","GAT":"D","GAC":"D","AAA":"K","AAG":"K","CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGA":"R","AGG":"R","TAA":"-","TAG":"-","TGA":"-"}

# Translation table for reverse complementing
translation_table = maketrans('tcgaTCGA','agctAGCT')

## Functions ##

class KeysNotIdentical(Exception):
    pass

# dicts must have identical keys
def dicts_2_table(labels, *dicts):
    text = '\t' + ('\t').join(labels) + '\n'
    for i in range(1, len(dicts)):
        if sorted(dicts[i-1].keys()) != sorted(dicts[i].keys()):
            raise KeysNotIdentical('All dictionaries passed must have identical keys!')
    keys = dicts[0].keys()
    for key in keys:
        text += ('\t').join([key] + [str(x[key]) for x in dicts]) + '\n'
    return text

# Convert a dictionary to a string
def dict_to_string(dic):
    return ('\n').join(['{} = {}'.format(key, dic[key]) for key in dic])

# Compute statistics from an array
# Requires numpy as np
# Requires Counter
def compute_stats(numbers):
    _count = len(numbers)
    _sum = sum(numbers)
    _min = min(numbers)
    _lq = np.percentile(numbers, 25)
    _mean = float(_sum) / _count
    _median = np.median(numbers)
    _mode = Counter(numbers).most_common(1)
    _uq = np.percentile(numbers, 75)
    _max = max(numbers)
    return {
        'count': _count,
        'sum': _sum,
        'min': _min,
        'lower_quartile': _lq,
        'mean': _mean,
        'median': _median,
        'mode': _mode,
        'upper_quartile': _uq,
        'max': _max
    }

# Reverse complement a sequence
def reverse_complement(sequence, translation_table):
    return str.translate(sequence[::-1], translation_table)

# Complement a sequence
def complement(seq, tt):
    return str.translate(seq, tt)

# Reverse complement function partial
rev_comp = partial(reverse_complement, translation_table = translation_table)
comp = partial(complement, tt = translation_table)

# Translate a sequence
def translate(dna, codon_dict, truncate = True):
    aa = ''
    seqlen = len(dna)
    if seqlen % 3 != 0:
        logging.warning('Sequence of length {} is not divisible by 3'.format(seqlen))
        if truncate:
            dna = dna[:-(seqlen % 3)]
            seqlen = len(dna)
        else:
            return None
    for i in range(0, seqlen, 3):
        try:
            aa += codon_dict[dna[i:i+3]]
        except:
            aa += 'X'
    return aa

# Translate function partial
trans = partial(translate, codon_dict = codon_dict)

# Create a generator for a fasta file
def iterate_fasta(fasta_file):
    with open(fasta_file, 'r') as f:
        iterator = (x[1] for x in groupby(f, lambda x: x[0] == '>'))
        for header in iterator:
            header = header.next().strip()
            seq = ('').join(s.strip() for s in iterator.next())
            yield header, seq

def to_string(_bytes):
    return _bytes.decode('utf-8')

def do_nothing(_string):
    return _string

def parse_fasta(fasta_file, gzipped):
    opener = open
    mode = 'r'
    fxn = do_nothing
    fxn2 = do_nothing
    if gzipped:
        opener = gzip.open
        mode = 'rb'
        fxn = to_string
        fxn2 = chr
    with opener(fasta_file, mode) as f:
        iterator = (x[1] for x in groupby(f, lambda x: fxn2(x[0]) == '>'))
        for header in iterator:
            header = fxn(next(header)).strip()
            seq = ('').join((fxn(s).strip() for s in next(iterator)))
            yield header, seq

def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

def parse_fastq(_file, gzipped):
    opener = open
    mode = 'r'
    fxn = do_nothing
    if gzipped:
        opener = gzip.open
        mode = 'rb'
        fxn = to_string
    with opener(_file, mode) as f:
        for identifier, sequence, extra, quality in grouper(f, 4):
            yield fxn(identifier).strip(), fxn(sequence).strip(), fxn(extra).strip(), fxn(quality).strip()

def get_file_type(_file):
    gzipped = False
    mode = 'r'
    opener = open
    fxn = do_nothing
    if os.path.splitext(_file)[1] == '.gz':
        gzipped = True
        mode = 'rb'
        opener = gzip.open
        fxn = chr
    with opener(_file, mode) as f:
        line = f.readline()
        if fxn(line[0]) == '>':
            return 'fasta', gzipped
        elif fxn(line[0]) == '@':
            return 'fastq', gzipped

def parse_gmap_gff(_file):
    with open(_file, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            cols = line.strip().split('\t')
            attributes = cols[8]
            start = int(cols[3])
            end = int(cols[4])
            feature = cols[2]
            attributes = dict([x.split('=') for x in attributes.split(';')])
            yield cols
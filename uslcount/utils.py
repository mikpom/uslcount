from __future__ import division
from collections import defaultdict
import numpy as np

def jct_freq(jct_pos, trt_len, rl, overhang=3):
    """
    Return how many of N random attempts returned 
    reads spanning splice junction.

    Parameters
    ----------
    jct_pos : list of int positions where splicing occurs. 
              This is 0-base number of the first nucleotide of the next exon.
    trt_len : int, transcript length 
    rl      : int, read length
    overhang  : int, length of overhang a read should have.
    N       : int, number of random permutations
    """
    if rl < overhang:
        raise ValueError('Read length should be larger that overhang')
    if trt_len < rl:
        raise ValueError('Transcript length should be larger than read length')
    if len(jct_pos)==0:
        raise ValueError('At least one junction position should be provided')
    
    # dict with position keys and {0 or 1} values depending whether read
    # starting at this position this will cross a splice junction.  If
    # overhang is less than overhang it is 0.
    _struct = defaultdict(lambda: 0)
    for pos in jct_pos:
        for pos in range(max(0, pos-rl+overhang), max(0, pos-overhang)):
            _struct[pos] = 1

    jctn_pos = sum(_struct.values())
    return jctn_pos / (trt_len - rl)

    # poss = np.random.randint(trt_len-rl+1, size=N)
    # spliced = np.sum([_struct[p] for p in poss])
    # return spliced / N


def jct_freq_pe(jct_pos, trt_len, rl, lsizes, overhang=3, N=1000):
    """
    Return how many of N random attempts returned 
    reads spanning splice junction (paired-end version).

    Parameters
    ----------
    jct_pos : list of int positions where splicing occurs. 
              This is 0-base number of the first nucleotide of the next exon.
    trt_len : int, transcript length 
    rl      : int, read length
    lsizes  : array of ints, possible linker sizes to pick from
    overhang  : int, length of overhang a read should have.
    N       : int, number of random permutations
    """
    
    if rl < overhang:
        raise ValueError('Read length should be larger that overhand')
    if len(jct_pos)==0:
        raise ValueError('At least one junction position should be provided')
    
    # picking random linker sizes
    _lsizes = np.random.choice(lsizes, N)

    # creating random left read and right read positions
    left = np.random.randint(max(1, trt_len-rl+1), size=N)
    right = left+rl+_lsizes

    # picking only those which are within transcript range
    mask = right<=(trt_len-rl)
    tot = np.sum(mask)
    if tot == 0:
        return 0.0
    
    left = left[mask]
    right = right[mask]

    # Initializing masks indicatin whether left or right read
    # will overlap splice junction position
    left_mask = np.zeros(shape=(len(jct_pos), tot), dtype=np.bool)
    right_mask = np.zeros(shape=(len(jct_pos), tot), dtype=np.bool)

    
    i = 0
    for pos in jct_pos:
        left_mask[i, :] = (left>=max(0, pos-rl+overhang))&(left<=max(0, pos-overhang))
        right_mask[i, :] = (right>=max(0, pos-rl+overhang))&(right<=max(0, pos-overhang))
        i+=1

    spliced = np.logical_or(np.any(left_mask, axis=0), np.any(right_mask, axis=0))    
    
    return np.sum(spliced) / tot

def dict2txt(d, fl):
    with open(fl, 'wt') as fh:
        for k, v in d.items():
            fh.write(str(k)+'\t'+str(v)+'\n')

def txt2dict(fl, ktype=str, vtype=str):
    d = {}
    with open(fl, 'rt') as fh:
        for line in fh:
            k,v = line.strip().split('\t')
            k = ktype(k)
            v = vtype(v)
            d[k] = v

    return d

import sys
from bx.intervals.intersection import Interval as bxiv

def inner_complement(ii):
    """
    Function taking list or iterable of Intervals
    and returning intervals representing
    space between them.

    Parameters
    ----------
    ii : list of intervals objects
    """
    srt_feats = sorted(ii, key=lambda f: f.start)
    if len(srt_feats) == 1:
        raise ValueError('At least two intervals should be given')
    _ii = []
    for ind in range(len(srt_feats) - 1):
        # Following Python conventions, i.e. start is 0-based
        # and end is one more than the last position belonging to interval
        ii = bxiv(start=srt_feats[ind].end, end=srt_feats[ind+1].start)
        _ii.append(ii)
    return _ii

# Union of intervals
def union(ii):
    ii = sorted(ii, key=lambda i: i.start)
    ret=[ ii[0] ]
    for x in ii[1:]:
        if ret[-1].end < x.start:
            ret.append(x)
        elif ret[-1].end == x.start:
            ret[-1].end = x.end
        if x.end > ret[-1].end:
            ret[-1].end = x.end
    return ret

# Function returning steps of interval tree values
def steps(tree, start, end):
    ii = tree.find(start, end) # all intervals in (start, end)
    steps = {}
    if ii:
        # Only set of interval ends found in (start, end) have to be checked
        margins = [i.start for i in ii] + [i.end for i in ii] + [start, end]
        pp = sorted(set(filter(lambda m: m>=start and m<=end, margins)))
        for i in range(len(pp)-1):
            vals = {_i.value for _i in tree.find(pp[i], pp[i+1])}
            steps[(pp[i], pp[i+1])] = vals
        if pp[0] > start:
            steps[(0, pp[0])] = set()
        if pp[-1] < end:
            steps[(pp[-1], end)] = set()
    else:
        steps[(start, end)] = set()
    return steps

class ContigousFeat(object):
    """
    Simple interval-like genomic feature
    """
    def __init__(self, id, chrom, start, end, strand):
        self.id = id
        self.chrom = chrom
        self.strand = strand
        self.start = start
        self.end = end

    def length(self):
        return self.end - self.start

    def __repr__(self):
        c = type(self).__name__
        return '<{:s} {:s} {:s}:{:d}-{:d}>'\
            .format(c, self.id, self.chrom, self.start, self.end)

class TranscriptFeat(object):
    def __init__(self, id, chrom, start, end, strand, exons):
        self.id = id
        self.chrom = chrom
        self.strand = strand
        self.start = start
        self.end = end
        self.exons = sorted(exons, key=lambda em: em.start)

    def __repr__(self):
        return '<TranscriptFeat {:s} {:s}:{:d}-{:d}>'\
            .format(self.id, self.chrom, self.start, self.end)

    def get_jct_positions(self):
        """
        Return 0-based locations of splice junctions.
        These are locations of the first base of the next exon.
        """
        pos = []
        if len(self.exons) <= 1:
            return pos
        cur_pos = 0
        for exon in self.exons[:-1]:
            cur_pos += exon.length()
            pos.append(cur_pos)

        return pos

    def get_introns(self):
        if len(self.exons) > 1:
            introns = inner_complement([bxiv(e.start, e.end) for e in self.exons])
        else:
            introns = []
        return introns

    def cum_length(self):
        return sum([e.length() for e in self.exons])

    def length(self):
        return self.end - self.start

class GeneFeat(object):
    def __init__(self, id, chrom, start, end, strand, transcripts, struct=None):
        self.id = id
        self.chrom = chrom
        self.strand = strand
        self.start = start
        self.end = end
        self.transcripts = transcripts
        self._struct = struct

    @property
    def struct(self):
        return self._struct

    @struct.setter
    def struct(self, value):
        self._struct = value

    def __repr__(self):
        return '<GeneFeat {:s} {:s}:{:d}-{:d}>'\
            .format(self.id, self.chrom, self.start, self.end)

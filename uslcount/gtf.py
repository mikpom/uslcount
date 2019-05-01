from __future__ import print_function
import sys
from collections import defaultdict
import gzip
import re
from .genomic import GeneFeat, TranscriptFeat, ContigousFeat

GTF_HEADER  = ['seqname', 'source', 'feature', 'start', 'end', 'score',
               'strand', 'frame']
R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA     = re.compile(r'\s*,\s*')
R_KEYVALUE  = re.compile(r'(\s+|\s*=\s*)')

class Feature(object):
    def __init__(self, data):
        self.chrom = data['seqname']
        self.source = data['source']
        self.featuretype = data['feature']
        self.start = int(data['start']) - 1
        self.end = int(data['end'])
        self.score = data['score']
        self.strand = data['strand']
        self.frame = data['frame']
        self.attr = data['attr']

def parse_gtf(filename):
    """
    Open an optionally gzipped GTF file and generate a dict for each line.
    """
    fn_open = gzip.open if filename.endswith('.gz') else open

    with fn_open(filename, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            else:
                yield Feature(parse(line))

def parse(line):
    """
    Parse a single GTF line and return a dict.
    """
    result = {}

    fields = line.rstrip().split('\t')

    for i, col in enumerate(GTF_HEADER):
        result[col] = _get_value(fields[i])

    # INFO field consists of "key1=value;key2=value;...".
    infos = [x for x in re.split(R_SEMICOLON, fields[8]) if x.strip()]

    attr = {}
    for i, info in enumerate(infos, 1):
        # It should be key="value".
        try:
            key, _, value = re.split(R_KEYVALUE, info, 1)
        # But sometimes it is just "value".
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Ignore the field if there is no value.
        if value:
            attr[key] = _get_value(value)

    result['attr'] = attr
    return result


def _get_value(value):
    if not value:
        return None

    # Strip double and single quotes.
    value = value.strip('"\'')

    # Return a list if the value has a comma.
    if ',' in value:
        value = re.split(R_COMMA, value)
    # These values are equivalent to None.
    elif value in ['', '.', 'NA']:
        return None

    return value

def read_features_gtf(gtf, ftypes=('gene', 'transcript', 'exon'), quiet=False):
    """
    Read GTF file to dictionary {id -> Feature}

    Parameters
    ----------
    gtf    : str, GTF file to read
    ftypes : tuple of str, feature types to read other types will be ignored
    quite  : bool, if True do not print progress
    """
    def _print_progress(cnt, quiet):
        if (cnt % 100000 == 0) and (not quiet):
            print(str(cnt)+' GFF lines processed')
        
    gtf = parse_gtf(gtf)
    feats = {}
    seen = {}
    cnt = 0
    for feat in gtf:
        cnt += 1
        if not feat.featuretype in ftypes:
            _print_progress(cnt, quiet)
            continue
        id_field = feat.featuretype + '_id'
        try:
            feat_id = feat.attr[id_field]
        except KeyError:
            if feat.featuretype == 'exon':
                feat_id = 'exon'+str(cnt)
        if feat_id in seen:
            seen[feat_id] += 1
            feat_id += '_'+str(seen[feat_id])
        else:
            seen[feat_id] = 1
        feats[feat_id] = feat
        _print_progress(cnt, quiet)
    return feats

def get_rels(feats):
    """
    Get relationships in the form a dictionary 
    {gene --> {trt --> [exon Features]}}

    Parameters
    ----------
    feats : iterable of HTSeq Feature objects
    """

    rels = defaultdict(lambda: defaultdict(list))
    for feat_id in feats:
        feat=feats[feat_id]
        if feat.featuretype == 'exon':
            rels[feat.attr['gene_id']]\
                [feat.attr['transcript_id']].append(feat_id)
    rels = dict(rels)
    for k in rels:
        rels[k] = dict(rels[k])

    return rels

def gtf2features(gtf, quiet=False):
    gtf_feats = read_features_gtf(gtf, quiet=quiet)
    rels = get_rels(gtf_feats)
    genes = {}
    for gid in rels:
        transcripts = []
        for trt_id in rels[gid]:
            exons = []
            for exon_id in rels[gid][trt_id]:
                exon = gtf_feats[exon_id]
                exons.append(ContigousFeat(id=exon_id, chrom=exon.chrom,
                                           start=exon.start, end=exon.end,
                                           strand=exon.strand))
            if not trt_id in gtf_feats: # have to infer trt feature ourselves
                chrom = exons[0].chrom
                strand = exons[0].strand
                start = min([e.start for e in exons])
                end = max([e.end for e in exons])
            else: # reusing trt feature from GTF
                trt = gtf_feats[trt_id]
                chrom = trt.chrom
                strand = trt.strand
                start = trt.start
                end = trt.end
            transcripts.append(TranscriptFeat(id=trt_id, chrom=chrom,
                                              start=start, end=end,
                                              strand=strand, exons=exons))
        if not gid in gtf_feats: # have to infer gene feature ourselves
            chrom = transcripts[0].chrom
            strand = transcripts[0].strand
            start = min([t.start for t in transcripts])
            end = max([t.end for t in transcripts])
        else: # reusing gene feature from GTF
            gn = gtf_feats[gid]
            start = gn.start
            end = gn.end
            chrom = gn.chrom
            strand = gn.strand
        genes[gid] = GeneFeat(id=gid, chrom=chrom,
                              start=start, end=end,
                              strand=strand, transcripts=transcripts)
    return genes
    

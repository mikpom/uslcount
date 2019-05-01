from __future__ import print_function, division
import os
from collections import defaultdict, OrderedDict
import json
from uslcount.gtf import gtf2features
from .genomic import inner_complement, union, steps, \
    GeneFeat, TranscriptFeat, ContigousFeat
from .utils import jct_freq, jct_freq_pe, txt2dict, dict2txt
from bx.intervals.intersection import Interval as bxiv
from bx.intervals.intersection import IntervalTree

def write_features(genes, fh):
    """
    Function exports genomic features to text file.
    Only genes, transcripts and exons are considered.
    """

    def _get_base_attr(feat):
        base_attr = 'id', 'chrom', 'start', 'end', 'strand'
        return list(map(lambda attr: str(getattr(feat, attr)), base_attr))

    for gid, gfeat in genes.items():
        for trt in gfeat.transcripts:
            for exon in trt.exons:
                fields = ['exon']
                fields +=  _get_base_attr(exon)
                fh.write('\t'.join(fields)+'\n')
            fields = ['transcript']
            fields += _get_base_attr(trt)
            fields.append('|'.join([e.id for e in trt.exons]))
            fh.write('\t'.join(fields)+'\n')
        fields = ['gene']
        fields += _get_base_attr(gfeat)
        fields.append('|'.join([t.id for t in gfeat.transcripts]))
        fields.append(gfeat.struct)
        fh.write('\t'.join(fields)+'\n')

def read_features(fh):
    """
    Function imports genomic features from text file.
    Only genes, transcripts and exons are considered.
    """

    cur_exons = {}
    cur_transcripts = {}
    feats = OrderedDict()
    for line in fh:
        if line.startswith('exon'):
            id, chrom, start, end, strand = line.strip().split('\t')[1:]
            cur_exons[id] = ContigousFeat(id, chrom, int(start), int(end), strand)
        elif line.startswith('transcript'):
            id, chrom, start, end, strand, ee = line.strip().split('\t')[1:]
            exons = [cur_exons[i] for i in ee.split('|')]
            cur_transcripts[id] = TranscriptFeat(id, chrom, int(start), int(end), strand, exons)
        elif line.startswith('gene'):
            id, chrom, start, end, strand, tt, struct = line.strip().split('\t')[1:]
            transcripts = [cur_transcripts[i] for i in tt.split('|')]
            feats[id] = GeneFeat(id, chrom, int(start), int(end), strand, transcripts, struct)
            cur_exons = {}
            cur_transcripts = {}

    return feats

def genomic_intervals(genes, chroms):
    """
    Function builds bx-python Interval Trees from
    dictionary of genes -> transcripts -> exons
    """
    # Exon intervals by strand
    exons_st = {c:{'+':IntervalTree(), '-':IntervalTree()} for c in chroms}
    # Unstranded exon intervals
    exons_ust = {c:IntervalTree() for c in chroms}
    # Intronic intervals (those where not a single transcript has and exon)
    introns = {c:IntervalTree() for c in chroms}

    for gid in genes:
        gn = genes[gid]
        gn_ii = []
        for trt in gn.transcripts:
            for exon in trt.exons:
                gn_ii.append(bxiv(exon.start, exon.end, gid))

        # Projection (union) of all exons
        _proj = union(gn_ii)

        # Adding exon intervals to the trees
        for uiv in _proj:
            exons_st[gn.chrom][gn.strand].add_interval(uiv)
            exons_ust[gn.chrom].add_interval(uiv)

        # Adding intronic intervals to the trees
        if len(_proj)>1:
            for interiv in inner_complement(_proj):
                introns[gn.chrom].add_interval(bxiv(interiv.start, interiv.end, gid))
    return exons_st, exons_ust, introns

class GenomicData(object):
    def __init__(self, exons_st, exons_ust, introns, intr2gn,
                 strands, len_ust_tot, len_ust_unq, len_st_unq, len_ust_intr,
                 genes, chroms, sj_overhang, data_dir):

        self.exons_st = exons_st
        self.exons_ust = exons_ust
        self.introns = introns
        self.intr2gn = dict(intr2gn)
        self.strands = strands
        self.len_ust_tot = len_ust_tot
        self.len_ust_unq = len_ust_unq
        self.len_st_unq = len_st_unq
        self.len_ust_intr = len_ust_intr
        self.genes = genes
        self.chroms = chroms
        self.sj_overhang = sj_overhang
        self.data_dir = data_dir

    @classmethod
    def load(cls, data_dir, quiet=False):
        if not quiet:
            print('Loading genome...', end='');

        # read simple dictionaries
        strands = txt2dict(os.path.join(data_dir, 'strands'))
        chroms = txt2dict(os.path.join(data_dir, 'chroms'))
        len_ust_tot = txt2dict(os.path.join(data_dir, 'len_ust_tot'), vtype=int)
        len_ust_unq = txt2dict(os.path.join(data_dir, 'len_ust_unq'), vtype=int)
        len_st_unq = txt2dict(os.path.join(data_dir, 'len_st_unq'), vtype=int)
        len_ust_intr = txt2dict(os.path.join(data_dir, 'len_ust_intr'), vtype=int)

        # read introns
        intr2gn = {}
        with open(os.path.join(data_dir, 'introns'), 'rt') as fh:
            intr2gn = {}
            for line in fh:
                chrom, start, end, gg = line.strip().split('\t')
                intr2gn[(chrom, int(start), int(end))] = gg.split('|')

        # read gene features
        with open(os.path.join(data_dir, 'genes'), 'rt') as fh:
            genes = read_features(fh)

        # read building params
        build_param = {}
        with open(os.path.join(data_dir, 'param.txt'), 'rt') as param_file:
            for line in param_file:
                param, value = line.strip().split()
                build_param[param] = value
        sj_overhang = int(build_param['sj_overhang'])
        data_dir = build_param['data_dir']

        # create Interval Trees
        exons_st, exons_ust, introns = genomic_intervals(genes, chroms)
        if not quiet:
            print('done')

        return cls(exons_st=exons_st, exons_ust=exons_ust, introns=introns,
                   intr2gn=intr2gn, strands=strands, len_ust_tot=len_ust_tot,
                   len_ust_unq=len_ust_unq, len_st_unq=len_st_unq,
                   len_ust_intr=len_ust_intr, genes=genes, chroms=chroms,
                   sj_overhang=sj_overhang, data_dir=data_dir)

    def save(self, data_dir):
        if not os.path.isdir(data_dir):
            os.mkdir(data_dir)

        dict2txt(self.strands, os.path.join(data_dir, 'strands'))
        dict2txt(self.chroms, os.path.join(data_dir, 'chroms'))
        dict2txt(self.len_ust_tot, os.path.join(data_dir, 'len_ust_tot'))
        dict2txt(self.len_ust_unq, os.path.join(data_dir, 'len_ust_unq'))
        dict2txt(self.len_st_unq, os.path.join(data_dir, 'len_st_unq'))
        dict2txt(self.len_ust_intr, os.path.join(data_dir, 'len_ust_intr'))

        with open(os.path.join(data_dir, 'introns'), 'wt') as fh:
            for intr, gg in self.intr2gn.items():
                fh.write('\t'.join(map(str, intr)))
                fh.write('\t'+'|'.join(gg)+'\n')

        with open(os.path.join(data_dir, 'genes'), 'wt') as fh:
            write_features(self.genes, fh)

        with open(os.path.join(data_dir, 'param.txt'), 'wt') as param_file:
            param_file.write('sj_overhang\t'+str(self.sj_overhang)+'\n')
            param_file.write('data_dir\t'+data_dir+'\n')

    @classmethod
    def build(cls, fl, out_dir=None, rconf=None, sj_overhang=3, quiet=False):
        if not quiet:
            print('Start read gene models')

        _genes = gtf2features(fl)
        genes = OrderedDict()
        for gid in sorted(_genes, key=lambda g: (_genes[g].chrom, _genes[g].start)):
            genes[gid] = _genes[gid]

        if not quiet:
            print('Done reading gene models')
            print('Reading exonic data')

        chroms = defaultdict(int)
        strands = {}
        intr2gn = defaultdict(set)
        for gm in genes.values():
            strands[gm.id] = gm.strand
            if gm.end > chroms[gm.chrom]:
                chroms[gm.chrom] = gm.end
            for trt in gm.transcripts:
                for intron in trt.get_introns():
                    intr2gn[(gm.chrom, intron.start, intron.end)].add(gm.id)
        chroms = dict(chroms)

        # Create Interval Trees
        exons_st, exons_ust, introns = genomic_intervals(genes, chroms)

        # Compute gene lengths
        len_ust_tot = {g:0 for g in genes}
        len_ust_unq = {g:0 for g in genes}
        len_st_unq = {g:0 for g in genes}
        len_ust_intr = {g:0 for g in genes}

        for chrom in chroms:
            for iv, fs in steps(exons_ust[chrom], 0, chroms[chrom]).items():
                if len(fs) == 1:
                    gn = list(fs)[0]
                    len_ust_unq[gn] += (iv[1]-iv[0])
                    len_ust_tot[gn] += (iv[1]-iv[0])
                elif len(fs) > 1:
                    for gn in fs:
                        len_ust_tot[gn] += (iv[1]-iv[0])

        for chrom in chroms:
            for strand in '+', '-':
                for iv, fs in steps(exons_st[chrom][strand], 0, chroms[chrom]).items():
                    if len(fs) == 1:
                        gn = list(fs)[0]
                        len_st_unq[gn] += (iv[1]-iv[0])

        for chrom in chroms:
            for iv, fs in steps(introns[chrom], 0, chroms[chrom]).items():
                for gn in fs:
                    len_ust_intr[gn] += (iv[1]-iv[0])

        with_intr = {g for gset in intr2gn.values() for g in gset}
        for gid, gfeat in genes.items():
            if gid in with_intr:
                if len_ust_intr[gid] == 0:
                    gfeat.struct = 'multicovered'
                else:
                    gfeat.struct = 'multiexonic'
            else:
                gfeat.struct = 'monoexonic'

        gdat = cls(exons_st=exons_st, exons_ust= exons_ust, introns=introns,
                   intr2gn=intr2gn, strands=strands,
                   len_ust_tot=len_ust_tot, len_ust_unq=len_ust_unq, len_st_unq=len_st_unq,
                   len_ust_intr=len_ust_intr, genes=genes, chroms=chroms,
                   sj_overhang=sj_overhang, data_dir=out_dir)

        if not out_dir is None:
            gdat.save(out_dir)
        if not quiet:
            print('Done!')
        return gdat

    def genes_rgn(self, chrom, start, end):
        """
        Return set of genes in the region (ignoring strand)
        """
        gn_rgn = set()
        for i in self.exons_ust[chrom].find(start, end):
            gn_rgn.add(i.value)
        return gn_rgn

    def jct_freq(self, rl, lsizes=None, sj_overhang=3, N=10000, quiet=False):
        """
        Return dictionary of expected splice junction frequency for all genes
        """
        freq = {}
        cnt = 0
        if not lsizes is None:
            minl = min(lsizes)

        for gid, gn in self.genes.items():
            cnt += 1
            trt_lens = {t.id:t.cum_length() for t in gn.transcripts}
            tot_trt_len = sum(trt_lens.values())
            _freq = {} # per transcript junction observing probabilities
            for trt in gn.transcripts:
                jct_pos = trt.get_jct_positions()
                if len(jct_pos)>0:
                    if lsizes is None: # single-end reads
                        if trt_lens[trt.id]<rl:
                            pr = 0
                        else:
                            pr = jct_freq(jct_pos, trt_lens[trt.id], rl,
                                          overhang=sj_overhang)
                    else: # paired-end reads
                        if (trt_lens[trt.id] < 2*rl + minl) or \
                           (trt_lens[trt.id] < rl):
                            pr = 0
                        else:
                            pr = jct_freq_pe(jct_pos, trt_lens[trt.id], rl,
                                             lsizes, overhang=sj_overhang, N=N)
                    _freq[trt.id] = pr
                else:
                    _freq[trt.id] = 0.0
            freq[gid] = 0
            for trt in gn.transcripts:
                freq[gid] += _freq[trt.id]*trt_lens[trt.id]/tot_trt_len

            if cnt % 1000 == 0 and not quiet:
                print(str(cnt)+ ' genes processed', end='\r')

        return freq

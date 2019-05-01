import unittest
from uslcount.gdata import GenomicData
from uslcount.utils import read_pairs
from pkg_resources import resource_filename as pkg_file
from HTSeq import GenomicInterval
from bx.intervals.intersection import IntervalTree
from bx.intervals.intersection import Interval as bxiv
import numpy as np
from uslcount.main import count_bam

gtf_file = pkg_file('uslcount', 'tests/data/gencode.v28.annotation.chr15.gtf.gz')
tspan3 = 'ENSG00000140391.14'

class test_GenomicData(unittest.TestCase):
    gdat = GenomicData.build(gtf_file, quiet=True)
    def test_data(self):
        gm = self.gdat.gms[tspan3]
        proj = gm.get_projection()
        self.assertEqual(len(proj.intervals), 7)

    #@profile
    def test_genes_region(self):
        for i in range(10000):
            iv = GenomicInterval(chrom='chr15', start=77040229,
                                 end=77067600, strand='.')
            gn_iv = self.gdat.genes_rgn(iv)
            self.assertSetEqual(gn_iv, {tspan3, 'ENSG00000259652.1',
                                        'ENSG00000278991.1'})

            
def union(ivs):            
    ii = sorted(ivs, key=lambda i: i.start)
    ret=[ ii[0] ]
    for x in ii[1:]:
        if ret[-1].end < x.start:
            ret.append(x)
        elif ret[-1].end == x.start:
            ret[-1].end = x.end
        if x.end > ret[-1].end:
            ret[-1].end = x.end
    return ret


def steps(tree, start, end):
    ivs = tree.find(start, end)
    pp = sorted([i.start for i in ivs] + [i.end for i in ivs])
    steps = {}
    for i in range(len(pp)-1):
        vals = {_i.value for _i in tree.find(pp[i], pp[i+1])}
        steps[(pp[i], pp[i+1])] = vals
    if pp[0] > start:
        steps[(0, pp[0])] = set()
    if pp[-1] < end:
        steps[(pp[-1], end)] = set()
    return steps
            
class test_bx_steps(unittest.TestCase):
    #@profile
    def test_strats(self):
        gdir = pkg_file('uslcount', 'tests/data/tmp_strats')
        gdat = GenomicData.load(gdir)

        exon_ivs = {c:IntervalTree() for c in gdat.exons_st.chrom_vectors}
        for gid in gdat.gms:
            gn = gdat.gms[gid]
            ii = []
            for trt in gn.transcripts:
                for exon in trt.exons:
                    ii.append(bxiv(exon.start, exon.end, value=gid))
            for uiv in union(ii):
                exon_ivs[gn.chrom].add_interval(uiv)

        N=10000
        poss = np.random.randint(20000000, size=N)
        for i in range(N):
            pos = poss[i]
            ivs = exon_ivs['chr15'].find(pos, pos+50)
            ii = steps(exon_ivs['chr15'], pos, pos+50)

        for i in range(N):
            pos = poss[i]
            iv = GenomicInterval(chrom='chr15', start=pos, end=pos+50)
            ii = list(gdat.exons_ust[iv].steps())

class test_steps_implementations(unittest.TestCase):
    #@profile
    def test_equality(self):
        gdir = pkg_file('uslcount', 'tests/data/tmp_strats')
        gdat = GenomicData.load(gdir)
        
        exon_ivs = {c:IntervalTree() for c in gdat.exons_st.chrom_vectors}
        for gid in gdat.gms:
            gn = gdat.gms[gid]
            ii = []
            for trt in gn.transcripts:
                for exon in trt.exons:
                    ii.append(bxiv(exon.start, exon.end, value=gid))
            for uiv in union(ii):
                exon_ivs[gn.chrom].add_interval(uiv)

        send = 77050749
        chrom = 'chr15'
        bx_steps = steps(exon_ivs[chrom], 0, send)

        htseq_steps = gdat.exons_ust[GenomicInterval('chr15', 0, send)].steps()
        for iv, fs in htseq_steps:
            self.assertEqual(bx_steps[(iv.start, iv.end)], fs)

class test_bam_pairs(unittest.TestCase):
    def test_pairs(self):
        bam = '/media/hdd/bio/uslcount/mappings/HPA/ERR315425_vs_hg38_gencode28.Aligned.sortedByCoord.out.bam'
        cnt = 0
        for rp in read_pairs(bam):
            cnt += 1
        print(cnt, 'alignments processed')
        
# class test_simple_count(unittest.TestCase):
#     def test_speed(self):
#         gdat = GenomicData.build(gtf_file)
#         bam_f = pkg_file('uslcount', 'tests/data/45_E6_2005_RH_1_B_CELL_NAIVE_vs_hg38.chr15.bam')
#         #knt = count_bam_fast(bam_f, gdat, 'r')
#         knt = count_bam_fast(bam_f, gdat, 'R', filt_kwargs={'ignore_vendor_failed':False})
#         self.assertEqual(knt['ENSG00000140575.12'], 6746)
#         self.assertEqual(knt['ENSG00000269974.1'], 9)
        

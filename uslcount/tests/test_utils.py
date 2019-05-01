import unittest
from uslcount import tests
from bx.intervals.intersection import Interval as bxiv
from uslcount.gtf import parse_gtf
from uslcount.genomic import inner_complement

class test_inner_complement(unittest.TestCase):
    def test_inner_complement_toy_example(self):
        i1 = bxiv(100, 200)
        i2 = bxiv(300, 400)
        ivs = inner_complement([i1, i2])
        self.assertEqual(len(ivs), 1)
        iv = ivs[0]
        self.assertEqual(iv.start, 200)
        self.assertEqual(iv.end, 300)
         
    # TODO merged (normalization) functionality
    # def test_inter_iv_merged(self):
    #     i1 = GenomicInterval(chrom='chr', start=100, end=200, strand='+')
    #     i2 = GenomicInterval(chrom='chr', start=200, end=300, strand='+')
    #     i3 = GenomicInterval(chrom='chr', start=400, end=500, strand='+')
    #     ivs = inter_iv([i1, i2, i3])
    #     print(ivs)
    #     self.assertEqual(len(ivs), 1)
    #     iv = ivs[0]
    #     self.assertEqual(iv.length, 100)
    #     self.assertEqual(iv.start, 300)
    #     self.assertEqual(iv.end, 400)
        
    def test_inner_complement_IL24(self):
        gtf = parse_gtf(tests.gtf_file1)
        exon_feats = filter(lambda f: f.featuretype == 'exon' and \
                             f.attr['transcript_id']=='ENST00000611909.4', gtf)
        exon_feats = list(exon_feats)
        self.assertEqual(len(exon_feats), 5)
        introns = inner_complement(exon_feats)
        self.assertEqual(len(introns), 4)
        intron1 = introns[0]
        self.assertEqual(intron1.end - intron1.start, 115)

    def test_inner_complement_FCMR(self):
        gtf = parse_gtf(tests.gtf_file1)
        # picking one transcript of IL24
        exon_feats = filter(lambda f: f.featuretype == 'exon' and \
                             f.attr['transcript_id']=='ENST00000530505.1', gtf)
        exon_feats = list(exon_feats)
        introns = inner_complement(exon_feats)
        self.assertEqual(len(introns), 2)
        

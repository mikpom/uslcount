import unittest
from uslcount.gtf import read_features_gtf, gtf2features, get_rels
from uslcount import tests
from uslcount.tests import IL24, FCMR

class test_reading(unittest.TestCase):
    def test_read_features_gtf(self):
        feats = read_features_gtf(tests.gtf_file1)
        self.assertEqual(feats[IL24].attr['gene_name'], 'IL24')
        # exon ENSE00002241804.1 should appears in 4 transcript so should be 4 features
        eid = 'ENSE00002241804.1'
        ee = list(filter(lambda fid: fid.startswith(eid), feats))
        self.assertEqual(len(ee), 5)

    def test_get_relations(self):
        feats = read_features_gtf(tests.gtf_file1)
        rels = get_rels(feats)
        self.assertEqual(len(rels[IL24]['ENST00000391929.7']), 7)
        
class test_genomic_model_IL24(unittest.TestCase):
    gmodels = gtf2features(tests.gtf_file1)
    gm = gmodels[IL24]

    def test_reading(self):
        self.assertEqual(len(self.gm.transcripts), 7)

    def test_start(self):
        tstart = min([t.start for t in self.gm.transcripts])
        gstart = self.gm.start
        self.assertEqual(tstart, 206897442)
        self.assertEqual(gstart, 206897442)

    def test_repr(self):
        self.assertEqual(repr(self.gm), '<GeneFeat ENSG00000162892.15 chr1:206897442-206904139>')
        
    def test_introns(self):
        # Getting longest transcript
        trt = max(self.gm.transcripts, key=lambda t: t.length())
        introns = trt.get_introns()
        self.assertEqual(introns[0].start, 206897615)

    def test_jct_positions(self):
        trt = list(filter(lambda t: t.id=='ENST00000367093.3',
                          self.gm.transcripts))[0]
        self.assertListEqual(trt.get_jct_positions(), [172, 318, 517, 580, 655])

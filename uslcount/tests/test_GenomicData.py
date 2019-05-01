import unittest
from pkg_resources import resource_filename as pkg_file
from uslcount.gdata import GenomicData
from uslcount.genomic import steps
from uslcount.tests import IL24, FCMR, gtf_file1, gtf_file2

class test_GenomicData(unittest.TestCase):
    gdat = GenomicData.build(gtf_file1, quiet=True)
    gdat2 = GenomicData.build(gtf_file2, quiet=True)
    def test_data(self):
        self.assertListEqual(list(self.gdat.genes),
                             [IL24, FCMR, 'ENSG00000271680.1', 'ENSG00000226945.1'])
        self.assertEqual(self.gdat.strands[IL24], '+')
        self.assertEqual(self.gdat.intr2gn[('chr1', 206897615, 206897730)], {IL24})
        self.assertEqual(self.gdat.len_ust_unq[IL24], 1979)
        self.assertEqual(self.gdat.len_ust_intr[IL24], 4718)
        self.assertEqual(self.gdat.genes[IL24].struct, 'multiexonic')
        self.assertEqual(self.gdat.genes['ENSG00000226945.1'].struct, 'monoexonic')

    def test_genes_region(self):
        self.assertSetEqual(self.gdat.genes_rgn('chr1', 206903668, 206906348),
                            {IL24, 'ENSG00000271680.1', FCMR})

    def test_steps(self):
        st = steps(self.gdat2.exons_st['chr12']['+'], 120460759, 120461061)
        self.assertSetEqual(st[(120460759, 120461061)], {'ENSG00000257218.5'})

    # TODO more strict test
    def test_jct_freq(self):
        freqs = self.gdat.jct_freq(rl=101, lsizes=[50])
        self.assertGreater(freqs[IL24], 0.1)

class test_reading_gdata(unittest.TestCase):
    gdir = pkg_file('uslcount', 'tests/data/tmp_test_reading_gdata')
    gdat = GenomicData.build(gtf_file1, out_dir=gdir, quiet=True)
    def test_loading(self):
        gdat = GenomicData.load(self.gdir, quiet=True)

    def test_genes_rgn(self):
        gdat = GenomicData.load(self.gdir, quiet=True)
        gg = gdat.genes_rgn('chr1', 206903668, 206906348)
        self.assertSetEqual(gg, {IL24, 'ENSG00000271680.1', FCMR})
        

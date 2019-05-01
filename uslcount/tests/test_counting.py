import unittest
from uslcount.gdata import GenomicData
from uslcount.tests import IL24, FCMR, GATC
from uslcount.main import count_bam
from uslcount import tests

gdat1 = GenomicData.build(tests.gtf_file1, quiet=True)
gdat2 = GenomicData.build(tests.gtf_file2, quiet=True)

class test_counting(unittest.TestCase):
    def test_count_IL24(self):
        cnt = count_bam(tests.bam_file1, gdat1, strand='R', quiet=True)
        self.assertEqual(cnt[cnt['gid']=='ENSG00000226945.1']['st'], 5)
        self.assertEqual(cnt[cnt['gid']==IL24]['st'], 56)
        
    def test_count_IL24_mapq0(self):
        cnt = count_bam(tests.bam_file1, gdat1,
                        strand='R', quiet=True, mapq=0)
        self.assertEqual(cnt[cnt['gid']=='ENSG00000226945.1']['st'], 19)
        
    def test_count_IL24_with_secondary(self):
        cnt = count_bam(tests.bam_file1, gdat1, ignore_secondary=False,
                        strand='R', quiet=True, mapq=0)
        self.assertEqual(cnt[cnt['gid']=='ENSG00000226945.1']['st'], 19+1089)

    def test_count_IL24_with_qcfailed(self):
        cnt = count_bam(tests.bam_file1, gdat1, ignore_qc_failed=False,
                        strand='R', quiet=True, mapq=10)
        self.assertEqual(cnt[cnt['gid']=='ENSG00000226945.1']['st'], 17)
        
    def test_count_dyn(self):
        cnt = count_bam(tests.bam_file2, gdat2, mode='dev', strand='R', quiet=True)
        rec = cnt[cnt['gid']=='ENSG00000219355.2'][0]
        self.assertEqual(rec['st'], 1)
        self.assertEqual(rec['ust'], 9)
        self.assertEqual(rec['ust'], 9)
        for rec in cnt:
            self.assertGreaterEqual(rec['ust'], rec['ust'])

    def test_count_dyn_subsample(self):
        cnt = count_bam(tests.bam_file2_5pct, gdat2, mode='dev', strand='R', quiet=True)
        data = cnt[cnt['gid']==GATC][0]
        self.assertEqual(data['st'], 25)
        self.assertEqual(data['stno'], 11)
        self.assertEqual(data['ustno'], 15)
        self.assertEqual(data['ust'], 59)

    def test_count_dyn_st_pe(self):
        cnt = count_bam(tests.bam_file8, gdat2, paired=True, strand='R', quiet=True)
        data = cnt[cnt['gid']==GATC][0]
        self.assertEqual(data['st'], 44)

class test_counting_sam(unittest.TestCase):
    def test_count_dyn_subsample(self):
        cnt = count_bam(tests.sam_file2_5pct, gdat2, mode='dev', strand='R', quiet=True)
        data = cnt[cnt['gid']==GATC][0]
        self.assertEqual(data['st'], 25)
        self.assertEqual(data['ustno'], 15)
    

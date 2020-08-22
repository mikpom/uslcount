from io import StringIO
import numpy as np
import unittest
from uslcount.gdata import GenomicData
from uslcount.main import count_bam, collect_aln_stats, analyze_bam, \
    write_analyze_res
from uslcount.tests import IL24, FCMR, GATC
from uslcount import tests
from uslcount.bamutils import EmptyBamError
                
class test_aln_analysis(unittest.TestCase):
    gdat2 = GenomicData.build(tests.gtf_file2, quiet=True)
    def test_metrics(self):
        cnt = count_bam(tests.bam_file2, self.gdat2, quiet=True, mode='analyze')
        stats = collect_aln_stats(cnt, self.gdat2, self.gdat2.jct_freq(50))
        data = stats[stats['gid']==GATC][0]
        self.assertAlmostEqual(data['exo2intr'], 4.842997815190262)
        self.assertLess(abs(-4.0543-data['jctF']), 0.2)

class test_full_analysis(unittest.TestCase):
    gdat1 = GenomicData.build(tests.gtf_file1, quiet=True)
    gdat2 = GenomicData.build(tests.gtf_file2, quiet=True)
    gdat3 = GenomicData.build(tests.gtf_file3, quiet=True)
    gdat4 = GenomicData.build(tests.gtf_file4, quiet=True)
    gdat6 = GenomicData.build(tests.gtf_file6, quiet=True)

    def test_empty_bam(self):
        with self.assertRaises(EmptyBamError):
            res = analyze_bam(tests.bam_file9, self.gdat1, strand='R', quiet=True)
    
    def test_IL24_example(self):
        res = analyze_bam(tests.bam_file1, self.gdat1, strand='R', quiet=True)

        self.assertEqual(np.sum(res['ust']>=0), res.shape[0])
        # Extracting record for IL24
        rec = res[res['gid']==IL24][0]
        self.assertEqual(rec['jctn'], 4)
        self.assertGreater(rec['prob'], 0.16)

        out = StringIO()
        write_analyze_res(res, out)
        out.seek(0)
        header = out.readline().strip().split('\t')
        self.assertEqual(header[-1], 'confidence_score')
        for line in out:
            if line.startswith(IL24):
                print(line)
                fields = line.split('\t')
                self.assertIn('4', fields)
                self.assertLessEqual(float(fields[-1]), 0.9)

    def test_IL24_pe(self):
        res = analyze_bam(tests.bam_file5, self.gdat1, quiet=True)
        rec = res[res['gid']==IL24][0]
        self.assertEqual(rec['ustno'], 10)
        self.assertGreater(rec['prob'], 0.2)

        # Check junction counts for FCMR
        rec2 = res[res['gid']==FCMR][0]
        self.assertEqual(rec2['jctn'], 20)

    def test_unsorted_paired(self):
        with self.assertRaises(ValueError):
            analyze_bam(tests.bam_file7, self.gdat1, quiet=True)

    # def test_annotation_mismatch(self):
    #     gdat = GenomicData.build(tests.gtf_file5, quiet=True)
    #     with self.assertWarns(UserWarning):
    #         analyze_bam(tests.bam_file1, gdat, quiet=True)
        
    def test_IL24_pe_spleen(self):
        res = analyze_bam(tests.bam_file6, self.gdat1, quiet=True)
        rec = res[res['gid']==IL24][0]
        self.assertEqual(rec['ustno'], 630)
        self.assertGreater(rec['prob'], 0.19)
        
    def test_IL24_hg19_example(self):
        res = analyze_bam(tests.bam_file3, self.gdat3, quiet=True)
        self.assertEqual(np.sum(res['ust']>=0), res.shape[0])
        rec = res[res['gid']=='IL24'][0]
        #self.assertEqual(rec['opposite_strand']['FAIM3'], 854)
        self.assertGreater(rec['prob'], 0.05)

        rec2 = res[res['gid']=='FAIM3'][0]
        self.assertLess(rec2['prob'], 0.1)

    def test_non_negativity(self):
        res = analyze_bam(tests.bam_file4, self.gdat4, quiet=True)
        rec = res[res['gid']=='ENSG00000169045.17'][0]
        self.assertGreater(res['ust'][0], 0)

    def test_dyn_example(self):
        res = analyze_bam(tests.bam_file2, self.gdat2, quiet=True)
        self.assertEqual(np.sum(res['ust']>=0), res.shape[0])
        rec = res[res['gid']==GATC][0]
        self.assertGreater(rec['prob'], 0.0)

    def test_for_zero_division(self):
        res = analyze_bam(tests.bam_file10, self.gdat6, quiet=True)
        self.assertEqual(res.shape[0], 1)

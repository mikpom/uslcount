import pysam
import unittest
from uslcount import tests
from uslcount.bamutils import read_pairs, get_lsize, get_skipped, lib_param

class test_bam_parsing(unittest.TestCase):
    def test_positions(self):
        infile = pysam.AlignmentFile(tests.bam_file1, 'r')
        for aln in infile:
            if aln.query_name == 'HWI-D00361:324:C981HANXX:2:2302:17748:35698':
                aln_pos = aln.get_blocks()
                self.assertListEqual(aln_pos, [(206902064, 206902072), (206902975, 206903015)])
                skipped = get_skipped(aln)
                self.assertListEqual(skipped, [(206902072,206902975)])

    def test_read_pairs(self):
        _rp = filter(lambda p: p[0] and p[0].query_name=='ERR315425.3152069',
                    read_pairs(tests.bam_file5, quiet=True))
        _rp = list(_rp)
        self.assertEqual(len(_rp), 1)
        rp = _rp[0]
        self.assertEqual(rp[1].reference_end, 206904532)

    def test_read_pairs2(self):
        _rp = filter(lambda p: p[0].query_name=='SRR5424812.15061297',
                    read_pairs(tests.bam_file8, quiet=True))
        _rp = list(_rp)
        self.assertEqual(len(_rp), 1)
        rp = _rp[0]
        self.assertEqual(rp[0].reference_end, 120461323)
                
    def test_pairing_of_secondary_aln(self):
        cnt = 0
        for rp in read_pairs(tests.bam_file6, quiet=True):
            if (rp[0] and rp[0].query_name == 'ERR315338.5627716'):
                cnt += 1
        self.assertEqual(cnt, 3)

    # def test_lib_param(self):
    #     rl, paired = lib_param(tests.bam_file9)
    #     self.assertEqual(rl, 100)

class test_pe_inserts(unittest.TestCase):

    def test_linker_size_bam5(self):
        for aln1, aln2 in read_pairs(tests.bam_file6, quiet=True):
            if aln1 is None:
                continue
            elif aln1.query_name == 'ERR315425.4408087':
                self.assertEqual(get_lsize(aln1, aln2), 259)
            elif aln1.query_name == 'ERR315425.106023':
                self.assertEqual(get_lsize(aln1, aln2), 214)
            elif aln1.query_name == 'ERR315425.8104616':
                self.assertEqual(get_lsize(aln1, aln2), -64)
        
    def test_linker_size_bam6(self):
        for aln1, aln2 in read_pairs(tests.bam_file6, quiet=True):
            if aln1.query_name == 'ERR315338.3417836':
                self.assertEqual(get_lsize(aln1, aln2), -96)
                break
        

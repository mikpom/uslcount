import unittest
from pkg_resources import resource_filename as pkg_file
from uslcount.gdata import GenomicData
from uslcount.main import analyze_bam

gtf_file = pkg_file('uslcount', 'tests/data/gencode.v28.annotation.chr12.gtf.gz')
bam_file_se = '/media/hdd/bio/R24/mappings/45_E6_2005_RH_1_CD4_NAIVE/'\
              +'45_E6_2005_RH_1_CD4_NAIVE_GRCh38_NCBI_analysis_set_gencode_v28.Aligned.sortedByCoord.out.chr12.bam'

bam_file_pe = '/media/hdd/bio/uslcount/HPA/mappings/ERR315335_adrenal_gland/'\
              +'ERR315335_adrenal_gland_GRCh38_NCBI_analysis_set_gencode_v28.Aligned.sortedByCoord.out.chr12.bam'

class test_analysis_se(unittest.TestCase):
    gdat = GenomicData.build(gtf_file)

    def test_analyze(self):
        res = analyze_bam(bam_file_se, self.gdat)
        rec = res[res['gid']=='ENSG00000219355.2'][0]
        print(rec)
        self.assertEqual(rec['ust'], 9)

class test_analysis_pe(unittest.TestCase):
    gdat = GenomicData.build(gtf_file)

    def test_analyze(self):
        res = analyze_bam(bam_file_pe, self.gdat)
        
    

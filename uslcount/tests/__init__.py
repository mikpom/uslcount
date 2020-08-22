import unittest
from pkg_resources import resource_filename as pkg_file

# Nice gene names
IL24 = 'ENSG00000162892.15'
FCMR = 'ENSG00000162894.11'
GATC = 'ENSG00000257218.5'

# GTF file slices
gtf_file1 = pkg_file('uslcount', 'tests/data/IL24_FCMR.gtf')
gtf_file2 = pkg_file('uslcount', 'tests/data/dynein_L_locus.gtf')
gtf_file3 = pkg_file('uslcount', 'tests/data/IL24_FCMR_hg19.gtf')
gtf_file4 = pkg_file('uslcount', 'tests/data/HNRNPH1.gtf')

# containing only chrY annotation for warning tests
gtf_file5 = pkg_file('uslcount', 'tests/data/chrY_annots.gtf')

gtf_file6 = pkg_file('uslcount', 'tests/data/test.zero.gtf')

# samtools view -hb source.bam chr1:206896367-206924281 > IL24_FCMR.bam
# source bam is 51_C7_2007_RH_1_CD4_STIM
bam_file1 = pkg_file('uslcount', 'tests/data/bam/IL24_FCMR.bam')

# samtools view -hb source.bam chr12:120436262-120500073 > dynein_L_locus.bam
# source bam is is CD4_NAIVE RNA-seq data aligned versus
# GRCh38_NCBI_analysis_set with Gencode v28 annotations
bam_file2 = pkg_file('uslcount', 'tests/data//bam/dynein_L_locus.bam')
# 5% subsample of previous file 
bam_file2_5pct = pkg_file('uslcount', 'tests/data/bam/dynein_L_locus.5pct.bam')
sam_file2_5pct = pkg_file('uslcount', 'tests/data/bam/dynein_L_locus.5pct.sam')


# samtools view -bh source.bam chr1:207069771-207096773 > IL24_FCMR_hg19.bam
# source bam is CD4 activated  cells RNA-seq data aligned versus hg19 genome with UCSC annotations
bam_file3 = pkg_file('uslcount', 'tests/data/bam/IL24_FCMR_hg19.bam')

# samtools view -bh source.bam chr5:179614177-179634784 > HNRNPH1.bam
# source bam is Th17 RNA-seq data aligned versus
# GRCh38_NCBI_analysis_set with Gencode v28 annotations
bam_file4 = pkg_file('uslcount', 'tests/data/bam/HNRNPH1.bam')

# samtools view -hb souce.bam chr1:206895688-206924622 > IL24_FCMR_pe.bam
# source bam is ERR315425 aligned versus
# GRCh38_NCBI_analysis_set with Gencode v28 annotations
bam_file5 = pkg_file('uslcount', 'tests/data/bam/IL24_FCMR_pe.bam')

# samtools view -hb souce.bam chr1:206895688-206924622 > IL24_FCMR_pe_spleen.bam
# source bam is ERR315338 aligned versus
# GRCh38_NCBI_analysis_set with Gencode v28 annotations
bam_file6 = pkg_file('uslcount', 'tests/data/bam/IL24_FCMR_pe_spleen.bam')

# Unsorted bam file derived by aligning ERR315338
bam_file7 = pkg_file('uslcount', 'tests/data/bam/ERR315338.unsorted.bam')

# samtools view -s 2.03 -bh source.bam chr12:120436262-120500073 > dynein_L_st_pe.bam
# source bam is SRR5424812 aligned versus
# GRCh38_NCBI_analysis_set with Gencode v28 annotations
bam_file8 = pkg_file('uslcount', 'tests/data/bam/dynein_L_st_pe.bam')

# empty
bam_file9 = pkg_file('uslcount', 'tests/data/bam/empty.bam')

# empty
bam_file10 = pkg_file('uslcount', 'tests/data/bam/simulated.bam')

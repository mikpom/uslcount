import unittest
from .tests import test_GenomicData, test_counting, test_gene_models, \
    test_jctns, test_unstranded_unalysis, test_utils, test_bamutils

tloader = unittest.TestLoader()
main_tests = unittest.TestSuite()
for test_module in test_GenomicData, test_counting, test_gene_models, \
                   test_jctns, test_unstranded_unalysis, test_utils, \
                   test_bamutils:

    suit = tloader.loadTestsFromModule(test_module)
    main_tests.addTest(suit)

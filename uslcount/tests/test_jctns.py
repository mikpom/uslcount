from __future__ import division
import unittest
import numpy as np
from uslcount.utils import jct_freq, jct_freq_pe

class test_spliced_reads_simulation(unittest.TestCase):
    def test_jct_freq_se(self):
        trt_len = 200
        rl = 50
        intron_pos = 100
        sj_overhang=3
        freq = jct_freq([intron_pos], trt_len, rl, overhang=sj_overhang)
        expected = (rl-2*sj_overhang+1) / (trt_len-rl+1)
        delta = abs(freq - expected)
        self.assertLess(delta, 0.01)

    def test_jct_freq_pe(self):
        trt_len = 300
        rl = 50
        sj_overhang=3
        isizes = np.asarray([50, -25])
        intron_pos = 100
        N=100000

        # Get expected probabilities
        p1 = (rl-2*sj_overhang) / (trt_len-(2*rl+isizes[0])+1)
        _rl = 2*rl+isizes[1]
        p2 = (_rl-2*sj_overhang) / (trt_len-_rl+1)
        p0 = 0.5*p1 + 0.5*p2

        freq1 = jct_freq_pe([intron_pos], trt_len, rl, isizes[[0]],
                                       overhang=sj_overhang, N=N)
        self.assertLess(abs( (freq1-p1)/p1 ), 0.05)
        freq2 = jct_freq_pe([intron_pos], trt_len, rl, isizes[[1]],
                                       overhang=sj_overhang, N=N)
        freq2se = jct_freq([intron_pos], trt_len, 2*rl+isizes[1],
                           overhang=sj_overhang)
        self.assertLess(abs( (freq2-p2)/freq2 ), 0.05)
        self.assertLess(abs( (freq2-freq2se)/freq2se ), 0.05)
        
        freq0 = jct_freq_pe([intron_pos], trt_len, rl, isizes,
                                       overhang=sj_overhang, N=N)
        self.assertLess(abs( (freq0-p0)/freq0 ), 0.05)

    def test_jct_freq_pe_short(self):
        trt_len = 50
        rl = 101
        sj_overhang=3
        isizes = np.asarray([50, -25])
        intron_pos = 25
        N=1000

        freq = jct_freq_pe([intron_pos], trt_len, rl, isizes,
                           overhang=sj_overhang, N=N)
        self.assertEqual(freq, 0.0)

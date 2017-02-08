#!/bin/env python3
# kSpecAnal - A Spectrum Analyser
# v20170208_1247, HanishKVC
#

import rtlsdr

sdr = rtlsdr.RtlSdr()

sdr.sample_rate = 1e6
sdr.center_freq = 88e6
#sdr.freq_correction = 0
sdr.gain = 0
print(dir(sdr))

data = sdr.read_samples(1024)
print(data)


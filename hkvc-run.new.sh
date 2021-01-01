#!/bin/sh
python/kspecanal.py SCAN startFreq 80e6 endFreq 120e6
python/kspecanal.py SCAN startFreq 30000000 endFreq 200000000 samplingRate 2000000 gain 48.0 window hanning
python/kspecanal.py ZEROSPAN centerFreq 91000000 samplingRate 2e6 gain 19.1
python/kspecanal.py zeroSpan centerFreq 30000000

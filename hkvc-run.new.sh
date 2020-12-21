#!/bin/sh
python/kspecanal.py SCAN startFreq 30000000 endFreq 200000000 samplingRate 2000000 gain 48.0 window true cumuMode Max
python/kspecanal.py SCAN startFreq 80e6 endFreq 120e6 samplingRate 2e6 gain 48.0
python/kspecanal.py ZEROSPAN centerFreq 91000000 samplingRate 2000000 gain 48.0
python/kspecanal.py ZEROSPAN centerFreq 30000000 samplingRate 2000000 gain 48.0

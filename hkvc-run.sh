#!/bin/sh
python/kspecanal.py ZERO_SPAN 91000000 2000000 48.0
python/kspecanal.py ZERO_SPAN 30000000 2000000 48.0
python/kspecanal.py SCAN 30000000 200000000 2000000 48.0
python/kspecanal.py FULL_SPAN
python/kspecanal.py SCAN 80e6 120e6 2e6 48.0

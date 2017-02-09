#!/bin/sh

for f in 85000000 91100000; do
	#for g in 0.0 0.9 3.7 16.6 28.0 40.2; do
	for g in 3.7 16.6 28.0 40.2; do
		rtl_sdr -f $f -n 1024000 -g $g test-$f-g$g.bin
	done
done


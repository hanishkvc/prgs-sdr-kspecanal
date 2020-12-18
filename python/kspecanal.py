#!/usr/bin/env python3
# kSpecAnal - A simple spectrum analyser using RtlSdr
# HanishKVC, v20201218IST1054
#


import sys
import numpy as np
import matplotlib.pyplot as plt
import rtlsdr


gNonOverlap = 0.1
gCenterFreq = 91.7e6
gSamplingRate = 2.4e6
gFftSize = 2**16
gFullSize = gFftSize*4
gGain = 7.1


def sdr_setup(sdr, fC, fS, gain):
    '''
    Setup rtlsdr.
    Also skip few samples to avoid any junk, when it is settling.
    '''
    sdr.sample_rate = fS
    sdr.center_freq = fC
    sdr.gain = gain
    samples = sdr.read_samples(16*1024)
    print("SetupSDR: fC[{}] fS[{}] gain[{}]".format(fC, fS, gain))


def sdr_curscan(sdr):
    numLoops = int(gFullSize/(gFftSize*gNonOverlap))
    print("curscan: numLoops[{}] fullSize[{}]".format(numLoops, gFullSize))
    samples = sdr.read_samples(gFullSize)
    fftAll = np.zeros(gFftSize)
    for i in range(numLoops):
        iStart = int(i*gFftSize*gNonOverlap)
        iEnd = iStart + gFftSize
        tSamples = samples[iStart:iEnd]
        if len(tSamples) < gFftSize:
            break
        fftN = 2*abs(np.fft.fft(tSamples))/len(tSamples)
        fftAll = (fftAll + fftN)/2
    fftAll[0] = 0
    return fftAll


def zero_span(sdr, d):
    sdr_setup(sdr, d['centerFreq'], d['samplingRate'], d['gain'])
    freqs = np.fft.fftfreq(gFftSize,1/d['samplingRate']) + d['centerFreq']
    print("ZeroSpan: min[{}] max[{}]".format(min(freqs), max(freqs)))
    while True:
        #print(".")
        fftAll = sdr_curscan(sdr)
        plt.cla()
        plt.plot(freqs, fftAll)
        plt.show(block=False)
        plt.pause(0.001)


def scan_range(sdr, d):
    freqSpan = d['samplingRate']
    curFreq = d['startFreq'] + freqSpan/2
    dataFAll = np.array([])
    freqsAll = np.array([])
    while curFreq < d['endFreq']:
        sdr_setup(sdr, curFreq, d['samplingRate'], d['gain'])
        freqs = np.fft.fftfreq(gFftSize,1/d['samplingRate']) + curFreq
        dataF = sdr_curscan(sdr)
        dataFAll = np.append(dataFAll, dataF)
        freqsAll = np.append(freqsAll, freqs)
        plt.cla()
        plt.plot(freqsAll, dataFAll)
        plt.pause(0.001)
        curFreq += freqSpan


gD = {}
gD['samplingRate'] = gSamplingRate
gD['gain'] = gGain

curMode = sys.argv[1].upper()
if curMode == 'SCAN':
    gD['startFreq'] = float(sys.argv[2])
    gD['endFreq'] = float(sys.argv[3])
    gD['centerFreq'] = gD['startFreq'] + ((gD['endFreq'] - gD['startFreq'])/2)
else:
    gD['centerFreq'] = float(sys.argv[2])
    gD['startFreq'] = gD['centerFreq'] - gD['samplingRate']/2
    gD['endFreq'] = gD['centerFreq'] + gD['samplingRate']/2
print("INFO: startFreq[{}] centerFreq[{}] endFreq[{}], samplingRate[{}]".format(gD['startFreq'], gD['centerFreq'], gD['endFreq'], gD['samplingRate']))

sdr = rtlsdr.RtlSdr()
if curMode == 'SCAN':
    scan_range(sdr, gD)
else:
    zero_span(sdr, gD)
sdr.close()

input("Press any key to quit...")

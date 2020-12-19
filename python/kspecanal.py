#!/usr/bin/env python3
# kSpecAnal - A simple spectrum analyser using RtlSdr
# HanishKVC, v20201218IST1054
#


import sys
import numpy as np
import matplotlib.pyplot as plt
import rtlsdr


gNonOverlap = 0.1
gCenterFreq = 92e6
gSamplingRate = 2.4e6
gFftSize = 2**14
gFullSize = gFftSize*8
gGain = 19.1
gWindow = True
gZeroSpanFftDispProcMode = 'LogNoGain'
gScanRangeFftDispProcMode = 'LogNoGain'
gScanRangeClipProcMode = 'HistLowClip'
gScanRangeClipProcMode = 'Clip2MinAmp'


def data_proc(d, vals, dataProc):
    '''
    Process the passed array of values in different ways.

    HistLowClip: Takes a histogram of passed data and resets all values below the 2nd bin to equal 2nd bin.
    Clip2MinAmp: Clip values below a predefined value.
    Log: Convert to dB.
    LogNoGain: Convert to dB and substract the gain applied.
    '''
    if dataProc == 'HistLowClip':
        hist = np.histogram(vals)
        vals[vals[:]<hist[1][1]] = hist[1][1]
    elif dataProc == 'Clip2MinAmp':
        vals = np.clip(vals, (1/256)*0.33, None)
    elif dataProc == 'Log':
        vals = 10*np.log10(vals)
    elif dataProc == 'LogNoGain':
        vals = 10*np.log10(vals)-d['gain']
    return vals


def fftvals_dispproc(d, vals, fftDispProcMode):
    '''
    Process fft value wrt displaying it.
    '''
    if fftDispProcMode == 'Raw':
        return vals
    if fftDispProcMode.startswith('LogNoGain'):
        if fftDispProcMode == 'LogNoGainHistLowClip':
            vals = data_proc(d, vals, 'HistLowClip')
        valLogs = data_proc(d, vals, 'LogNoGain')
    return valLogs


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
    '''
    Scan the currently set freq band (upto max sampling rate supported).
    Inturn normalise the fft on captured samples.

    It does a overlapped sliding over the captured samples, with a per
    fft window of gFftSize.

    gFullSize is the amount of data captured over which overlapped sliding
    is done.

    Based on gWindow, hanning window may be applied to the data before fft.
    The result is compensated wrt windowing related loss of amplitude.

    As IQ data is what is got from the hardware, so both +ve and -ve freqs
    are used to get the embedded signals in the sample data.
    '''
    numLoops = int(gFullSize/(gFftSize*gNonOverlap))
    print("curscan: numLoops[{}] fullSize[{}]".format(numLoops, gFullSize))
    samples = sdr.read_samples(gFullSize)
    fftAll = np.zeros(gFftSize)
    if gWindow:
        win = np.hanning(gFftSize)
    else:
        win = np.ones(gFftSize)
    winAdj = len(win)/np.sum(win)
    for i in range(numLoops):
        iStart = int(i*gFftSize*gNonOverlap)
        iEnd = iStart + gFftSize
        tSamples = samples[iStart:iEnd]
        if len(tSamples) < gFftSize:
            break
        fftN = winAdj*2*abs(np.fft.fft(tSamples*win))/len(tSamples)
        fftAll = (fftAll + fftN)/2
    #fftAll[0] = 0
    return fftAll


def zero_span(sdr, d):
    '''
    Repeatadly keep scanning a specified freq band, which is configured
    by default to be the max sampling rate supported by the hardware.
    '''
    sdr_setup(sdr, d['centerFreq'], d['samplingRate'], d['gain'])
    freqs = np.fft.fftfreq(gFftSize,1/d['samplingRate']) + d['centerFreq']
    freqs = np.fft.fftshift(freqs)
    print("ZeroSpan: min[{}] max[{}]".format(min(freqs), max(freqs)))
    while True:
        fftCur = sdr_curscan(sdr)
        fftCur = np.fft.fftshift(fftCur)
        fftPr = fftvals_dispproc(d, fftCur, gZeroSpanFftDispProcMode)
        plt.cla()
        plt.plot(freqs, fftPr)
        plt.pause(0.001)


def _scan_range(sdr, d, freqsAll, fftAll):
    '''
    Scan a specified range, this can be larger than the freq band
    that can be sampled/scanned by the hardware in one go, in which
    case it will do multiple scans to cover the full specified range.
    '''
    freqSpan = d['samplingRate']
    curFreq = d['startFreq'] + freqSpan/2
    print("_scanRange: start:{} end:{} samplingRate:{}".format(d['startFreq'], d['endFreq'], d['samplingRate']))
    if type(freqsAll) == type(None):
        totalFreqs = d['endFreq'] - d['startFreq']
        numGroups = (int(totalFreqs/freqSpan) + 1)
        totalEntries = numGroups * gFftSize
        print("_scanRange: totalFreqs:{} numGroups:{} totalEntries:{}".format(totalFreqs, numGroups, totalEntries))
        fftAll = np.zeros(totalEntries)
        freqsAll = np.fft.fftshift(np.fft.fftfreq(totalEntries, 1/(numGroups*freqSpan)) + d['startFreq'] + (numGroups*freqSpan)/2)
    i = 0
    while curFreq < d['endFreq']:
        iStart = i*gFftSize
        iEnd = iStart+gFftSize
        sdr_setup(sdr, curFreq, d['samplingRate'], d['gain'])
        freqs = np.fft.fftfreq(gFftSize,1/d['samplingRate']) + curFreq
        freqs = np.fft.fftshift(freqs)
        freqsAll[iStart:iEnd] = freqs
        fftCur = sdr_curscan(sdr)
        fftCur = data_proc(d, fftCur, gScanRangeClipProcMode)
        fftCur = np.fft.fftshift(fftCur)
        fftAll[iStart:iEnd] = fftCur
        fftPr = fftvals_dispproc(d, np.copy(fftAll), gScanRangeFftDispProcMode)
        fftPr[np.isinf(fftPr)] = 0
        plt.cla()
        plt.plot(freqsAll, fftPr)
        plt.pause(0.001)
        curFreq += freqSpan
        i += 1
    return freqsAll, fftAll


def scan_range(sdr, d):
    freqs = None
    ffts = None
    while True:
        freqs, ffts = _scan_range(sdr, d, freqs, ffts)


gD = {}
gD['samplingRate'] = gSamplingRate
gD['gain'] = gGain

if len(sys.argv) < 2:
    curMode = 'ZERO_SPAN'
else:
    curMode = sys.argv[1].upper()
if curMode == 'SCAN':
    gD['startFreq'] = float(sys.argv[2])
    gD['endFreq'] = float(sys.argv[3])
    gD['centerFreq'] = gD['startFreq'] + ((gD['endFreq'] - gD['startFreq'])/2)
else:
    if len(sys.argv) < 3:
        gD['centerFreq'] = gCenterFreq
    else:
        gD['centerFreq'] = float(sys.argv[2])
    gD['startFreq'] = gD['centerFreq'] - gD['samplingRate']/2
    gD['endFreq'] = gD['centerFreq'] + gD['samplingRate']/2
print("INFO: startFreq[{}] centerFreq[{}] endFreq[{}], samplingRate[{}]".format(gD['startFreq'], gD['centerFreq'], gD['endFreq'], gD['samplingRate']))

plt.show(block=False)
sdr = rtlsdr.RtlSdr()
if curMode == 'SCAN':
    scan_range(sdr, gD)
else:
    zero_span(sdr, gD)
sdr.close()

input("Press any key to quit...")

#!/usr/bin/env python3
# kSpecAnal - A simple spectrum analyser using RtlSdr
# HanishKVC, v20201218IST1054
#


import sys
import numpy as np
import matplotlib.pyplot as plt
import rtlsdr


gbPltHeatMap = True
gbPltLevels = True

gNonOverlap = 0.1
gCenterFreq = 92e6
gSamplingRate = 2.4e6
gFftSize = 2**14
gFft2FullMult = 8
gGain = 19.1
gWindow = False
gMinAmp4Clip = (1/256)*0.33
gZeroSpanFftDispProcMode = 'LogNoGain'
gScanRangeFftDispProcMode = 'LogNoGain'
gScanRangeClipProcMode = 'HistLowClip'
gScanRangeClipProcMode = 'Clip2MinAmp'
gCumuMode = 'AVG'


PLTFIG_LEVELS = "Levels"
PLTFIG_HEATMAP = "Heatmap"



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
        vals = np.clip(vals, d['minAmp4Clip'], None)
    elif dataProc == 'Log':
        vals = 10*np.log10(vals)
    elif dataProc == 'LogNoGain':
        vals = 10*np.log10(vals)-d['gain']
    return vals


def data_cumu(d, mode, curVals, cStart, cEnd, newVals, nStart, nEnd):
    '''
    Cumulate data from different buffers using one of possible logics

    Copy: copy values in newVals buffer into curVals buffer at specified offsets.
    Avg: average the values between curVals and newVals buffer and store into curVals.
    Max: copy the larger value between curVals and newVals into curVals.
    '''
    if mode == 'COPY':
        curVals[cStart:cEnd] = newVals[nStart:nEnd]
    elif mode == 'AVG':
        curVals[cStart:cEnd] += newVals[nStart:nEnd]
        curVals[cStart:cEnd] /= 2
    elif mode == 'MAX':
        curVals[cStart:cEnd] = np.max([newVals[nStart:nEnd], curVals[cStart:cEnd]], axis=0)
    else:
        msg = "ERROR: Unknown cumuMode [{}], Quiting...".format(mode)
        prg_quit(d, msg)
    return curVals


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


def sdr_curscan(sdr, d):
    '''
    Scan the currently set freq band (upto max sampling rate supported).
    Inturn normalise the fft on captured samples.

    It does a overlapped sliding over the captured samples, with a per
    fft window of d['fftSize'].

    d['fullSize'] is the amount of data captured over which overlapped sliding
    is done.

    Based on gWindow, hanning window may be applied to the data before fft.
    The result is compensated wrt windowing related loss of amplitude.

    As IQ data is what is got from the hardware, so both +ve and -ve freqs
    are used to get the embedded signals in the sample data.
    '''
    numLoops = int(d['fullSize']/(d['fftSize']*d['nonOverlap']))
    #print("curscan: numLoops[{}] fullSize[{}]".format(numLoops, d['fullSize']))
    samples = sdr.read_samples(d['fullSize'])
    fftAll = np.zeros(d['fftSize'])
    if gWindow:
        win = np.hanning(d['fftSize'])
    else:
        win = np.ones(d['fftSize'])
    winAdj = len(win)/np.sum(win)
    for i in range(numLoops):
        iStart = int(i*d['fftSize']*d['nonOverlap'])
        iEnd = iStart + d['fftSize']
        tSamples = samples[iStart:iEnd]
        if len(tSamples) < d['fftSize']:
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
    freqs = np.fft.fftfreq(d['fftSize'],1/d['samplingRate']) + d['centerFreq']
    freqs = np.fft.fftshift(freqs)
    print("ZeroSpan: min[{}] max[{}]".format(min(freqs), max(freqs)))
    if gbPltHeatMap:
        maxHM = 128
        fftHM = np.zeros((maxHM, d['fftSize']))
        indexHM = 0
        plt.figure(PLTFIG_HEATMAP)
        hm = plt.imshow(fftHM, extent=(0,1, 0,1))
    while True:
        print("INFO:ZeroSpan:", indexHM, fftHM[indexHM,100])
        fftCur = sdr_curscan(sdr, d)
        fftCur = np.fft.fftshift(fftCur)
        fftPr = fftvals_dispproc(d, fftCur, gZeroSpanFftDispProcMode)
        if gbPltHeatMap:
            plt.figure(PLTFIG_HEATMAP)
            fftHM[indexHM,:] = fftPr
            hm.set_data(fftHM)
            hm.autoscale()
            plt.pause(0.001)
            indexHM = (indexHM + 1) % maxHM
        if gbPltLevels:
            plt.figure(PLTFIG_LEVELS)
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
        totalEntries = numGroups * d['fftSize']
        print("_scanRange: totalFreqs:{} numGroups:{} totalEntries:{}".format(totalFreqs, numGroups, totalEntries))
        fftAll = np.zeros(totalEntries)
        freqsAll = np.fft.fftshift(np.fft.fftfreq(totalEntries, 1/(numGroups*freqSpan)) + d['startFreq'] + (numGroups*freqSpan)/2)
    i = 0
    while curFreq < d['endFreq']:
        iStart = i*d['fftSize']
        iEnd = iStart+d['fftSize']
        sdr_setup(sdr, curFreq, d['samplingRate'], d['gain'])
        freqs = np.fft.fftfreq(d['fftSize'],1/d['samplingRate']) + curFreq
        freqs = np.fft.fftshift(freqs)
        freqsAll[iStart:iEnd] = freqs
        fftCur = sdr_curscan(sdr, d)
        fftCur = data_proc(d, fftCur, gScanRangeClipProcMode)
        fftCur = np.fft.fftshift(fftCur)
        fftAll = data_cumu(d, d['cumuMode'], fftAll, iStart, iEnd, fftCur, 0, len(fftCur))
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
        print_info(d)
        freqs, ffts = _scan_range(sdr, d, freqs, ffts)


def handle_args(d):
    d['prgMode'] = 'ZEROSPAN'
    d['samplingRate'] = gSamplingRate
    d['gain'] = gGain
    d['centerFreq'] = gCenterFreq
    d['fftSize'] = gFftSize
    d['nonOverlap'] = gNonOverlap
    d['window'] = gWindow
    d['minAmp4Clip'] = gMinAmp4Clip
    d['cumuMode'] = gCumuMode
    iArg = 1
    while iArg < len(sys.argv):
        curArg = sys.argv[iArg].upper()
        if (curArg == 'ZEROSPAN') or (curArg == 'SCAN'):
            d['prgMode'] = curArg
        elif (curArg == 'CENTERFREQ'):
            iArg += 1
            d['centerFreq'] = float(sys.argv[iArg])
        elif (curArg == 'STARTFREQ'):
            iArg += 1
            d['startFreq'] = float(sys.argv[iArg])
        elif (curArg == 'ENDFREQ'):
            iArg += 1
            d['endFreq'] = float(sys.argv[iArg])
        elif (curArg == 'SAMPLINGRATE'):
            iArg += 1
            d['samplingRate'] = float(sys.argv[iArg])
        elif (curArg == 'GAIN'):
            iArg += 1
            d['gain'] = float(sys.argv[iArg])
        elif (curArg == 'MINAMP4CLIP'):
            iArg += 1
            d['minAmp4Clip'] = float(sys.argv[iArg])
        elif (curArg == 'NONOVERLAP'):
            iArg += 1
            d['nonOverlap'] = float(sys.argv[iArg])
        elif (curArg == 'FFTSIZE'):
            iArg += 1
            d['fftSize'] = int(sys.argv[iArg])
        elif (curArg == 'CUMUMODE'):
            iArg += 1
            d['cumuMode'] = sys.argv[iArg].upper()
        elif (curArg == 'WINDOW'):
            iArg += 1
            if sys.argv[iArg].upper() == 'TRUE':
                d['window'] = True
            else:
                d['window'] = False
        iArg += 1
    if d['prgMode'] == 'SCAN':
        d['centerFreq'] = d['startFreq'] + ((d['endFreq'] - d['startFreq'])/2)
    else:
        d['startFreq'] = d['centerFreq'] - d['samplingRate']/2
        d['endFreq'] = d['centerFreq'] + d['samplingRate']/2
    d['fullSize'] = d['fftSize'] * gFft2FullMult


def print_info(d):
    print("INFO: startFreq[{}] centerFreq[{}] endFreq[{}]".format(d['startFreq'], d['centerFreq'], d['endFreq']))
    print("INFO: samplingRate[{}], gain[{}]".format(d['samplingRate'], d['gain']))
    print("INFO: fullSize[{}], fftSize[{}], nonOverlap[{}], window[{}]".format(d['fullSize'], d['fftSize'], d['nonOverlap'], d['window']))
    print("INFO: minAmp4Clip[{}], cumuMode[{}]".format(d['minAmp4Clip'], d['cumuMode']))


def prg_quit(d, msg = None):
    if type(msg) != type(None):
        print(msg)
    quit()


def plt_figures():
    plt.ion()
    plt.figure(PLTFIG_LEVELS)
    plt.figure(PLTFIG_HEATMAP)



gD = {}
handle_args(gD)
print_info(gD)
plt_figures()
sdr = rtlsdr.RtlSdr()
if gD['prgMode'] == 'SCAN':
    scan_range(sdr, gD)
else:
    zero_span(sdr, gD)
sdr.close()



input("Press any key to quit...")


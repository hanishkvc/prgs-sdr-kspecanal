#!/usr/bin/env python3
# kSpecAnal - A simple spectrum analyser using RtlSdr
# HanishKVC, v20201218IST1054
#


import sys
import time
import signal
import numpy as np
import matplotlib.pyplot as plt
import rtlsdr


gbPltHeatMap = True
gbPltLevels = True

gNonOverlap = 0.1
gCenterFreq = 92e6
gSamplingRate = 2.4e6
gFftSize = 2**14
gFft2FullMult4Less = 8
gFft2FullMult4More = 2
gGain = 19.1
gWindow = False
gMinAmp4Clip = (1/256)*0.33
gZeroSpanFftDispProcMode = 'LogNoGain'
gScanRangeFftDispProcMode = 'LogNoGain'
gScanRangeClipProcMode = 'HistLowClip'
gScanRangeClipProcMode = 'Clip2MinAmp'
gCumuMode = 'AVG'
gScanRangeNonOverlap = 0.75


PRGMODE_SCAN = 'SCAN'
PRGMODE_ZEROSPAN = 'ZEROSPAN'

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


gSdrReadUnit = 2**18
def sdr_read(sdr, length):
    '''
    Read from the rtlsdr.
    If the given length is larger than gSdrReadUnit, then it is broken down
    into multiple smaller reads, to avoid libusb io errors.

    If the length to be read is larger than gSdrReadUnit, then it requires
    to be a multiple of gSdrReadUnit.
    '''
    if length > gSdrReadUnit:
        loopCnt = length//gSdrReadUnit
        readLength = gSdrReadUnit
    else:
        loopCnt = 1
        readLength = length
    samples = np.zeros(length, dtype=complex)
    for i in range(loopCnt):
        samples[i*readLength:(i+1)*readLength] = sdr.read_samples(readLength)
    return samples


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
    samples = sdr_read(sdr, d['fullSize'])
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

    Display the instanteneous signal levels as well as
    history of signal levels as a heat map.
    '''
    sdr_setup(sdr, d['centerFreq'], d['samplingRate'], d['gain'])
    freqs = np.fft.fftfreq(d['fftSize'],1/d['samplingRate']) + d['centerFreq']
    freqs = np.fft.fftshift(freqs)
    print("ZeroSpan: min[{}] max[{}]".format(min(freqs), max(freqs)))
    if d['bPltHeatMap']:
        maxHM = 128
        fftHM = np.zeros((maxHM, d['fftSize']))
        indexHM = 0
        plt.figure(PLTFIG_HEATMAP)
        hm = plt.imshow(fftHM, extent=(0,1, 0,1))
        plt.xticks([0, 0.5, 1], [d['startFreq'], d['centerFreq'], d['endFreq']])
        plt.xlabel("Freqs")
        plt.ylabel("ScanNum")
    prevTime = time.time()
    while True:
        curTime = time.time()
        print("ZeroSpan:{}".format(curTime-prevTime))
        prevTime = curTime
        fftCur = sdr_curscan(sdr, d)
        fftCur = np.fft.fftshift(fftCur)
        fftPr = fftvals_dispproc(d, fftCur, gZeroSpanFftDispProcMode)
        if d['bPltHeatMap']:
            plt.figure(PLTFIG_HEATMAP)
            fftHM[indexHM,:] = fftPr
            hm.set_data(fftHM)
            hm.autoscale()
            plt.draw()
            indexHM = (indexHM + 1) % maxHM
        if d['bPltLevels']:
            plt.figure(PLTFIG_LEVELS)
            plt.cla()
            plt.plot(freqs, fftPr)
            plt.draw()
        if d['bPltHeatMap'] or d['bPltLevels']:
            plt.pause(0.0001)


def _scan_range(sdr, d, freqsAll, fftAll):
    '''
    Scan a specified range, this can be larger than the freq band
    that can be sampled/scanned by the hardware in one go, in which
    case it will do multiple scans to cover the full specified range.
    '''
    freqSpan = d['samplingRate']
    if (((freqSpan*d['scanRangeNonOverlap'])%1) != 0) or ((d['fftSize']*d['scanRangeNonOverlap'])%1 != 0):
        msg = "ERROR: freqSpan [{}] or fftSize[{}] x scanRangeNonOverlap [{}] is not int".format(freqSpan, d['fftSize'], d['scanRangeNonOverlap'])
        prg_quit(d, msg)
    curFreq = d['startFreq'] + freqSpan/2
    print("_scanRange: start:{} end:{} samplingRate:{}".format(d['startFreq'], d['endFreq'], d['samplingRate']))
    if type(freqsAll) == type(None):
        totalFreqs = d['endFreq'] - d['startFreq']
        numGroups = (int(totalFreqs/freqSpan) + 1)
        totalEntries = numGroups * d['fftSize']
        print("_scanRange: totalFreqs:{} numGroups:{} totalEntries:{}".format(totalFreqs, numGroups, totalEntries))
        fftAll = np.zeros(totalEntries)
        freqsAll = np.fft.fftshift(np.fft.fftfreq(totalEntries, 1/(numGroups*freqSpan)) + d['startFreq'] + (numGroups*freqSpan)/2)
        if d['bPltHeatMap']:
            d['fftCursMax'] = 128
            d['fftCursIndex'] = 0
            d['fftCurs'] = np.zeros([d['fftCursMax'], totalEntries])
    i = 0
    while curFreq < d['endFreq']:
        iStart = int(i*d['fftSize']*d['scanRangeNonOverlap'])
        iEnd = iStart+d['fftSize']
        sdr_setup(sdr, curFreq, d['samplingRate'], d['gain'])
        freqs = np.fft.fftfreq(d['fftSize'],1/d['samplingRate']) + curFreq
        freqs = np.fft.fftshift(freqs)
        freqsAll[iStart:iEnd] = freqs
        fftCur = sdr_curscan(sdr, d)
        fftCur = data_proc(d, fftCur, gScanRangeClipProcMode)
        fftCur = np.fft.fftshift(fftCur)
        if d['bPltHeatMap']:
            d['fftCurs'][d['fftCursIndex'], iStart:iEnd] = fftCur
        fftAll = data_cumu(d, d['cumuMode'], fftAll, iStart, iEnd, fftCur, 0, len(fftCur))
        fftPr = fftvals_dispproc(d, np.copy(fftAll), gScanRangeFftDispProcMode)
        fftPr[np.isinf(fftPr)] = 0
        if d['bPltLevels']:
            plt.figure(PLTFIG_LEVELS)
            plt.cla()
            plt.plot(freqsAll, fftPr)
            plt.pause(0.001)
        curFreq += freqSpan*d['scanRangeNonOverlap']
        i += 1
    return freqsAll, fftAll


def scan_range(sdr, d):
    freqs = None
    ffts = None
    if d['bPltHeatMap']:
        plt.figure(PLTFIG_HEATMAP)
        hm = plt.imshow(np.zeros([3,3]), extent=(0,1, 0,1))
    while True:
        print_info(d)
        freqs, ffts = _scan_range(sdr, d, freqs, ffts)
        if d['bPltHeatMap']:
            plt.figure(PLTFIG_HEATMAP)
            hm.set_data(d['fftCurs'])
            hm.autoscale()
            plt.pause(0.001)
            d['fftCursIndex'] = (d['fftCursIndex'] + 1) % d['fftCursMax']


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
    d['bPltHeatMap'] = gbPltHeatMap
    d['bPltLevels'] = gbPltLevels
    d['scanRangeNonOverlap'] = gScanRangeNonOverlap
    iArg = 1
    while iArg < len(sys.argv):
        curArg = sys.argv[iArg].upper()
        if (curArg == PRGMODE_ZEROSPAN) or (curArg == PRGMODE_SCAN):
            d['prgMode'] = curArg
        elif (curArg == 'QUICKFULLSCAN'):
            d['prgMode'] = PRGMODE_SCAN
            d['startFreq'] = 30e6
            d['endFreq'] = 1.5e9
            d['fftSize'] = 256
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
        elif (curArg == 'SCANRANGENONOVERLAP'):
            iArg += 1
            d['scanRangeNonOverlap'] = float(sys.argv[iArg])
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
        elif (curArg == 'BPLTHEATMAP'):
            iArg += 1
            if sys.argv[iArg].upper() == 'TRUE':
                d['bPltHeatMap'] = True
            else:
                d['bPltHeatMap'] = False
        elif (curArg == 'BPLTLEVELS'):
            iArg += 1
            if sys.argv[iArg].upper() == 'TRUE':
                d['bPltLevels'] = True
            else:
                d['bPltLevels'] = False
        else:
            msg = "ERROR:handle_args: Unknown argument [{}]".format(curArg)
            prg_quit(d, msg)
        iArg += 1
    if d['prgMode'] == PRGMODE_SCAN:
        d['centerFreq'] = d['startFreq'] + ((d['endFreq'] - d['startFreq'])/2)
    else:
        d['startFreq'] = d['centerFreq'] - d['samplingRate']/2
        d['endFreq'] = d['centerFreq'] + d['samplingRate']/2
    if d['fftSize'] < (d['samplingRate']//8):
        d['fullSize'] = d['fftSize'] * gFft2FullMult4Less
    else:
        d['fullSize'] = d['fftSize'] * gFft2FullMult4More
    if (d['fullSize'] > gSdrReadUnit) and ((d['fullSize'] % gSdrReadUnit) != 0):
        prg_quit(d, "ERROR:fullSize[{}] Not multiple of gSdrReadUnit[{}]".format(d['fullSize'], gSdrReadUnit))


def print_info(d):
    print("INFO: startFreq[{}] centerFreq[{}] endFreq[{}]".format(d['startFreq'], d['centerFreq'], d['endFreq']))
    print("INFO: samplingRate[{}], gain[{}]".format(d['samplingRate'], d['gain']))
    print("INFO: fullSize[{}], fftSize[{}], cumuMode[{}], window[{}]".format(d['fullSize'], d['fftSize'], d['cumuMode'], d['window']))
    print("INFO: minAmp4Clip[{}], nonOverlap[{}], scanRangeNonOverlap[{}]".format(d['minAmp4Clip'], d['nonOverlap'], d['scanRangeNonOverlap']))
    print("INFO: prgMode [{}], bPltLevels[{}],  bPltHeatMap[{}]".format(d['prgMode'], d['bPltLevels'], d['bPltHeatMap']))


def prg_quit(d, msg = None):
    if type(msg) != type(None):
        print(msg)
    sys.exit()


def plt_figures(d):
    plt.ion()
    if d['bPltLevels']:
        plt.figure(PLTFIG_LEVELS)
    if d['bPltHeatMap']:
        plt.figure(PLTFIG_HEATMAP)


def handle_sigint(signum, stack):
    prg_quit(gD, "INFO:sigint: quiting on user request...")


def handle_signals(d):
    signal.signal(signal.SIGINT, handle_sigint)



gD = {}
handle_args(gD)
print_info(gD)
handle_signals(gD)
plt_figures(gD)
sdr = rtlsdr.RtlSdr()
if gD['prgMode'] == PRGMODE_SCAN:
    scan_range(sdr, gD)
else:
    zero_span(sdr, gD)
sdr.close()



input("Press any key to quit...")


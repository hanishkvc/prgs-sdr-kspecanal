#!/usr/bin/env python3
# kSpecAnal - A simple spectrum analyser using RtlSdr
# HanishKVC, v20201226IST1854
#


import sys
import time
import signal
import numpy as np
import matplotlib.pyplot as plt
#import rtlsdr
import testfft as rtlsdr



PRGMODE_SCAN = 'SCAN'
PRGMODE_ZEROSPAN = 'ZEROSPAN'
PRGMODE_ALIAS_FMSCAN = 'FMSCAN'
PRGMODE_ALIAS_QUICKFULLSCAN = 'QUICKFULLSCAN'
PLTFIG_LEVELS = "Levels"
PLTFIG_HEATMAP = "Heatmap"
PLTCOMPRESS_MAX = 'MAX'
PLTCOMPRESS_AVG = 'AVG'
PLTCOMPRESS_RAW = 'RAW'
CUMUMODE_MAX = 'MAX'
CUMUMODE_AVG = 'AVG'
CUMUMODE_RAW = 'RAW'



gPrgModeDefault = PRGMODE_ALIAS_FMSCAN
gbPltHeatMap = True
gbPltLevels = True
gCurScanNonOverlap = 0.1
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
gPrgLoopCnt = 8192
gXRes = 2048
gPltCompress = PLTCOMPRESS_AVG



def data_proc(d, vals, dataProc, infTo = None):
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
        if infTo != None:
            vals[np.isinf(vals)] = infTo
    elif dataProc == 'LogNoGain':
        #print("min", np.min(vals), "max", np.max(vals))
        #print("DBUG: Zeros", np.argwhere(vals == 0))
        vals = 10*np.log10(vals)-d['gain']
        #print("DBUG: -Inf", np.argwhere(vals == - np.Inf))
        if infTo != None:
            vals[np.isinf(vals)] = infTo
    return vals


def data_cumu(d, mode, curVals, cStart, cEnd, newVals, nStart, nEnd):
    '''
    Cumulate data from different buffers using one of possible logics

    Raw: copy values in newVals buffer into curVals buffer at specified offsets.
    Avg: average the values between curVals and newVals buffer and store into curVals.
    Max: copy the larger value between curVals and newVals into curVals.
    '''
    if mode == CUMUMODE_RAW:
        curVals[cStart:cEnd] = newVals[nStart:nEnd]
    elif mode == CUMUMODE_AVG:
        curVals[cStart:cEnd] += newVals[nStart:nEnd]
        curVals[cStart:cEnd] /= 2
    elif mode == CUMUMODE_MAX:
        curVals[cStart:cEnd] = np.max([newVals[nStart:nEnd], curVals[cStart:cEnd]], axis=0)
    else:
        msg = "ERROR: Unknown cumuMode [{}], Quiting...".format(mode)
        prg_quit(d, msg)
    return curVals


def fftvals_dispproc(d, vals, fftDispProcMode, infTo=None):
    '''
    Process fft value wrt displaying it.
    '''
    if fftDispProcMode == 'Raw':
        return vals
    if fftDispProcMode.startswith('LogNoGain'):
        if fftDispProcMode == 'LogNoGainHistLowClip':
            vals = data_proc(d, vals, 'HistLowClip')
        valLogs = data_proc(d, vals, 'LogNoGain', infTo)
    return valLogs


def data_plotcompress(d, xData, yData):
    '''
    Reduce the amount of data, while still maintaining any significant values in the data.
    Length and XRes are assumed to belong to the powers of 2 series.
    '''
    if d['pltCompress'] == PLTCOMPRESS_RAW:
        return xData, yData
    xLen = len(xData)
    #xReduce = np.ceil(xLen / d['xRes'])
    xReduce = int(xLen / d['xRes'])
    rows = int(xLen/xReduce)
    cols = xReduce
    xTData = xData.reshape(rows, cols)
    xVals = np.average(xTData, axis=1)
    yTData = yData.reshape(rows, cols)
    if d['pltCompress'] == PLTCOMPRESS_MAX:
        yVals = np.max(yTData, axis=1)
    elif d['pltCompress'] == PLTCOMPRESS_AVG:
        yVals = np.average(yTData, axis=1)
    else:
        prg_quit(d, "ERROR:pltCompress:1D: Unknown mode [{}]".format(d['pltCompress']))
    return xVals, yVals


def data_2d_plotcompress(d, data):
    '''
    Reduce the number of elements in a 2D data set,
    by merging adjacent cols of each row, into a smaller subset of cols.
    '''
    if d['pltCompress'] == PLTCOMPRESS_RAW:
        return data
    yLen,xLen = data.shape
    xReduce = int(xLen/d['xRes'])
    rows = int(xLen/xReduce)
    cols = xReduce
    newData = np.zeros((yLen, rows))
    for y in range(yLen):
        xData = data[y,:]
        xTData = xData.reshape(rows, cols)
        if d['pltCompress'] == PLTCOMPRESS_MAX:
            xVals = np.max(xTData, axis=1)
        elif d['pltCompress'] == PLTCOMPRESS_AVG:
            xVals = np.average(xTData, axis=1)
        else:
            prg_quit(d, "ERROR:pltCompress:2D: Unknown mode [{}]".format(d['pltCompress']))
        newData[y,:] = xVals
    return newData


gPltHighsDelta4Marking = 0.025
gPltHighsNumMarkers = 5
gPltHighsPause = False
def plot_highs(d, freqs, levels):
    d['AxFreqs'].clear()
    d['AxFreqs'].set_xlabel("Freqs[MHz] - HighSigLvl")
    d['AxFreqs'].set_xticks([])
    d['AxFreqs'].set_yticks([])
    freqRange = freqs[-1] - freqs[0]
    delta4Marking = d['pltHighsDelta4Marking']*freqRange
    print("PlotHighs: Freqs {} to {} : delta4Marking {}".format(freqs[0], freqs[-1], delta4Marking))
    ordered = levels.argsort()
    marked = np.array([])
    cntMarked = 0
    for i in np.arange(-1,-len(freqs),-1):
        #print("PlotHighs:MarkedList:", marked)
        curFreq = freqs[ordered[i]]
        curLevel = levels[ordered[i]]
        matched = marked[abs(marked - curFreq) < delta4Marking]
        if len(matched) == 0:
            print("plotHighs:Marked: {}, {}".format(curFreq, curLevel))
            d['AxLevels'].plot(curFreq, curLevel, "o", label=curFreq)
            d['AxFreqs'].text(0.1,1.0-0.1*(cntMarked+1),curFreq/1e6)
            marked = np.append(marked, curFreq)
            cntMarked += 1
            if cntMarked >= d['pltHighsNumMarkers']:
                break
        else:
            #print("plotHighs:Skipped: {}, {}".format(curFreq, curLevel))
            pass
    d['AxLevels'].legend()
    if d['pltHighsPause']:
        input("PltHighsPause: Press any key to continue...")


def sdr_setup(sdr, fC, fS, gain):
    '''
    Setup rtlsdr.
    Also skip few samples to avoid any junk, when it is settling.

    If there is a failure in setting up the rtlsdr, then it closes
    and reopens a new instance of rtlsdr. The same is returned to
    the caller, along with a success or failure boolean.
    '''
    try:
        sdr.sample_rate = fS
        sdr.center_freq = fC
        sdr.gain = gain
        bOk = True
        samples = sdr.read_samples(16*1024)
    except:
        print("WARN:SetupSDR:FAILED: fC[{}] fS[{}] gain[{}]".format(fC, fS, gain))
        sdr.close()
        sdr = rtlsdr.RtlSdr()
        bOk = False
    print("SetupSDR:{}: fC[{}] fS[{}] gain[{}]".format(bOk, fC, fS, gain))
    return sdr, bOk


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


def sdr_curscan(d):
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
    numLoops = int(d['fullSize']/(d['fftSize']*d['curScanNonOverlap']))
    #print("curscan: numLoops[{}] fullSize[{}]".format(numLoops, d['fullSize']))
    samples = sdr_read(d['sdr'], d['fullSize'])
    fftAll = np.zeros(d['fftSize'])
    if gWindow:
        win = np.hanning(d['fftSize'])
    else:
        win = np.ones(d['fftSize'])
    winAdj = len(win)/np.sum(win)
    for i in range(numLoops):
        iStart = int(i*d['fftSize']*d['curScanNonOverlap'])
        iEnd = iStart + d['fftSize']
        tSamples = samples[iStart:iEnd]
        if len(tSamples) < d['fftSize']:
            break
        fftN = winAdj*2*abs(np.fft.fft(tSamples*win))/len(tSamples)
        fftAll = (fftAll + fftN)/2
    #fftAll[0] = 0
    return fftAll


def zero_span(d):
    '''
    Repeatadly keep scanning a specified freq band, which is configured
    by default to be the max sampling rate supported by the hardware.

    Display the instanteneous signal levels as well as
    history of signal levels as a heat map.
    '''
    d['sdr'], bSdrSetup = sdr_setup(d['sdr'], d['centerFreq'], d['samplingRate'], d['gain'])
    freqs = np.fft.fftfreq(d['fftSize'],1/d['samplingRate']) + d['centerFreq']
    freqs = np.fft.fftshift(freqs)
    print("ZeroSpan: min[{}] max[{}]".format(min(freqs), max(freqs)))
    if d['bPltHeatMap']:
        maxHM = 128
        fftHM = np.zeros((maxHM, d['fftSize']))
        indexHM = 0
        hm = d['AxHeatMap'].imshow(fftHM, extent=(0,1, 0,1), aspect='auto')
        d['AxHeatMap'].set_xticks([0, 0.5, 1])
        d['AxHeatMap'].set_xticklabels([d['startFreq'], d['centerFreq'], d['endFreq']])
        d['AxHeatMap'].set_xlabel("Freqs")
        d['AxHeatMap'].set_ylabel("ScanHistory")
    prevTime = time.time()
    for i in range(d['prgLoopCnt']):
        if d['cmd.stop']:
            break
        curTime = time.time()
        print("ZeroSpan:{}:{}".format(i, curTime-prevTime))
        prevTime = curTime
        fftCur = sdr_curscan(d)
        fftCur = np.fft.fftshift(fftCur)
        fftPr = fftvals_dispproc(d, fftCur, gZeroSpanFftDispProcMode)
        if d['bPltHeatMap']:
            fftHM[indexHM,:] = fftPr
            hm.set_data(fftHM)
            hm.autoscale()
            plt.draw()
            indexHM = (indexHM + 1) % maxHM
        if d['bPltLevels']:
            d['AxLevels'].cla()
            xFreqs, yLvls = data_plotcompress(d, freqs, fftPr)
            d['AxLevels'].plot(xFreqs, yLvls)
            plt.draw()
            plot_highs(d, xFreqs, yLvls)
        plt.pause(0.0001)


def _scan_range(d, freqsAll, fftAll):
    '''
    Scan a specified range, this can be larger than the freq band
    that can be sampled/scanned by the hardware in one go, in which
    case it will do multiple scans to cover the full specified range.

    It uses cumuMode to decide how to merge data from across multiple
    scans of the same freq band.
    '''
    freqSpan = d['samplingRate']
    if (((freqSpan*d['scanRangeNonOverlap'])%1) != 0) or ((d['fftSize']*d['scanRangeNonOverlap'])%1 != 0):
        msg = "ERROR: freqSpan [{}] or fftSize[{}] x scanRangeNonOverlap [{}] is not int".format(freqSpan, d['fftSize'], d['scanRangeNonOverlap'])
        prg_quit(d, msg)
    curFreq = d['startFreq'] + freqSpan/2
    startFreq = curFreq - freqSpan/2
    print("_scanRange: start:{} end:{} samplingRate:{}".format(d['startFreq'], d['endFreq'], d['samplingRate']))
    totalFreqs = d['endFreq'] - d['startFreq']
    numGroups = int(totalFreqs/freqSpan)
    totalEntries = numGroups * d['fftSize']
    print("_scanRange: totalFreqs:{} numGroups:{} totalEntries:{}".format(totalFreqs, numGroups, totalEntries))
    if type(freqsAll) == type(None):
        fftAll = np.ones(totalEntries) * d['minAmp4Clip']
        freqsAll = np.fft.fftshift(np.fft.fftfreq(totalEntries, 1/(numGroups*freqSpan)) + d['startFreq'] + (numGroups*freqSpan)/2)
        if d['bPltHeatMap']:
            d['fftCursMax'] = 128
            d['fftCursIndex'] = 0
            d['fftCurs'] = np.ones([d['fftCursMax'], totalEntries]) * d['minAmp4Clip']
    i = 0
    while startFreq < d['endFreq']:
        iStart = int(i*d['fftSize']*d['scanRangeNonOverlap'])
        iEnd = iStart+d['fftSize']
        sStart = 0
        if iEnd > totalEntries:
            sEnd = iEnd - iStart - (iEnd - totalEntries)
        else:
            sEnd = iEnd - iStart
        d['sdr'], bSdrSetup = sdr_setup(d['sdr'], curFreq, d['samplingRate'], d['gain'])
        freqs = np.fft.fftfreq(d['fftSize'],1/d['samplingRate']) + curFreq
        freqs = np.fft.fftshift(freqs)
        #print("_scanRange: iStart {}-{}, iEnd {}-{}, freqsMin {}, freqsMax {}, freqsLen {}".format(iStart, sStart, iEnd, sEnd, np.min(freqs), np.max(freqs), len(freqs)))
        freqsAll[iStart:iEnd] = freqs[sStart:sEnd]
        if bSdrSetup:
            fftCur = sdr_curscan(d)
        else:
            print("WARN:_scanRange: Dummy data for {} to {}".format(startFreq, startFreq+freqSpan))
            fftCur = np.ones(d['fftSize'])
        fftCur = data_proc(d, fftCur, gScanRangeClipProcMode)
        fftCur = np.fft.fftshift(fftCur)
        if d['bPltHeatMap']:
            d['fftCurs'][d['fftCursIndex'], iStart:iEnd] = fftCur[sStart:sEnd]
        fftAll = data_cumu(d, d['cumuMode'], fftAll, iStart, iEnd, fftCur, sStart, sEnd)
        fftPr = fftvals_dispproc(d, np.copy(fftAll), gScanRangeFftDispProcMode, infTo=0)
        if d['bPltLevels']:
            xFreqs, yLvls = data_plotcompress(d, freqsAll, fftPr)
            d['AxLevels'].clear()
            d['AxLevels'].plot(xFreqs, yLvls)
        curFreq += freqSpan*d['scanRangeNonOverlap']
        startFreq = curFreq - freqSpan/2
        plt.pause(0.0001)
        i += 1
    if d['bPltLevels']:
        fftPr = fftvals_dispproc(d, np.copy(fftAll), gScanRangeFftDispProcMode, infTo=0)
        xFreqs, yLvls = data_plotcompress(d, freqsAll, fftPr)
        plot_highs(d, xFreqs, yLvls)
    return freqsAll, fftAll


def scan_range(d):
    freqs = None
    ffts = None
    # Adjust endFreq such that start-end is a multiple of samplingRate, if required
    freqBands = (d['endFreq'] - d['startFreq'])/d['samplingRate']
    if (freqBands % 1) != 0:
        d['orig.EndFreq'] = d['endFreq']
        d['endFreq'] = d['startFreq'] + np.ceil(freqBands)*d['samplingRate']
        print("WARN:scanRange:Adjusting endFreq: orig [{}] adjusted [{}], so that fullRange is Multiple of samplingRate/freqBand [{}]".format(d['orig.EndFreq'], d['endFreq'], d['samplingRate']))
        #input("Press any key to continue...")
    if d['bPltHeatMap']:
        hm = d['AxHeatMap'].imshow(np.zeros([3,3]), extent=(0,1, 0,1), aspect='auto')
        centerFreq = d['startFreq'] + (d['endFreq'] - d['startFreq'])/2
        d['AxHeatMap'].set_xticks([0, 0.5, 1])
        d['AxHeatMap'].set_xticklabels([d['startFreq'], centerFreq, d['endFreq']])
        d['AxHeatMap'].set_xlabel("Freqs")
        d['AxHeatMap'].set_ylabel("ScanHistory")
    prevTime = time.time()
    for i in range(d['prgLoopCnt']):
        if d['cmd.stop']:
            break
        curTime = time.time()
        print("ZeroSpan:{}:{}".format(i, curTime-prevTime))
        prevTime = curTime
        #print_info(d)
        if (i % 2) == 0:
            d['AxLevels'].cla()
        freqs, ffts = _scan_range(d, freqs, ffts)
        if d['bPltHeatMap']:
            #print("DBUG:scanRange: min[{}] max[{}]".format(np.min(d['fftCurs']), np.max(d['fftCurs'])))
            hmData = data_2d_plotcompress(d, d['fftCurs'])
            hm.set_data(hmData)
            hm.autoscale()
            plt.pause(0.001)
            d['fftCursIndex'] = (d['fftCursIndex'] + 1) % d['fftCursMax']


def _arg_boolean(value):
    if value.upper() == "TRUE":
        return True
    else:
        return False


def handle_args(d):
    d['prgMode'] = gPrgModeDefault
    d['samplingRate'] = gSamplingRate
    d['gain'] = gGain
    d['centerFreq'] = gCenterFreq
    d['fftSize'] = gFftSize
    d['curScanNonOverlap'] = gCurScanNonOverlap
    d['window'] = gWindow
    d['minAmp4Clip'] = gMinAmp4Clip
    d['cumuMode'] = gCumuMode
    d['bPltHeatMap'] = gbPltHeatMap
    d['bPltLevels'] = gbPltLevels
    d['scanRangeNonOverlap'] = gScanRangeNonOverlap
    d['prgLoopCnt'] = gPrgLoopCnt
    d['xRes'] = gXRes
    d['pltCompress'] = gPltCompress
    d['pltHighsNumMarkers'] = gPltHighsNumMarkers
    d['pltHighsDelta4Marking'] = gPltHighsDelta4Marking
    d['pltHighsPause'] = gPltHighsPause
    iArg = 1
    while iArg < len(sys.argv):
        curArg = sys.argv[iArg].upper()
        if (curArg == PRGMODE_ZEROSPAN) or (curArg == PRGMODE_SCAN) or (curArg == PRGMODE_ALIAS_FMSCAN) or (curArg == PRGMODE_ALIAS_QUICKFULLSCAN):
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
        elif (curArg == 'CURSCANNONOVERLAP'):
            iArg += 1
            d['curScanNonOverlap'] = float(sys.argv[iArg])
        elif (curArg == 'SCANRANGENONOVERLAP'):
            iArg += 1
            d['scanRangeNonOverlap'] = float(sys.argv[iArg])
        elif (curArg == 'FFTSIZE'):
            iArg += 1
            d['fftSize'] = int(sys.argv[iArg])
        elif (curArg == 'XRES'):
            iArg += 1
            d['xRes'] = int(sys.argv[iArg])
        elif (curArg == 'CUMUMODE'):
            iArg += 1
            d['cumuMode'] = sys.argv[iArg].upper()
        elif (curArg == 'PLTCOMPRESS'):
            iArg += 1
            d['pltCompress'] = sys.argv[iArg].upper()
        elif (curArg == 'WINDOW'):
            iArg += 1
            d['window'] = _arg_boolean(sys.argv[iArg])
        elif (curArg == 'BPLTHEATMAP'):
            iArg += 1
            d['bPltHeatMap'] = _arg_boolean(sys.argv[iArg])
        elif (curArg == 'BPLTLEVELS'):
            iArg += 1
            d['bPltLevels'] = _arg_boolean(sys.argv[iArg])
        elif (curArg == 'PRGLOOPCNT'):
            iArg += 1
            d['prgLoopCnt'] = int(sys.argv[iArg])
        elif (curArg == 'PLTHIGHSNUMMARKERS'):
            iArg += 1
            d['pltHighsNumMarkers'] = int(sys.argv[iArg])
        elif (curArg == 'PLTHIGHSDELTA4MARKING'):
            iArg += 1
            d['pltHighsDelta4Marking'] = float(sys.argv[iArg])
        elif (curArg == 'PLTHIGHSPAUSE'):
            iArg += 1
            d['pltHighsPause'] = _arg_boolean(sys.argv[iArg])
        else:
            msg = "ERROR:handle_args: Unknown argument [{}]".format(curArg)
            prg_quit(d, msg)
        iArg += 1
    if (d['prgMode'] == PRGMODE_ALIAS_FMSCAN):
        d['prgMode'] = PRGMODE_SCAN
        d['startFreq'] = 88e6
        d['endFreq'] = 108e6
    elif (d['prgMode'] == PRGMODE_ALIAS_QUICKFULLSCAN):
        d['prgMode'] = PRGMODE_SCAN
        d['startFreq'] = 30e6
        d['endFreq'] = 1.5e9
        #d['fftSize'] = 256
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
    print("INFO: fullSize[{}], fftSize[{}], cumuMode[{}], window[{}], xRes[{}], pltCompress[{}]".format(d['fullSize'], d['fftSize'], d['cumuMode'], d['window'], d['xRes'], d['pltCompress']))
    print("INFO: minAmp4Clip[{}], curScanNonOverlap[{}], scanRangeNonOverlap[{}]".format(d['minAmp4Clip'], d['curScanNonOverlap'], d['scanRangeNonOverlap']))
    print("INFO: prgMode [{}], prgLoopCnt[{}], bPltLevels[{}],  bPltHeatMap[{}]".format(d['prgMode'], d['prgLoopCnt'], d['bPltLevels'], d['bPltHeatMap']))
    print("INFO: pltHighsNumMarkers[{}], pltHighsDelta4Marking[{}], pltHighsPause[{}]".format(d['pltHighsNumMarkers'], d['pltHighsDelta4Marking'], d['pltHighsPause']))


def prg_quit(d, msg = None, tryExit=True):
    if type(msg) != type(None):
        print(msg)
    d['cmd.stop'] = True
    if tryExit:
        sys.exit()


def event_pause(event):
    if gD['pltHighsPause']:
        gD['pltHighsPause'] = False
        gD['BtnPause'].label.set_text("Pause[ ]")
    else:
        gD['pltHighsPause'] = True
        gD['BtnPause'].label.set_text("Pause[x]")


def event_levels(event):
    if gD['bPltLevels']:
        gD['bPltLevels'] = False
        gD['BtnLevels'].label.set_text("Levels[ ]")
    else:
        gD['bPltLevels'] = True
        gD['BtnLevels'].label.set_text("Levels[x]")


def event_heatmap(event):
    if gD['bPltHeatMap']:
        gD['bPltHeatMap'] = False
        gD['BtnHeatMap'].label.set_text("HeatMap[ ]")
    else:
        gD['bPltHeatMap'] = True
        gD['BtnHeatMap'].label.set_text("HeatMap[x]")


def event_quit(event):
    gD['BtnQuit'].label.set_text("QuitWait")
    prg_quit(gD, "INFO:QuitClick: Quiting on user request...", False)


def plt_figures(d):
    plt.ion()
    # 4,5 => [[2,4],[2,1]], [[2,5]]
    # 4,5 => [[4,4],[4,1]]
    # 4,5 => [[4,5]]
    f = plt.figure("kSpecAnal", figsize=(12, 8), constrained_layout=True)
    gs = f.add_gridspec(nrows=8, ncols=5)
    d['AxLevels'] = f.add_subplot(gs[:4,:4])
    d['AxFreqs'] = f.add_subplot(gs[:4,4])
    d['AxFreqs'].set_xlabel("Freqs - HighSigLvl")
    d['AxHeatMap'] = f.add_subplot(gs[4:8,:4])
    d['AxBtnPause'] = f.add_subplot(gs[4,4])
    d['AxBtnLevels'] = f.add_subplot(gs[5,4])
    d['AxBtnHeatMap'] = f.add_subplot(gs[6,4])
    d['AxBtnQuit'] = f.add_subplot(gs[7,4])
    d['AxFreqs'].set_xticks([])
    d['AxFreqs'].set_yticks([])
    d['BtnPause'] = plt.Button(d['AxBtnPause'], "Pause")
    d['BtnPause'].on_clicked(event_pause)
    d['BtnLevels'] = plt.Button(d['AxBtnLevels'], "Levels")
    d['BtnLevels'].on_clicked(event_levels)
    d['BtnHeatMap'] = plt.Button(d['AxBtnHeatMap'], "HeatMap")
    d['BtnHeatMap'].on_clicked(event_heatmap)
    d['BtnQuit'] = plt.Button(d['AxBtnQuit'], "Quit")
    d['BtnQuit'].on_clicked(event_quit)


def handle_sigint(signum, stack):
    prg_quit(gD, "INFO:sigint: quiting on user request...")


def handle_signals(d):
    signal.signal(signal.SIGINT, handle_sigint)



gD = {}
gD['cmd.stop'] = False
handle_args(gD)
print_info(gD)
handle_signals(gD)
plt_figures(gD)
gD['sdr'] = rtlsdr.RtlSdr()
if gD['prgMode'] == PRGMODE_SCAN:
    scan_range(gD)
else:
    zero_span(gD)
gD['sdr'].close()



gD['BtnQuit'].label.set_text("QuitPress")
input("Press any key to quit...")


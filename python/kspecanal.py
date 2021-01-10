#!/usr/bin/env python3
# kSpecAnal - A simple spectrum analyser using RtlSdr
# HanishKVC, v20201226IST1854
#


import sys
import time
import signal
import numpy as np
import matplotlib.pyplot as plt
import pickle
import rtlsdr
#import testfft as rtlsdr



PRGMODE_SCAN = 'SCAN'
PRGMODE_ZEROSPAN = 'ZEROSPAN'
PRGMODE_ZEROSPANSAVE = 'ZEROSPANSAVE'
PRGMODE_ZEROSPANPLAY = 'ZEROSPANPLAY'
PRGMODE_ALIAS_FMSCAN = 'FMSCAN'
PRGMODE_ALIAS_QUICKFULLSCAN = 'QUICKFULLSCAN'
PLTFIG_LEVELS = "Levels"
PLTFIG_HEATMAP = "Heatmap"
PLTCOMPRESS_MAX = 'MAX'
PLTCOMPRESS_MIN = 'MIN'
PLTCOMPRESS_AVG = 'AVG'
PLTCOMPRESS_RAW = 'RAW'
PLTCOMPRESS_CONV = 'CONV'
CUMUMODE_MAX = 'MAX'
CUMUMODE_MIN = 'MIN'
CUMUMODE_AVG = 'AVG'
CUMUMODE_RAW = 'RAW'
WINDOW_HANNING = 'WIN.HANNING'
WINDOW_ONES = 'WIN.ONES'
WINDOW_HAMMING = 'WIN.HAMMING'
WINDOW_KAISER = 'WIN.KAISER'


## User controllable from commandline
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
gWindow = WINDOW_ONES
gMinAmp4Clip = (1/256)*0.00001
gScanRangeNonOverlap = 0.5
gPrgLoopCnt = 8192
gXRes = 512
gPltCompress = PLTCOMPRESS_AVG
gCurScanCumuMode = CUMUMODE_AVG
gbGrid = True


## Needs to be edited directly below, if need to change
gZeroSpanFftDispProcMode = 'LogNoGain'
gScanRangeFftDispProcMode = 'LogNoGain'
gScanRangeClipProcMode = 'HistLowClip'
gScanRangeClipProcMode = 'Clip2MinAmp'
gPltCompressHM = PLTCOMPRESS_MAX


## Controlled from GUI
gbDataMin = True
gbDataMax = True
gbDataAvg = True
gbDataCur = True


'''
conv = [0.25, 0.5, 0.25]
conv = [0.05, 0.1, 0.2, 0.3, 0.2, 0.1, 0.05]
conv = [0.1, 0.2, 0.4, 0.2, 0.1]
conv = np.kaiser(512,64)
conv = np.kaiser(256,64)
conv = np.kaiser(128,64)
conv = np.kaiser(64,64)
DataProcEndAdj = 1.5
'''
DataProcConv = np.kaiser(128,64)
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
    elif dataProc == 'Conv':
        conv = DataProcConv
        #winAdj = len(conv)/np.sum(conv)
        #winAdj = 10*np.log10(winAdj)
        vals = np.convolve(vals, conv, mode='same')
        avg = np.average(vals)
        vals[:12] = avg
        vals[-12:] = avg
    return vals


def data_cumu(d, mode, curVals, cStart, cEnd, newVals, nStart, nEnd):
    '''
    Cumulate data from different buffers using one of possible logics

    Raw: copy values in newVals buffer into curVals buffer at specified offsets.
    Avg: average the values between curVals and newVals buffer and store into curVals.
    Max: copy the larger value between curVals and newVals into curVals.
    Min: copy the smaller value between curVals and newVals into curVals.
    '''
    if type(curVals) == type(None):
        return np.copy(newVals)
    if mode == CUMUMODE_RAW:
        curVals[cStart:cEnd] = newVals[nStart:nEnd]
    elif mode == CUMUMODE_AVG:
        curVals[cStart:cEnd] += newVals[nStart:nEnd]
        curVals[cStart:cEnd] /= 2
    elif mode == CUMUMODE_MAX:
        curVals[cStart:cEnd] = np.max([newVals[nStart:nEnd], curVals[cStart:cEnd]], axis=0)
    elif mode == CUMUMODE_MIN:
        curVals[cStart:cEnd] = np.min([newVals[nStart:nEnd], curVals[cStart:cEnd]], axis=0)
    else:
        msg = "ERROR: Unknown cumuMode [{}], Quiting...".format(mode)
        prg_quit(d, msg)
    return curVals


def fftvals_dispproc(d, vals, fftDispProcMode, infTo=None):
    '''
    Process fft value wrt displaying it.
    '''
    modes = fftDispProcMode.split('.')
    for mode in modes:
        if mode == 'Raw':
            continue
        elif mode == 'LogNoGain':
            vals = data_proc(d, vals, 'LogNoGain', infTo)
        elif mode == 'HistLowClip':
            vals = data_proc(d, vals, 'HistLowClip')
        else:
            msg = "ERROR: Unknown DispProcMode [{}], Quitting...".format(mode)
            prg_quit(d, msg)
    return vals


def _data_plotcompress(d, data, mode):
    '''
    Process the given data in the specified way.

    Raw: Dont modify the data, return it as it is.
    Conv: Convolve the data before returning it.
    Max: Divide the data into d['xRes'] number of groups,
         and return the max from each group.
    Min: Divide the data into d['xRes'] number of groups,
         and return the min from each group.
    Avg: Divide the data into d['xRes'] number of groups,
         and return the avg of each group.

    data length should be divisible into equal size groups using d['xRes']
    for Max/Min/Avg.
    '''
    if mode == PLTCOMPRESS_RAW:
        return data
    elif mode == PLTCOMPRESS_CONV:
        return data_proc(d, data, 'Conv')
    elif (mode == PLTCOMPRESS_MAX) or (mode == PLTCOMPRESS_AVG):
        rows = d['xRes']
        cols = len(data)//rows
        if cols == 0:
            return data
        tData = data.reshape(rows, cols)
        if mode == PLTCOMPRESS_MAX:
            data = np.max(tData, axis=1)
        if mode == PLTCOMPRESS_MIN:
            data = np.min(tData, axis=1)
        elif mode == PLTCOMPRESS_AVG:
            data = np.average(tData, axis=1)
        return data
    else:
        prg_quit(d, "ERROR:_data_plotcompress: Unknown mode [{}]".format(mode))


def data_plotcompress(d, xData, yData, mode=None):
    '''
    Process and or Reduce the amount of data, while still trying to maintain any
    significant values in the data, if that is what user wants.

    xData is blindly averaged, while yData is processed as specified by user.
    '''
    if type(mode) == type(None):
        mode = d['pltCompress']
    if mode == PLTCOMPRESS_RAW:
        return xData, yData
    if mode == PLTCOMPRESS_CONV:
        yData = _data_plotcompress(d, yData, mode)
        return xData, yData
    xData = _data_plotcompress(d, xData, PLTCOMPRESS_AVG)
    yData = _data_plotcompress(d, yData, mode)
    return xData, yData


def data_2d_plotcompress(d, data, mode=None):
    '''
    If requested, then reduce the number of elements in a 2D data set,
    by merging adjacent cols of each row, into a smaller subset of cols.
    '''
    if type(mode) == type(None):
        mode = d['pltCompressHM']
    if mode == PLTCOMPRESS_RAW:
        return data
    rows,cols = data.shape
    newData = []
    for r in range(rows):
        newData.append(_data_plotcompress(d, data[r,:], mode))
    return np.array(newData)


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
    print("PlotHighs: Freqs {} to {} : delta4Marking {} : min {} max {}".format(freqs[0], freqs[-1], delta4Marking, np.min(levels), np.max(levels)))
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
            d['AxFreqs'].text(0.1,1.0-0.1*(cntMarked+1), "{}:{}".format(round(curFreq/1e6,8), round(curLevel,2)))
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

    The logic supports reading non multiples of gSdrReadUnit. However rtlsdr
    doesnt seem to support reading non power of 2 buffer sizes over libusb,
    SO ensure to read buffers of length which satisfy this requirement, ELSE
    the logic will read the nearest higher power of 2 amount of samples and
    discard the data towards the end, which is not needed.
    '''
    if length > gSdrReadUnit:
        loopCnt = length//gSdrReadUnit
        remaining = length % gSdrReadUnit
        readLength = gSdrReadUnit
    else:
        remaining = length
        loopCnt = 0
        readLength = 0
    samples = np.zeros(length, dtype=complex)
    for i in range(loopCnt):
        iStart = i*readLength
        iEnd = (i+1)*readLength
        samples[iStart:iEnd] = sdr.read_samples(readLength)
    if remaining > 0:
        iStart = gSdrReadUnit*loopCnt
        iEnd = iStart+remaining
        adjustedRead = 2**np.ceil(np.log2(remaining))
        if adjustedRead != remaining:
            print("WARN:WillDiscard:Reading {} for {} from {} to {}".format(adjustedRead, remaining, iStart, iEnd))
        samples[iStart:iEnd] = sdr.read_samples(adjustedRead)[0:remaining]
    return samples


gbUsePSD = False
def sdr_curscan(d):
    '''
    Scan the currently set freq band (upto max sampling rate supported).
    Inturn normalise the fft on captured samples.

    It does a overlapped sliding over the captured samples, with a per
    fft window of d['fftSize'].

    d['fullSize'] is the amount of data captured over which overlapped sliding
    is done.

    Based on d['window'], a user selected window function may be applied to the data
    before fft. The result is compensated wrt windowing related loss of amplitude.

    As IQ data is what is got from the hardware, so both +ve and -ve freqs
    are used to get the embedded signals in the sample data.
    '''
    numLoops = int(d['fullSize']/(d['fftSize']*d['curScanNonOverlap']))
    #print("curscan: numLoops[{}] fullSize[{}]".format(numLoops, d['fullSize']))
    samples = sdr_read(d['sdr'], d['fullSize'])
    fftAll = None
    win = d['theWin']
    winAdj = len(win)/np.sum(win)
    if d['bUsePSD']:
        noverlap = d['fftSize']*(1-d['curScanNonOverlap'])
        '''
        p = plt.specgram(samples, NFFT=d['fftSize'], window=win, noverlap=noverlap, Fs=2, mode='magnitude')
        plt.cla()
        fftAll = np.average(p[0], axis=1)
        '''
        p = plt.psd(samples, NFFT=d['fftSize'], window=win, noverlap=noverlap)
        plt.cla()
        fftAll = p[0]
        return fftAll
    for i in range(numLoops):
        iStart = int(i*d['fftSize']*d['curScanNonOverlap'])
        iEnd = iStart + d['fftSize']
        tSamples = samples[iStart:iEnd]
        if len(tSamples) < d['fftSize']:
            break
        fftN = winAdj*2*abs(np.fft.fft(tSamples*win))/len(tSamples)
        if i == 0:
            fftAll = fftN
        else:
            fftAll = data_cumu(d, d['curScanCumuMode'], fftAll, 0, len(fftAll), fftN, 0, len(fftN))
    fftAll = np.fft.fftshift(fftAll)
    return fftAll


def _adj_siglvls(d, fftPr):
    if d['AdjSigLvls'] != '':
        fftMax = d['Fft.Max'] - d['Fft.Adj']
        fftMin = d['Fft.Min'] - d['Fft.Adj']
        fftAvg = d['Fft.Avg'] - d['Fft.Adj']
        fftPrTmp = fftPr - d['Fft.Adj']
    else:
        fftMax = d['Fft.Max']
        fftMin = d['Fft.Min']
        fftAvg = d['Fft.Avg']
        fftPrTmp = fftPr
    return fftMax, fftMin, fftAvg, fftPrTmp


def _heatmap_create(d, data):
    #hm = d['AxHeatMap'].imshow(data, extent=(0,1, 0,1), aspect='auto')
    hm = d['AxHeatMap'].imshow(data, extent=(0,1, 0,1), aspect='auto', interpolation='bicubic', picker=True)
    d['AxHeatMap'].set_xticks([0, 0.25, 0.5, 0.75, 1])
    f25 = d['startFreq'] + (d['centerFreq'] - d['startFreq'])/2
    f75 = d['centerFreq'] + (d['endFreq'] - d['centerFreq'])/2
    d['AxHeatMap'].set_xticklabels([d['startFreq'], f25, d['centerFreq'], f75, d['endFreq']])
    d['AxHeatMap'].set_xlabel("Freqs")
    d['AxHeatMap'].set_ylabel("ScanHistory")
    return hm


def zero_span(d):
    '''
    Repeatadly keep scanning a specified freq band, which is configured
    by default to a relatively safe sampling rate supported by the hardware.

    Inturn the scanned result is processed into Max, Min, Avg and Cur buckets.

    Next these processed buckets are shown to the user as part of the SigLvls
    plot, provided the user as requested the same.

    Also a history of instanteneous(Cur) signal levels is shown as a heat map.
    '''
    d['timeWasStr'] = None
    d['Fft.Max'] = None
    d['Fft.Min'] = None
    d['Fft.Avg'] = None
    d['Fft.Cur'] = None
    d['sdr'], bSdrSetup = sdr_setup(d['sdr'], d['centerFreq'], d['samplingRate'], d['gain'])
    freqs = np.fft.fftfreq(d['fftSize'],1/d['samplingRate']) + d['centerFreq']
    freqs = np.fft.fftshift(freqs)
    print("ZeroSpan: min[{}] max[{}]".format(min(freqs), max(freqs)))
    if d['bPltHeatMap']:
        maxHM = 128
        if (d['pltCompressHM'] in [ PLTCOMPRESS_MAX, PLTCOMPRESS_MIN, PLTCOMPRESS_AVG]):
            if d['fftSize'] > d['xRes']:
                d['PltHeatMapWidth'] = d['xRes']
            else:
                d['PltHeatMapWidth'] = d['fftSize']
        else:
            d['PltHeatMapWidth'] = d['fftSize']
        fftHM = np.zeros((maxHM, d['PltHeatMapWidth']))
        indexHM = 0
        hm = _heatmap_create(d, fftHM)
    prevTime = time.time()
    for i in range(d['prgLoopCnt']):
        if d['cmd.stop']:
            break
        curTime = time.time()
        print("ZeroSpan:{}:{}".format(i, curTime-prevTime))
        prevTime = curTime
        fftCur = sdr_curscan(d)
        #print("DBUG:ZeroSpan:fftCur:0:{}:mid:{}".format(fftCur[0], fftCur[(len(fftCur)//2)-1:(len(fftCur)//2)+1]))
        #print("DBUG:ZeroSpan:fftCur:0:{}:mid:{}".format(fftCur[0], fftCur[(len(fftCur)//2)-1:(len(fftCur)//2)+1]))
        fftPr = fftvals_dispproc(d, fftCur, gZeroSpanFftDispProcMode)
        d['Fft.Cur'] = fftPr
        if d['bDataMax']:
            d['Fft.Max'] = data_cumu(d, CUMUMODE_MAX, d['Fft.Max'], 0, len(fftPr), fftPr, 0, len(fftPr))
        if d['bDataMin']:
            d['Fft.Min'] = data_cumu(d, CUMUMODE_MIN, d['Fft.Min'], 0, len(fftPr), fftPr, 0, len(fftPr))
        if d['bDataAvg']:
            d['Fft.Avg'] = data_cumu(d, CUMUMODE_AVG, d['Fft.Avg'], 0, len(fftPr), fftPr, 0, len(fftPr))
        #print("DBUG:ZeroSpan:fftPr:0:{}:mid:{}".format(fftPr[0], fftPr[(len(fftPr)//2)-1:(len(fftPr)//2)+1]))
        fftMax, fftMin, fftAvg, fftPrTmp = _adj_siglvls(d, fftPr)
        if d['bPltHeatMap']:
            fftHM[indexHM,:] = _data_plotcompress(d, fftPrTmp, d['pltCompressHM'])
            hm.set_data(fftHM)
            hm.autoscale()
            plt.draw()
            indexHM = (indexHM + 1) % maxHM
        if d['bPltLevels']:
            d['AxLevels'].cla()
            if d['bGrid']:
                d['AxLevels'].grid(True)
            if d['bDataMax']:
                xFreqs, yLvls = data_plotcompress(d, freqs, fftMax)
                d['AxLevels'].plot(xFreqs, yLvls, 'r')
            if d['bDataMin']:
                xFreqs, yLvls = data_plotcompress(d, freqs, fftMin)
                d['AxLevels'].plot(xFreqs, yLvls, 'y')
            if d['bDataAvg']:
                xFreqs, yLvls = data_plotcompress(d, freqs, fftAvg)
                d['AxLevels'].plot(xFreqs, yLvls, 'g')
            if d['bDataCur']:
                xFreqs, yLvls = data_plotcompress(d, freqs, fftPrTmp)
                d['AxLevels'].plot(xFreqs, yLvls, 'b')
            if d['timeWasStr'] != None:
                d['AxLevels'].set_xlabel(d['timeWasStr'])
            plt.draw()
            plot_highs(d, xFreqs, yLvls)
        plt.pause(0.0001)



gZeroSpanSaveFile = '/tmp/zerospan.save'
def zero_span_save(d):
    f = open(d['zeroSpanSaveFile'], "wb+")
    pickle.dump(d['centerFreq'], f)
    pickle.dump(d['samplingRate'], f)
    pickle.dump(d['gain'], f)
    d['sdr'], bSdrSetup = sdr_setup(d['sdr'], d['centerFreq'], d['samplingRate'], d['gain'])
    prevTime = time.time()
    for i in range(d['prgLoopCnt']):
        if d['cmd.stop']:
            break
        curTime = time.time()
        print("ZeroSpanSave:{}:{}".format(i, curTime-prevTime))
        prevTime = curTime
        fftCur = sdr_curscan(d)
        pickle.dump(curTime, f)
        pickle.dump(fftCur, f)
    f.close()



def zero_span_play(d):
    '''
    Returns the next time and fft result data from file which contains the saved data.
    Time is returned as seconds since epoch as well as in string format.
    '''
    d['timeWas'] = pickle.load(d['zeroSpanFile'])
    timeWasMilli = int((d['timeWas'] - int(d['timeWas']))*1000)
    timeWas = time.strftime("%Y%m%d%Z%H%M%S", time.gmtime(d['timeWas']))
    d['timeWasStr'] = "{}.{:03}".format(timeWas, timeWasMilli)
    print("INFO:zeroSpanPlay:timeWas:{}".format(d['timeWasStr']))
    return pickle.load(d['zeroSpanFile'])



gbScanRangeBaseDataIsRaw = False
def _scan_range(d, freqsAll, fftAll, runCount=-1):
    '''
    Scan a specified range, this can be larger than the freq band
    that can be sampled/scanned by the hardware in one go, in which
    case it will do multiple scans to cover the full specified range.

    These multiple scans to cover the full range, in turn can be either
    overlapped or not, as decided by scanRangeNonOverlap.

    It generates multiple data sets from the scanned signal levels
    like Max, Min, Avg and Cur.

    Cur is a average of the overlaped scanning during sliding window
    over the full frequency range.

    Max,Min and Avg is generated from Fft.Cur or fftPr.
    HeatMap is generated from Fft.Avg data.
    '''
    freqSpan = d['samplingRate']
    if ((freqSpan*d['scanRangeNonOverlap'])%1) != 0:
        msg = "ERROR: freqSpan [{}] x scanRangeNonOverlap [{}] is not int".format(freqSpan, d['scanRangeNonOverlap'])
        prg_quit(d, msg)
    if ((d['fftSize']*d['scanRangeNonOverlap'])%1) != 0:
        msg = "ERROR: fftSize[{}] x scanRangeNonOverlap [{}] is not int".format(d['fftSize'], d['scanRangeNonOverlap'])
        prg_quit(d, msg)
    curFreq = d['startFreq'] + freqSpan/2
    startFreq = curFreq - freqSpan/2
    #endFreq = curFreq + freqSpan/2
    print("_scanRange: start:{} end:{} samplingRate:{}".format(d['startFreq'], d['endFreq'], d['samplingRate']))
    totalFreqs = d['endFreq'] - d['startFreq']
    numGroups = int(totalFreqs/freqSpan)
    totalEntries = numGroups * d['fftSize']
    print("_scanRange: totalFreqs:{} numGroups:{} totalEntries:{}".format(totalFreqs, numGroups, totalEntries))
    if type(freqsAll) == type(None):
        fftAll = np.ones(totalEntries) * d['minAmp4Clip']
        d['Fft.Cur'] = fftvals_dispproc(d, fftAll, gScanRangeFftDispProcMode, infTo=0)
        d['Fft.Max'] = np.copy(d['Fft.Cur'])
        d['Fft.Avg'] = np.copy(d['Fft.Cur'])
        fftAll = np.ones(totalEntries)
        d['Fft.Min'] = fftvals_dispproc(d, fftAll, gScanRangeFftDispProcMode, infTo=0)
        freqsAll = np.fft.fftshift(np.fft.fftfreq(totalEntries, 1/(numGroups*freqSpan)) + d['startFreq'] + (numGroups*freqSpan)/2)
        if d['bPltHeatMap']:
            d['fftHMMax'] = 128
            d['fftHMIndex'] = 0
            hmData = np.ones([d['fftHMMax'], totalEntries]) * d['minAmp4Clip']
            d['fftHM'] = data_2d_plotcompress(d, hmData, d['pltCompressHM'])
    if runCount == 0:
        cumuMode4Avg = CUMUMODE_RAW
    else:
        cumuMode4Avg = CUMUMODE_AVG
    i = 0
    iOldEnd = 0
    while startFreq < d['endFreq']:
        iStart = int(i*d['fftSize']*d['scanRangeNonOverlap'])
        iEnd = iStart+d['fftSize']
        iDone = int((i+1)*d['fftSize']*d['scanRangeNonOverlap'])
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
        fftPr = fftvals_dispproc(d, np.copy(fftCur), gScanRangeFftDispProcMode, infTo=0)
        if d['bDataCur'] or True:
            sRawStart = sStart + (d['fftSize'] - (iEnd - iOldEnd))
            d['Fft.Cur'] = data_cumu(d, CUMUMODE_RAW, d['Fft.Cur'], iOldEnd, iEnd, fftPr, sRawStart, sEnd)
            if iOldEnd != 0:
                if iOldEnd > totalEntries:
                    iOldEnd = totalEntries
                sAvgEnd = sStart + (iOldEnd - iStart) # iOldEnd - iStart # Both are same
                d['Fft.Cur'] = data_cumu(d, CUMUMODE_AVG, d['Fft.Cur'], iStart, iOldEnd, fftPr, sStart, sAvgEnd)
            iOldEnd = iEnd
        if d['bScanRangeBaseDataIsRaw']:
            dstStart = iStart
            dstEnd = iEnd
            srcData = fftPr
            srcStart = sStart
            srcEnd = sEnd
        else:
            dstStart = iStart
            dstEnd = iDone
            srcData = d['Fft.Cur']
            srcStart = iStart
            srcEnd = iDone
        if d['bDataMax']:
            d['Fft.Max'] = data_cumu(d, CUMUMODE_MAX, d['Fft.Max'], dstStart, dstEnd, srcData, srcStart, srcEnd)
        if d['bDataMin']:
            d['Fft.Min'] = data_cumu(d, CUMUMODE_MIN, d['Fft.Min'], dstStart, dstEnd, srcData, srcStart, srcEnd)
        if d['bDataAvg'] or True:
            d['Fft.Avg'] = data_cumu(d, cumuMode4Avg, d['Fft.Avg'], dstStart, dstEnd, srcData, srcStart, srcEnd)
        fftMax, fftMin, fftAvg, fftCurAdj = _adj_siglvls(d, d['Fft.Cur'])
        if d['bPltLevels']:
            d['AxLevels'].clear()
            #d['AxLevels'].plot(startFreq, -30, "bo")
            #d['AxLevels'].plot(curFreq, -30, "ro")
            #d['AxLevels'].plot(endFreq, -30, "bo")
            if d['bGrid']:
                d['AxLevels'].grid(True)
            if d['bDataMax']:
                xFreqs, yLvls = data_plotcompress(d, freqsAll, fftMax)
                d['AxLevels'].plot(xFreqs, yLvls, 'r')
            if d['bDataMin']:
                xFreqs, yLvls = data_plotcompress(d, freqsAll, fftMin)
                d['AxLevels'].plot(xFreqs, yLvls, 'y')
            if d['bDataAvg']:
                xFreqs, yLvls = data_plotcompress(d, freqsAll, fftAvg)
                d['AxLevels'].plot(xFreqs, yLvls, 'g')
            if d['bDataCur']:
                xFreqs, yLvls = data_plotcompress(d, freqsAll, fftCurAdj)
                d['AxLevels'].plot(xFreqs, yLvls, 'b')
        curFreq += freqSpan*d['scanRangeNonOverlap']
        startFreq = curFreq - freqSpan/2
        #endFreq = curFreq + freqSpan/2
        plt.pause(0.0001)
        i += 1
    if d['bPltLevels']:
        plot_highs(d, xFreqs, yLvls)
    if d['bPltHeatMap']:
        d['fftHM'][d['fftHMIndex'], :] = _data_plotcompress(d, fftAvg, d['pltCompressHM'])
    return freqsAll, fftAll


def scan_range(d):
    freqs = None
    ffts = None
    # Adjust endFreq such that start-end is a multiple of samplingRate, if required
    freqBands = (d['endFreq'] - d['startFreq'])/d['samplingRate']
    if (freqBands % 1) != 0:
        d['orig.EndFreq'] = d['endFreq']
        d['endFreq'] = d['startFreq'] + np.ceil(freqBands)*d['samplingRate']
        d['centerFreq'] = d['startFreq'] + ((d['endFreq'] - d['startFreq'])/2)
        print("WARN:scanRange:Adjusting endFreq: orig [{}] adjusted [{}], so that fullRange is Multiple of samplingRate/freqBand [{}]".format(d['orig.EndFreq'], d['endFreq'], d['samplingRate']))
        #input("Press any key to continue...")
    if d['bPltHeatMap']:
        hm = _heatmap_create(d, np.zeros([3,3]))
    prevTime = time.time()
    for i in range(d['prgLoopCnt']):
        if d['cmd.stop']:
            break
        curTime = time.time()
        print("scanRange:{}:{}".format(i, curTime-prevTime))
        prevTime = curTime
        #print_info(d)
        freqs, ffts = _scan_range(d, freqs, ffts, i)
        if d['bPltHeatMap']:
            #print("DBUG:scanRange: min[{}] max[{}]".format(np.min(d['fftHM']), np.max(d['fftHM'])))
            hm.set_data(d['fftHM'])
            hm.autoscale()
            plt.pause(0.001)
            d['fftHMIndex'] = (d['fftHMIndex'] + 1) % d['fftHMMax']

gSaveSigLvls = ''
gAdjSigLvls = ''
def _save_siglvls(d):
    try:
        if d['SaveSigLvls'] != '':
            f = open(d['SaveSigLvls'], "wb+")
            pickle.dump(d['startFreq'], f)
            pickle.dump(d['endFreq'], f)
            pickle.dump(d['Fft.Avg'], f)
            f.close()
            print("INFO:_save_siglvls: success...", d['SaveSigLvls'])
        else:
            print("INFO:_save_siglvls: ignoring...")
    except:
        print("WARN:_save_siglvls: Failed...", d['SaveSigLvls'])


def _load_siglvls(d):
    try:
        if d['AdjSigLvls'] != '':
            f = open(d['AdjSigLvls'], "rb")
            startFreq = pickle.load(f)
            endFreq = pickle.load(f)
            d['Fft.Adj'] = pickle.load(f)
            f.close()
            if (startFreq == d['startFreq']) and (endFreq == d['endFreq']):
                print("INFO:_load_siglvls: success...", d['AdjSigLvls'])
            else:
                input("ERRR:_load_siglvls:{}:savedRange[{}-{}] curFreqRange[{}-{}]".format(d['AdjSigLvls'], startFreq, endFreq, d['startFreq'], d['endFreq']))
                d['AdjSigLvls'] = ''
        else:
            print("INFO:_load_siglvls: ignoring...")
    except:
        print("WARN:_load_siglvls: Failed...", d['AdjSigLvls'])
        d['AdjSigLvls'] = ''


def _arg_boolean(value):
    if value.upper() == "TRUE":
        return True
    else:
        return False


def handle_args(d):
    '''
    Initialise the global dictionary, as well as update it based on
    the arguments given by the user.
    '''
    d['prgMode'] = gPrgModeDefault
    d['samplingRate'] = gSamplingRate
    d['gain'] = gGain
    d['centerFreq'] = gCenterFreq
    d['fftSize'] = gFftSize
    d['curScanNonOverlap'] = gCurScanNonOverlap
    d['curScanCumuMode'] = gCurScanCumuMode
    d['window'] = gWindow
    d['minAmp4Clip'] = gMinAmp4Clip
    d['bPltHeatMap'] = gbPltHeatMap
    d['bPltLevels'] = gbPltLevels
    d['scanRangeNonOverlap'] = gScanRangeNonOverlap
    d['prgLoopCnt'] = gPrgLoopCnt
    d['xRes'] = gXRes
    d['pltCompress'] = gPltCompress
    d['pltHighsNumMarkers'] = gPltHighsNumMarkers
    d['pltHighsDelta4Marking'] = gPltHighsDelta4Marking
    d['pltHighsPause'] = gPltHighsPause
    d['pltCompressHM'] = gPltCompressHM
    d['SaveSigLvls'] = gSaveSigLvls
    d['AdjSigLvls'] = gAdjSigLvls
    d['bDataMin'] = gbDataMin
    d['bDataMax'] = gbDataMax
    d['bDataAvg'] = gbDataAvg
    d['bDataCur'] = gbDataCur
    d['bGrid'] = gbGrid
    d['bUsePSD'] = gbUsePSD
    d['bScanRangeBaseDataIsRaw'] = gbScanRangeBaseDataIsRaw
    d['zeroSpanSaveFile'] = gZeroSpanSaveFile
    d['zeroSpanPlayFile'] = gZeroSpanSaveFile
    iArg = 1
    while iArg < len(sys.argv):
        curArg = sys.argv[iArg].upper()
        if (curArg in [PRGMODE_ZEROSPAN, PRGMODE_ZEROSPANSAVE, PRGMODE_ZEROSPANPLAY, PRGMODE_SCAN, PRGMODE_ALIAS_FMSCAN, PRGMODE_ALIAS_QUICKFULLSCAN]):
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
        elif (curArg == 'CURSCANCUMUMODE'):
            iArg += 1
            d['curScanCumuMode'] = sys.argv[iArg].upper()
        elif (curArg == 'SCANRANGENONOVERLAP'):
            iArg += 1
            d['scanRangeNonOverlap'] = float(sys.argv[iArg])
        elif (curArg == 'FFTSIZE'):
            iArg += 1
            d['fftSize'] = int(sys.argv[iArg])
        elif (curArg == 'XRES'):
            iArg += 1
            d['xRes'] = int(sys.argv[iArg])
        elif (curArg == 'BDATAMIN'):
            iArg += 1
            d['bDataMin'] = _arg_boolean(sys.argv[iArg])
        elif (curArg == 'BDATAMAX'):
            iArg += 1
            d['bDataMax'] = _arg_boolean(sys.argv[iArg])
        elif (curArg == 'BDATAAVG'):
            iArg += 1
            d['bDataAvg'] = _arg_boolean(sys.argv[iArg])
        elif (curArg == 'BDATACUR'):
            iArg += 1
            d['bDataCur'] = _arg_boolean(sys.argv[iArg])
        elif (curArg == 'PLTCOMPRESS'):
            iArg += 1
            d['pltCompress'] = sys.argv[iArg].upper()
        elif (curArg == 'WINDOW'):
            iArg += 1
            d['window'] = "WIN.{}".format(sys.argv[iArg].upper())
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
        elif (curArg == 'SAVESIGLVLS'):
            iArg += 1
            d['SaveSigLvls'] = sys.argv[iArg]
        elif (curArg == 'ADJSIGLVLS'):
            iArg += 1
            d['AdjSigLvls'] = sys.argv[iArg]
        elif (curArg == 'BGRID'):
            iArg += 1
            d['bGrid'] = _arg_boolean(sys.argv[iArg])
        elif (curArg == 'BUSEPSD'):
            iArg += 1
            d['bUsePSD'] = _arg_boolean(sys.argv[iArg])
        elif (curArg == 'BSCANRANGEBASEDATAISRAW'):
            iArg += 1
            d['bScanRangeBaseDataIsRaw'] = _arg_boolean(sys.argv[iArg])
        elif (curArg == 'ZEROSPANSAVEFILE'):
            iArg += 1
            d['zeroSpanSaveFile'] = sys.argv[iArg]
        elif (curArg == 'ZEROSPANPLAYFILE'):
            iArg += 1
            d['zeroSpanPlayFile'] = sys.argv[iArg]
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
        d['fftSize'] = 64
        d['pltCompress'] = PLTCOMPRESS_RAW
    if d['prgMode'] == PRGMODE_SCAN:
        d['centerFreq'] = d['startFreq'] + ((d['endFreq'] - d['startFreq'])/2)
    else: # ZeroSpan or related i.e save or play
        d['startFreq'] = d['centerFreq'] - d['samplingRate']/2
        d['endFreq'] = d['centerFreq'] + d['samplingRate']/2
    if d['fftSize'] < (d['samplingRate']//8):
        d['fullSize'] = d['fftSize'] * gFft2FullMult4Less
    else:
        d['fullSize'] = d['fftSize'] * gFft2FullMult4More
    #if (d['fullSize'] > gSdrReadUnit) and ((d['fullSize'] % gSdrReadUnit) != 0):
    #    prg_quit(d, "ERROR:fullSize[{}] Not multiple of gSdrReadUnit[{}]".format(d['fullSize'], gSdrReadUnit))
    d['WIN.HAMMING'] = np.hamming(d['fftSize'])
    d['WIN.HANNING'] = np.hanning(d['fftSize'])
    d['WIN.KAISER'] = np.kaiser(d['fftSize'], 64)
    d['WIN.ONES'] = np.ones(d['fftSize'])
    d['theWin'] = d[d['window']]
    # Adjust xRes, if required
    if d['xRes'] > d['fftSize']:
        print("WARN:fftSize[{}] < xRes[{}], setting xRes to fftSize".format(d['fftSize'], d['xRes']))
        d['xRes'] = d['fftSize']
    else:
        minXRes = 300
        if d['fftSize'] % d['xRes'] != 0:
            for i in range(int(d['fftSize']/minXRes), 0, -1):
                if d['fftSize'] % i == 0:
                    newXRes = d['fftSize'] // i
                    input("WARN:fftSize[{}] NotMultipleOf xRes[{}], setting xRes to {}".format(d['fftSize'], d['xRes'], newXRes))
                    break
            d['xRes'] = newXRes



def print_info(d):
    print("INFO: startFreq[{}] centerFreq[{}] endFreq[{}]".format(d['startFreq'], d['centerFreq'], d['endFreq']))
    print("INFO: samplingRate[{}], gain[{}], bUsePSD[{}]".format(d['samplingRate'], d['gain'], d['bUsePSD']))
    print("INFO: fullSize[{}], fftSize[{}], curScanCumuMode[{}], window[{}]".format(d['fullSize'], d['fftSize'], d['curScanCumuMode'], d['window']))
    print("INFO: minAmp4Clip[{}], curScanNonOverlap[{}], scanRangeNonOverlap[{}], bScanRangeBaseDataIsRaw[{}]".format(
            d['minAmp4Clip'], d['curScanNonOverlap'], d['scanRangeNonOverlap'], d['bScanRangeBaseDataIsRaw']))
    print("INFO: prgMode [{}], prgLoopCnt[{}], bPltLevels[{}],  bPltHeatMap[{}]".format(d['prgMode'], d['prgLoopCnt'], d['bPltLevels'], d['bPltHeatMap']))
    print("INFO: pltHighsNumMarkers[{}], pltHighsDelta4Marking[{}], pltHighsPause[{}]".format(d['pltHighsNumMarkers'], d['pltHighsDelta4Marking'], d['pltHighsPause']))
    print("INFO: xRes [{}], bGrid [{}], pltCompress [{}], pltCompressHM [{}]".format(d['xRes'], d['bGrid'], d['pltCompress'], d['pltCompressHM']))
    print("INFO: SaveSigLvls [{}], AdjSigLvls [{}]; zeroSpanSaveFile[{}], zeroSpanPlayFile[{}]".format(d['SaveSigLvls'], d['AdjSigLvls'], d['zeroSpanSaveFile'], d['zeroSpanPlayFile']))
    print("INFO: bDataMax [{}], bDataMin [{}], bDataAvg[{}], bDataCur [{}]".format(d['bDataMax'], d['bDataMin'] , d['bDataAvg'], d['bDataCur']))



def prg_quit(d, msg = None, tryExit=True):
    if type(msg) != type(None):
        print(msg)
    d['cmd.stop'] = True
    if tryExit:
        sys.exit()


def _update_boolbtn(btn, btnBool, text):
    if btnBool:
        btn.label.set_text('{}[x]'.format(text))
    else:
        btn.label.set_text('{}[ ]'.format(text))


def update_boolbtns(d):
    if not (d['bDataMin'] or d['bDataMax'] or d['bDataAvg'] or d['bDataCur']):
        d['bDataAvg'] = True
    _update_boolbtn(d['BtnLevels'], d['bPltLevels'], 'Levels')
    _update_boolbtn(d['BtnHeatMap'], d['bPltHeatMap'], 'HeatMap')
    _update_boolbtn(d['BtnPause'], d['pltHighsPause'], 'Pause')
    _update_boolbtn(d['BtnMinLvls'], d['bDataMin'], 'MinLvls')
    _update_boolbtn(d['BtnMaxLvls'], d['bDataMax'], 'MaxLvls')
    _update_boolbtn(d['BtnAvgLvls'], d['bDataAvg'], 'AvgLvls')
    _update_boolbtn(d['BtnCurLvls'], d['bDataCur'], 'CurLvls')


def event_pause(event):
    if gD['pltHighsPause']:
        gD['pltHighsPause'] = False
    else:
        gD['pltHighsPause'] = True
    update_boolbtns(gD)


def event_levels(event):
    if gD['bPltLevels']:
        gD['bPltLevels'] = False
    else:
        gD['bPltLevels'] = True
    update_boolbtns(gD)


def event_minlvls(event):
    if gD['bDataMin']:
        gD['bDataMin'] = False
    else:
        gD['bDataMin'] = True
    update_boolbtns(gD)


def event_maxlvls(event):
    if gD['bDataMax']:
        gD['bDataMax'] = False
    else:
        gD['bDataMax'] = True
    update_boolbtns(gD)


def event_avglvls(event):
    if gD['bDataAvg']:
        gD['bDataAvg'] = False
    else:
        gD['bDataAvg'] = True
    update_boolbtns(gD)


def event_curlvls(event):
    if gD['bDataCur']:
        gD['bDataCur'] = False
    else:
        gD['bDataCur'] = True
    update_boolbtns(gD)


def event_heatmap(event):
    if gD['bPltHeatMap']:
        gD['bPltHeatMap'] = False
    else:
        gD['bPltHeatMap'] = True
    update_boolbtns(gD)


def event_quit(event):
    gD['BtnQuit'].label.set_text("QuitWait")
    prg_quit(gD, "INFO:QuitClick: Quiting on user request...", False)


def handle_pick(event):
    '''
    Show the freq corresponding to position clicked on the heatmap.

    BCas, once xticks are set, the default logic to show x and y loc
    on the figure window, no longer works for x.
    '''
    me = event.mouseevent
    freq = gD['startFreq'] + (gD['endFreq']-gD['startFreq'])*me.xdata
    print("INFO:PickEvent:HeatMap:Freq:", freq)
    '''
    #print(me.x, me.y, me.xdata, me.ydata, freq)
    try:
        if gD['HMFreqText'] != None:
            gD['HMFreqText'].remove()
    except:
        pass
    gD['HMFreqText'] = gD['AxHeatMap'].text(0,0, freq)
    '''
    gD['AxHeatMap'].set_xlabel("Freqs [ClickedFreq:{}]".format(freq))


def plt_figures(d):
    plt.ion()
    # 8,5 => [[2,4],[2,1]], [[2,5]]
    # 8,5 => [[4,4],[4,1]]
    # 8,5 => [[4,5]]
    f = plt.figure("kSpecAnal", figsize=(12, 8), constrained_layout=True)
    gs = f.add_gridspec(nrows=16, ncols=5)
    d['AxLevels'] = f.add_subplot(gs[:8,:4])
    d['AxFreqs'] = f.add_subplot(gs[:8,4])
    d['AxFreqs'].set_xlabel("Freqs - HighSigLvl")
    d['AxHeatMap'] = f.add_subplot(gs[8:16,:4])
    d['AxBtnLevels'] = f.add_subplot(gs[8,4])
    d['AxBtnHeatMap'] = f.add_subplot(gs[9,4])
    d['AxBtnMaxLvls'] = f.add_subplot(gs[10,4])
    d['AxBtnMinLvls'] = f.add_subplot(gs[11,4])
    d['AxBtnAvgLvls'] = f.add_subplot(gs[12,4])
    d['AxBtnCurLvls'] = f.add_subplot(gs[13,4])
    d['AxBtnPause'] = f.add_subplot(gs[14,4])
    d['AxBtnQuit'] = f.add_subplot(gs[15,4])
    d['AxFreqs'].set_xticks([])
    d['AxFreqs'].set_yticks([])
    d['BtnPause'] = plt.Button(d['AxBtnPause'], "Pause")
    d['BtnPause'].on_clicked(event_pause)
    d['BtnLevels'] = plt.Button(d['AxBtnLevels'], "Levels")
    d['BtnLevels'].on_clicked(event_levels)
    d['BtnHeatMap'] = plt.Button(d['AxBtnHeatMap'], "HeatMap")
    d['BtnHeatMap'].on_clicked(event_heatmap)
    d['BtnMinLvls'] = plt.Button(d['AxBtnMinLvls'], "MinLvls")
    d['BtnMinLvls'].on_clicked(event_minlvls)
    d['BtnMaxLvls'] = plt.Button(d['AxBtnMaxLvls'], "MaxLvls")
    d['BtnMaxLvls'].on_clicked(event_maxlvls)
    d['BtnAvgLvls'] = plt.Button(d['AxBtnAvgLvls'], "AvgLvls")
    d['BtnAvgLvls'].on_clicked(event_avglvls)
    d['BtnCurLvls'] = plt.Button(d['AxBtnCurLvls'], "CurLvls")
    d['BtnCurLvls'].on_clicked(event_curlvls)
    d['BtnQuit'] = plt.Button(d['AxBtnQuit'], "Quit")
    d['BtnQuit'].on_clicked(event_quit)
    f.canvas.mpl_connect('pick_event', handle_pick)
    update_boolbtns(d)


def handle_sigint(signum, stack):
    prg_quit(gD, "INFO:sigint: quiting on user request...")


def handle_signals(d):
    signal.signal(signal.SIGINT, handle_sigint)


def do_run(d):
    global sdr_curscan
    if d['prgMode'] == PRGMODE_SCAN:
        scan_range(d)
    elif d['prgMode'] == PRGMODE_ZEROSPANSAVE:
        zero_span_save(d)
    elif d['prgMode'] == PRGMODE_ZEROSPANPLAY:
        d['zeroSpanFile'] = f = open(d['zeroSpanPlayFile'], "rb")
        d['centerFreq'] = pickle.load(f)
        d['samplingRate'] = pickle.load(f)
        d['gain'] = pickle.load(f)
        print_info(d)
        sdr_curscan = zero_span_play
        zero_span(d)
        d['zeroSpanFile'].close()
    else:
        zero_span(d)


gD = {}
gD['cmd.stop'] = False
handle_args(gD)
_load_siglvls(gD)
print_info(gD)
handle_signals(gD)
plt_figures(gD)
gD['sdr'] = rtlsdr.RtlSdr()
do_run(gD)
gD['sdr'].close()



_save_siglvls(gD)
gD['BtnQuit'].label.set_text("QuitPress")
input("Press any key to quit...")


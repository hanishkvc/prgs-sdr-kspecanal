#!/bin/env python3
# kSpecAnal - A Spectrum Analyser
# v20170211_1342, HanishKVC
#

import rtlsdr
import cairo
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

gHeight = 720
gWidth = 1024

gDwellTime = 8e-3 # 20e-3
gFftSize = 2048 # 512
gbLivePlot = True
gbAdaptiveFixedYAxisZeroSpanPlot = False

gDebugLevel = 5
def dprint(dbgLvl, msg):
    if (dbgLvl <= gDebugLevel):
        print(msg)


def minmax_complex(data):
    tMinR = 0.0
    tMaxR = 0.0
    for i in data:
        if (tMinR > i.real):
            tMinR = i.real
        if (tMaxR < i.real):
            tMaxR = i.real
    return (tMinR, tMaxR)


def minmax_iq(dI, dQ):
    tMinR = min(min(dI), min(dQ))
    tMaxR = max(max(dI), max(dQ))
    return (tMinR, tMaxR)


def minmax(data):
    tMinR = min(data)
    tMaxR = max(data)
    return (tMinR, tMaxR)


def cairoplot_setup(sPlotSVGFile):
    pSetup = {}
    # Get Cairo Context
    pCRSurface = cairo.SVGSurface(sPlotSVGFile, gWidth, gHeight)
    pCR = cairo.Context(pCRSurface)
    pCR.set_source_rgb(100,100,200)
    pCR.paint()
    pSetup["CRSurface"] = pCRSurface
    pSetup["CRContext"] = pCR
    return pSetup


def cairoplot_dot(pSetup, x, y, bRect):
    pCR = pSetup["CRContext"]
    if (bRect):
        pCR.rectangle(x,y,x+1,y+1)
        pCR.fill()
        #pCR.stroke()
    else:
        pCR.move_to(x,y)
        pCR.show_text("*")


def cairoplot_xy(pSetup, x, y):
    pCR = pSetup["CRContext"]
    pYMin = pSetup["YMin"]
    pYMax = pSetup["YMax"]
    pYPerPixel = pSetup["YPerPixel"]
    cData = y-pYMin
    cY = cData / pYPerPixel
    #print("plot_xy: {},{} = {},{}".format(x,y,x,cY))
    pCR.set_source_rgb(100, 0, 0)
    cairoplot_dot(pSetup, x, cY, False)


def dprint_bytes(dbgLvl, dBytes):
    if dbgLvl > gDebugLevel:
        return
    for i in dBytes:
        print(i, end=" ")
    print()


def read_iq(sdr, numBytes):
    dBytes = sdr.read_bytes(numBytes)
    dprint(10,"min[{}], max[{}]".format(min(dBytes), max(dBytes)))
    dprint_bytes(10, dBytes)
    dI = np.array(dBytes[::2]) # dBytes[0::2]
    dQ = np.array(dBytes[1::2])
    dI = (dI - 127)/128
    dQ = (dQ - 127)/128
    dMinMax = minmax_iq(dI, dQ)
    return dI, dQ, dMinMax


def read_and_discard(sdr):
    data = sdr.read_samples(16*1024)
    #print(data)
    dMinMax = minmax_complex(data)
    return data, dMinMax


def rtlsdr_setup(sdr, centerFreq, sampleRate, gain):
    print("Setup: centerFreq[{}], sampleRate[{}], gain[{}]".format(centerFreq, sampleRate, gain))
    sdr.center_freq = centerFreq
    sdr.sample_rate = sampleRate
    sdr.gain = gain
    read_and_discard(sdr)
    return sdr


gArgs = {}
def handle_args():
    argCnt = len(sys.argv)
    if (argCnt >= 2):
        gArgs["mode"] = sys.argv[1]
    else:
        gArgs["mode"] = "ZERO_SPAN"
    if (gArgs["mode"] == "ZERO_SPAN"):
        if (argCnt >= 3):
            gArgs["centerFreq"] = sys.argv[2]
        else:
            gArgs["centerFreq"] = 91.1e6
        if (argCnt >= 4):
            gArgs["sampleRate"] = sys.argv[3]
        else:
            gArgs["sampleRate"] = 1e6
        #sdr.freq_correction = 0
        if (argCnt >= 5):
            gArgs["gain"] = float(sys.argv[4])
        else:
            gArgs["gain"] = 0
    elif (gArgs["mode"] == "SCAN"):
        if (argCnt >= 3):
            gArgs["startFreq"] = sys.argv[2]
        else:
            gArgs["startFreq"] = 30e6
        if (argCnt >= 4):
            gArgs["endFreq"] = sys.argv[3]
        else:
            gArgs["endFreq"] = 1e9
        if (argCnt >= 5):
            gArgs["sampleRate"] = sys.argv[4]
        else:
            gArgs["sampleRate"] = 1e6
        #sdr.freq_correction = 0
        if (argCnt >= 6):
            gArgs["gain"] = float(sys.argv[5])
        else:
            gArgs["gain"] = 0
    else:
        gArgs["mode"] = "FULL_SPAN"
        gArgs["startFreq"] = 28e6
        gArgs["endFreq"] = 1.7e9
        gArgs["sampleRate"] = 1e6
        gArgs["gain"] = 48.0


def rtlsdr_init():
    sdr = rtlsdr.RtlSdr()
    print(dir(sdr))
    print("GainValues/*DivBy10AndThenUse*/:{}".format(sdr.gain_values))
    return sdr


def rtlsdr_info(sdr):
    print("CenterFreq[{}], SamplingRate[{}], Gain[{}]".format(sdr.center_freq, sdr.sample_rate, sdr.gain))


def rtlsdr_curscan(sdr):
    totalSamples = sdr.sample_rate*gDwellTime
    decimationCnt = totalSamples/gFftSize
    if (decimationCnt < 1):
        print("WARN: dwellTime[{}] leads to totalSamples[{}] less than fftSize[{}], adjusting dwellTime to reach fftSize".format(gDwellTime, totalSamples, gFftSize))
    decimationCnt = int(np.ceil(decimationCnt))
    dprint(5,"decimationCnt[{}]".format(decimationCnt))
    dataFDC = 0
    dataF = np.zeros(gFftSize/2)
    for i in range(decimationCnt):
        dI, dQ, dMinMax = read_iq(sdr, gFftSize*2)
        data = dI + dQ
        dataFft = np.abs(np.fft.fft(data)/len(data))
        dataFft = dataFft[:len(dataFft)/2]*2
        dataFDC += dataFft[0]
        dataFft[0] = 1/255
        dataF += dataFft
    dataF = dataF/decimationCnt
    dataFDC = dataFDC/decimationCnt
    dMinMax = minmax(dataF)

    dprint(10, "DataFFT [{}]\n\tLength[{}]\n\tMinMax [{}]\n".format(dataF, len(dataF), dMinMax))
    dprint(2, "DataFFTDC [{}]\n\tLength[{}]\n\tMinMax [{}]\n".format(dataFDC, len(dataF), dMinMax))

    return data, dataF, dataFDC


def rtlsdr_scan(sdr, startFreq, endFreq, sampleRate, gain):
    if (gbLivePlot):
        pf = plt.figure()
        plt.show(block=False)
    freqSpan = sampleRate/2
    curFreq = startFreq
    dataFAll = np.array([])
    while curFreq < endFreq:
        rtlsdr_setup(sdr, curFreq, sampleRate, gain)
        data, dataF, dataFDC = rtlsdr_curscan(sdr)
        dataFAll = np.append(dataFAll, dataF)
        cairoplot_data(dataF, curFreq, sampleRate/2)
        if (gbLivePlot):
            plt.cla()
            freqAxis = np.linspace(startFreq, curFreq+freqSpan, len(dataFAll))
            plt.plot(freqAxis, dataFAll)
            plt.pause(0.001)
        curFreq += freqSpan
    if (gbLivePlot):
        input("Press any key...")
        plt.close(pf)
    return dataFAll


def rtlsdr_zerospan_repeat(sdr, centerFreq, sampleRate, gain):
    if (gbLivePlot):
        pf = plt.figure()
        plt.show(block=False)
    freqSpan = sampleRate/2
    startFreq = centerFreq - freqSpan/2
    endFreq = centerFreq + freqSpan/2
    dataFAll = np.array([])
    rtlsdr_setup(sdr, centerFreq, sampleRate, gain)
    iCnt = 0
    ymin = 1.0
    ymax = 0.0
    while True:
        data, dataF, dataFDC = rtlsdr_curscan(sdr)
        if (iCnt == 0):
            dataFAll = dataF
            iCnt += 1
        else:
            dataFAll += dataF
            dataFAll /= 2
        cairoplot_data(dataF, centerFreq, sampleRate/2)
        if (gbLivePlot):
            ymin = min(ymin,min(dataFAll))
            ymax = max(ymax,max(dataFAll))
            plt.cla()
            if (gbAdaptiveFixedYAxisZeroSpanPlot):
                plt.ylim(ymin, ymax)
            if (gbModeCentered):
                freqAxis = np.linspace(startFreq, endFreq, len(dataFAll))
            else:
                freqAxis = np.linspace(centerFreq, centerFreq+freqSpan, len(dataFAll))
            plt.plot(freqAxis, dataFAll)
            plt.pause(0.001)
    if (gbLivePlot):
        input("Press any key...")
        plt.close(pf)
    return dataFAll


bDisplaySVGFile = False
def cairoplot_data(dataF, freq, span):
    sPlotSVGFile = "/tmp/plot_{}.svg".format(freq)
    pSetup = cairoplot_setup(sPlotSVGFile)

    dMinMax = minmax(dataF)
    # Scale
    dMinMax = (dMinMax[0]*2, dMinMax[1]*2)
    # Calculate the Y Axis details
    pSetup["YMin"] = pYMin = dMinMax[0]
    pSetup["YMax"] = pYMax = dMinMax[1]
    pSetup["YPerPixel"] = (pYMax-pYMin)/gHeight
    print(pSetup)

    x=1
    for i in dataF:
        cairoplot_xy(pSetup, x, i)
        x+=1

    pSetup["CRSurface"].flush()
    pSetup["CRSurface"].finish()
    if (bDisplaySVGFile):
        os.system("display {}".format(sPlotSVGFile))


gbModeCentered = False
iDEBUG = 0
def plot_data(data, dataF, startOrCenterFreq, freqSpan, gain):
    if (data != None):
        plt.subplot(2,1,1)
        plt.plot(data)
        plt.subplot(2,1,2)
    fftBins = len(dataF)
    deltaFreq = freqSpan/fftBins
    if (gbModeCentered):
        startFreq = startOrCenterFreq - (freqSpan/2)
        centerFreq = startOfCenterFreq
        endFreq = startOrCenterFreq + (freqSpan/2)
    else:
        startFreq = startOrCenterFreq
        centerFreq = startOrCenterFreq + freqSpan/2
        endFreq = startOrCenterFreq + freqSpan
    freqAxis = np.linspace(startFreq, endFreq, fftBins)
    print("StartFreq[{}], CenterFreq[{}], EndFreq[{}], FreqSpan[{}]".format(startFreq, centerFreq, endFreq, freqSpan))
    print("\tNumOfFFTBins[{}], deltaFreq[{}]".format(fftBins, deltaFreq))
    print("\tNumOf XPoints[{}], YPoints[{}]".format(len(freqAxis), len(dataF)))
    if (iDEBUG == 0):
        plt.subplot(2,1,1)
        plt.plot(freqAxis, dataF)
        plt.title("dataF")

        plt.subplot(2,1,2)
    elif (iDEBUG == 1):
        plt.subplot(2,1,1)
        plt.semilogy(freqAxis, dataF)
        plt.subplot(2,1,2)
    elif (iDEBUG == 2):
        plt.subplot(2,2,1)
        plt.plot(freqAxis, dataF)

        plt.subplot(2,2,2)
        minAmp = 8*(1.0/255)
        #dprint(2,dataF)
        #dataF = np.clip(dataF, minAmp*0.25, 2.0)
        dataF1 = 10*np.log10(dataF/minAmp)
        #dprint(2,dataF)
        dataF1 = -gain + dataF1
        plt.plot(freqAxis, dataF1)
        plt.show()
    elif (iDEBUG == 3):
        plt.subplot(2,2,1)
        plt.plot(freqAxis, dataF)
        plt.title("dataF")

        plt.subplot(2,2,2)
        minAmp = (1.0/255)
        dataF1 = 10*np.log10(dataF/minAmp)
        plt.plot(freqAxis, dataF1)
        plt.title("10*log10 wrt 1/255")

        plt.subplot(2,2,3)
        minAmp = 4*(1.0/255)
        dataF1 = 10*np.log10(dataF/minAmp)
        plt.plot(freqAxis, dataF1)
        plt.title("10*log10 wrt 4/255")

        plt.subplot(2,2,4)
        minAmp = 16*(1.0/255)
        dataF1 = 10*np.log10(dataF/minAmp)
        plt.plot(freqAxis, dataF1)
        plt.title("10*log10 wrt 16/255")

        plt.show()
    elif (iDEBUG == 4):
        plt.subplot(2,1,1)
        plt.plot(freqAxis, dataF)
        plt.title("dataF")

        plt.subplot(2,1,2)

        # Because 8bit IQ data, so smallest value sampled is 2/(2**8)
        # ie 1VUnitpeak-peak i.e -1Vunit to +1Vunit or 1/128
        minAmp = (2.0/256)
        dataF = np.clip(dataF, minAmp*0.25, 2.0)
        dataF = 10*np.log10(dataF/minAmp)
        dataF = -gain + dataF
        plt.plot(freqAxis, dataF)
        plt.title("10*log10 wrt 2/256")
        plt.show()



    # Because 8bit IQ data, so smallest value sampled is 2/(2**8)
    # ie 1VUnitpeak-peak i.e -1Vunit to +1Vunit or 1/128
    minAmp = (2.0/256)
    dataF = dataF + minAmp*0.33
    dataF = 10*np.log10(dataF/minAmp)
    dataF = -gain + dataF
    plt.plot(freqAxis, dataF)
    plt.title("10*log10 wrt 2/256")
    plt.show()

    #cairoplot_data(dataF, centerFreq, freqSpan)


handle_args()
sdr = rtlsdr_init()
rtlsdr_info(sdr)


if (gArgs["mode"] == "FULL_SPAN") or (gArgs["mode"] == "SCAN"):
    startFreq = float(gArgs["startFreq"])
    endFreq = float(gArgs["endFreq"])
    sampleRate = float(gArgs["sampleRate"])
    gain = gArgs["gain"]
    print("Mode[{}], startFreq[{}], endFreq[{}], sampleRate[{}], gain[{}]".format(gArgs["mode"], startFreq, endFreq, sampleRate, gain))
    dataF = rtlsdr_scan(sdr, startFreq, endFreq, sampleRate, gain)
    if (gbModeCentered):
        freqArg =  (endFreq+startFreq)/2
    else:
        freqArg = startFreq
    plot_data(None, dataF, freqArg, (endFreq-startFreq), gain)
else:
    centerFreq = float(gArgs["centerFreq"])
    sampleRate = float(gArgs["sampleRate"])
    gain = gArgs["gain"]
    print("Mode[{}], centerFreq[{}], sampleRate[{}], gain[{}]".format(gArgs["mode"], centerFreq, sampleRate, gain))
    if (gbLivePlot):
        rtlsdr_zerospan_repeat(sdr, centerFreq, sampleRate, gain)
    else:
        rtlsdr_setup(sdr, centerFreq, sampleRate, gain)
        data, dataF, dataFDC = rtlsdr_curscan(sdr)
        plot_data(data, dataF, sdr.center_freq, sdr.sample_rate/2, gain)


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
    print("min[{}], max[{}]".format(min(dBytes), max(dBytes)))
    dprint_bytes(10, dBytes)
    dI = np.array(dBytes[::2]) # dBytes[0::2]
    dQ = np.array(dBytes[1::2])
    dMinMax = minmax_iq(dI, dQ)
    return dI, dQ, dMinMax


def read_and_discard(sdr):
    data = sdr.read_samples(4096)
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
    dI, dQ, dMinMax = read_iq(sdr, 2048)
    data = dI + dQ
    print("Data MinMax [{}]".format(dMinMax))

    dataF = np.abs(np.fft.fft(data)/len(data))
    dataF = dataF[:len(dataF)/2]*2
    dataFDC = dataF[0]
    dataF[0] = 0
    dMinMax = minmax(dataF)

    dprint(10, "DataFFT [{}]\n\tLength[{}]\n\tMinMax [{}]\n".format(dataF, len(dataF), dMinMax))
    dprint(2, "DataFFTDC [{}]\n\tLength[{}]\n\tMinMax [{}]\n".format(dataFDC, len(dataF), dMinMax))

    return data, dataF


def rtlsdr_scan(sdr, startFreq, endFreq, sampleRate, gain):
    freqSpan = sampleRate/2
    curFreq = startFreq
    dataFAll = np.array([])
    while curFreq < endFreq:
        rtlsdr_setup(sdr, curFreq, sampleRate, gain)
        data, dataF = rtlsdr_curscan(sdr)
        dataFAll = np.append(dataFAll, dataF)
        cairoplot_data(dataF, curFreq, sampleRate/2)
        curFreq += freqSpan
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
def plot_data(data, dataF, startOrCenterFreq, freqSpan):
    plt.subplot(2,1,1)
    if (data != None):
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
    #freqAxis = np.arange(startFreq, endFreq, deltaFreq)
    freqAxis = np.linspace(startFreq, endFreq, fftBins)
    print("StartFreq[{}], CenterFreq[{}], EndFreq[{}], FreqSpan[{}]".format(startFreq, centerFreq, endFreq, freqSpan))
    print("\tNumOfFFTBins[{}], deltaFreq[{}]".format(fftBins, deltaFreq))
    print("\tNumOf XPoints[{}], YPoints[{}]".format(len(freqAxis), len(dataF)))
    plt.plot(freqAxis, dataF)
    plt.show()

    cairoplot_data(dataF, centerFreq, freqSpan)


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
    plot_data(None, dataF, freqArg, (endFreq-startFreq))
else:
    centerFreq = gArgs["centerFreq"]
    sampleRate = float(gArgs["sampleRate"])
    gain = gArgs["gain"]
    print("Mode[{}], centerFreq[{}], sampleRate[{}], gain[{}]".format(gArgs["mode"], centerFreq, sampleRate, gain))
    rtlsdr_setup(sdr, centerFreq, sampleRate, gain)
    data, dataF = rtlsdr_curscan(sdr)
    plot_data(data, dataF, sdr.center_freq, sdr.sample_rate/2)


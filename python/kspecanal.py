#!/bin/env python3
# kSpecAnal - A Spectrum Analyser
# v20170208_1247, HanishKVC
#

import rtlsdr
import cairo
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

gHeight = 720
gWidth = 1024


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


def print_bytes(dBytes):
    for i in dBytes:
        print(i, end=" ")
    print()


def read_iq(sdr, numBytes):
    dBytes = sdr.read_bytes(numBytes)
    print("min[{}], max[{}]".format(min(dBytes), max(dBytes)))
    print_bytes(dBytes)
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
    return sdr


def rtlsdr_init():
    sdr = rtlsdr.RtlSdr()
    argCnt = len(sys.argv)

    if (argCnt >= 2):
        centerFreq = sys.argv[1]
    else:
        centerFreq = 91.1e6
    if (argCnt >= 3):
        sampleRate = sys.argv[2]
    else:
        sampleRate = 1e6
    #sdr.freq_correction = 0
    if (argCnt >= 4):
        gain = float(sys.argv[3])
    else:
        gain = 0
    rtlsdr_setup(sdr, centerFreq, sampleRate, gain)
    return sdr


def rtlsdr_info(sdr):
    print(dir(sdr))
    print("GainValues/*DivBy10AndThenUse*/:{}".format(sdr.gain_values))
    print("CenterFreq[{}], SamplingRate[{}], Gain[{}]".format(sdr.center_freq, sdr.sample_rate, sdr.gain))


def rtlsdr_curscan(sdr):
    dI, dQ, dMinMax = read_iq(sdr, 2048)
    data = dI + dQ
    print("Data MinMax [{}]".format(dMinMax))

    dataF = np.abs(np.fft.fft(data)/len(data))
    dataF = dataF[:len(dataF)/2]*2
    dMinMax = minmax(dataF)
    print("DataFFT [{}]\n\tLength[{}]\n\tMinMax [{}]\n".format(dataF, len(dataF), dMinMax))
    dataF[0] = 0
    dMinMax = minmax(dataF)
    print("DataFFT MinMax without DCComponent at 0 [{}]".format(dMinMax))

    return data, dataF


def rtlsdr_scan(sdr, startFreq, endFreq, sampleRate):
    gain = sdr.gain
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


bModeCentered = True
def plot_data(data, dataF, centerFreq, freqSpan):
    plt.subplot(2,1,1)
    if (data != None):
        plt.plot(data)
    plt.subplot(2,1,2)
    fftBins = len(dataF)
    deltaFreq = freqSpan/fftBins
    if (bModeCentered):
        startFreq = centerFreq - (freqSpan/2)
        endFreq = centerFreq + (freqSpan/2)
    else:
        startFreq = centerFreq
        endFreq = centerFreq + freqSpan #- deltaFreq
    freqAxis = np.arange(startFreq, endFreq, deltaFreq)
    print("StartFreq[{}], CenterFreq[{}], EndFreq[{}], FreqSpan[{}]".format(startFreq, centerFreq, endFreq, freqSpan))
    print("\tNumOfFFTBins[{}], deltaFreq[{}]".format(fftBins, deltaFreq))
    print("\tNumOf XPoints[{}], YPoints[{}]".format(len(freqAxis), len(dataF)))
    plt.plot(freqAxis, dataF)
    plt.show()

    cairoplot_data(dataF, centerFreq, freqSpan)


sdr = rtlsdr_init()
rtlsdr_info(sdr)

read_and_discard(sdr)
bFullScan = True
if (bFullScan):
    startFreq = 30e6
    endFreq = 1e9
    dataF = rtlsdr_scan(sdr, startFreq, endFreq, 1e6)
    plot_data(None, dataF, (endFreq+startFreq)/2, (endFreq-startFreq))
else:
    data, dataF = rtlsdr_curscan(sdr)
    plot_data(data, dataF, sdr.center_freq, sdr.sample_rate/2)


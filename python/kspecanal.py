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

def plot_setup(dMinMax):
    # Get Cairo Context
    pCRSurface = cairo.SVGSurface("/tmp/t1.svg", gWidth, gHeight)
    pCR = cairo.Context(pCRSurface)
    pCR.set_source_rgb(100,100,200)
    pCR.paint()
    # Calculate the Y Axis details
    pYMin = dMinMax[0]
    pYMax = dMinMax[1]
    pYPerPixel = (pYMax-pYMin)/gHeight
    return (pCRSurface, pCR, pYMin, pYMax, pYPerPixel)

def plot_dot(pSetup, x, y, bRect):
    pCR = pSetup[1]
    if (bRect):
        pCR.rectangle(x,y,x+1,y+1)
        pCR.fill()
        #pCR.stroke()
    else:
        pCR.move_to(x,y)
        pCR.show_text("*")


def plot_xy(pSetup, x, y):
    pCR = pSetup[1]
    pYMin = pSetup[2]
    pYMax = pSetup[3]
    pYPerPixel = pSetup[4]
    cData = y-pYMin
    cY = cData / pYPerPixel
    #print("plot_xy: {},{} = {},{}".format(x,y,x,cY))
    pCR.set_source_rgb(100, 0, 0)
    plot_dot(pSetup, x, cY, False)

def print_bytes(dBytes):
    for i in dBytes:
        print(i, end=" ")
    print()

def read_iq(numBytes):
    dBytes = sdr.read_bytes(numBytes)
    print("min[{}], max[{}]".format(min(dBytes), max(dBytes)))
    print_bytes(dBytes)
    dI = np.array(dBytes[::2]) # dBytes[0::2]
    dQ = np.array(dBytes[1::2])
    dMinMax = minmax_iq(dI, dQ)
    return dI, dQ, dMinMax

def read_and_discard():
    data = sdr.read_samples(4096)
    #print(data)
    dMinMax = minmax_complex(data)
    return data, dMinMax


sdr = rtlsdr.RtlSdr()
argCnt = len(sys.argv)

if (argCnt >= 2):
    sdr.center_freq = sys.argv[1]
else:
    sdr.center_freq = 91.1e6
if (argCnt >= 3):
    sdr.sample_rate = sys.argv[2]
else:
    sdr.sample_rate = 1e6
#sdr.freq_correction = 0
if (argCnt >= 4):
    sdr.gain = float(sys.argv[3])
else:
    sdr.gain = 0
print(dir(sdr))
print("GainValues/*DivBy10AndThenUse*/:{}".format(sdr.gain_values))
print("CenterFreq[{}], SamplingRate[{}], Gain[{}]".format(sdr.center_freq, sdr.sample_rate, sdr.gain))

read_and_discard()

dI, dQ, dMinMax = read_iq(2048)
data = dI + dQ
print("Data MinMax [{}]".format(dMinMax))

dataF = np.abs(np.fft.fft(data)/len(data))
dMinMax = minmax(dataF)
print("DataFFT [{}]\n\tLength[{}]\n\tMinMax [{}]\n".format(dataF, len(dataF), dMinMax))
dataF[0] = 0
dMinMax = minmax(dataF)
print("DataFFT MinMax without DCComponent at 0 [{}]".format(dMinMax))

plt.plot(dataF)
plt.show()

# Scale
dMinMax = (dMinMax[0]*2, dMinMax[1]*2)

pSetup = plot_setup(dMinMax)
print(pSetup)
x=1
for i in dataF:
    plot_xy(pSetup, x, i)
    x+=1

pSetup[0].flush()
pSetup[0].finish()
os.system("display /tmp/t1.svg")


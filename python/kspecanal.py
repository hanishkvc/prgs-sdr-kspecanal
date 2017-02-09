#!/bin/env python3
# kSpecAnal - A Spectrum Analyser
# v20170208_1247, HanishKVC
#

import rtlsdr
import cairo
import os
import sys

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

data = sdr.read_samples(1024)
#print(data)
dMinMax = minmax_complex(data)

dBytes = sdr.read_bytes(64*1024)
dBytes = sdr.read_bytes(1024)
print("min[{}], max[{}]".format(min(dBytes), max(dBytes)))
print_bytes(dBytes)
print(dMinMax)

# Scale
dMinMax = (dMinMax[0]*2, dMinMax[1]*2)

pSetup = plot_setup(dMinMax)
print(pSetup)
x=1
for i in data:
    plot_xy(pSetup, x, i.real)
    x+=1

pSetup[0].flush()
pSetup[0].finish()
os.system("display /tmp/t1.svg")


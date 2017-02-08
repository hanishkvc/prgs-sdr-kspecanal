#!/bin/env python3
# kSpecAnal - A Spectrum Analyser
# v20170208_1247, HanishKVC
#

import rtlsdr
import cairo
import os

gHeight = 720
gWidth = 1024

def minmax(data):
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


sdr = rtlsdr.RtlSdr()

sdr.sample_rate = 1e6
sdr.center_freq = 91.1e6
#sdr.freq_correction = 0
sdr.gain = 0
print(dir(sdr))

data = sdr.read_samples(1024)
print(data)
dMinMax = minmax(data)
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


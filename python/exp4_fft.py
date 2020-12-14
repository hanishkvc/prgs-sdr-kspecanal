#!/usr/bin/env python3
# Test impact of partial signal cycle and windowing
#

import numpy as np
import matplotlib.pyplot as plt


# Setup parameters
sr = 80
f1 = 2
f2 = 8
f3 = 38
startTime = 0
endTime = 2


# Generate time units
t = np.arange(startTime, endTime, 1/sr)
t = np.linspace(startTime, endTime, sr*(endTime-startTime))


# set plot params
nRows = 15
nCols = 4
plt.figure(figsize=(4*nCols, 3*nRows))


def fft_ex(s):
    '''
    Normalised Fft of raw sample (rectangle window)
    Normalised Fft of kaiser windowed sample

    Normalised to reach fft amplitude of 1.0, given input amplitude of 1.0 (to -1.0)
    '''
    fftN = np.abs(np.fft.fft(s)/len(s))*2
    win = np.kaiser(len(s),15)
    sW = s*win
    k2rWin = (len(sW)/win.sum())
    fftWN = np.abs(np.fft.fft(sW)/len(sW))*k2rWin*2
    print('rawFftMax', max(fftN), 'kaisFftMax', max(fftWN))
    return s, fftN, sW, fftWN


def plot_fft(fftN):
    '''
    Plot fft data, till Nyquist freq
    Also match it to its corresponding frequencies
    '''
    fftN = fftN[:int(len(fftN)/2)]
    freqs = np.linspace(0, sr/2, len(fftN))
    plt.grid(True)
    plt.plot(freqs, fftN, ".-")


def plot_it(s, fftN, sW, fftWN, r=nRows, c=nCols, i=1):
    '''
    Plot fft data, till Nyquist freq
    Also match it to its corresponding frequencies
    '''
    print("plot_it", r, c, i)
    plt.subplot(r, c, i)
    plt.plot(s)
    plt.subplot(r, c, i+1)
    plot_fft(fftN)
    plt.subplot(r, c, i+2)
    plt.plot(sW)
    plt.subplot(r, c, i+3)
    plot_fft(fftWN)


def do_it(s, i=1, r=nRows, c=nCols, msg=None):
    if msg != None:
        print("{}:AreaUnderAmps[{}]".format(msg, np.sum(np.abs(s))))
    s,fftN, sW, fftWN = fft_ex(s)
    plot_it(s, fftN, sW, fftWN, r, c, i)
    return fftN, fftWN


def do_it_sliding(s, winSize=0.5, pi=1, r=nRows, c=nCols, msg=None):
    winSamples = int(winSize*sr)
    numLoops = int(len(s)/winSamples)
    fftNCum = np.zeros(winSamples)
    fftWNCum = np.zeros(winSamples)
    for i in range(numLoops):
        sT = s[i*winSamples:(i+1)*winSamples]
        sT,fftN, sW, fftWN = fft_ex(sT)
        fftNCum = fftNCum + np.abs(fftN)
        fftWNCum = fftWNCum + np.abs(fftWN)
    plot_it(s, fftNCum*1/numLoops, sW, fftWNCum*1/numLoops, r, c, pi)
    return fftNCum, fftWNCum


# Generate the signal
s1 = 1*np.sin(2*np.pi*f1*t)
do_it(s1,1, msg="s1")
s2 = 1*np.sin(2*np.pi*f2*t)
do_it(s2,5, msg="s2")
s3 = 2*np.sin(2*np.pi*f3*t)
do_it(s3,9, msg="s3")
s = s1 + s2 + s3
do_it(s,13)
print("NumOfSamples:", len(s))


# Plot the window functions
plt.subplot(nRows, nCols, 17)
plt.plot(np.kaiser(sr, 0))
plt.plot(np.hanning(sr))
plt.plot(np.kaiser(sr, 8))
plt.plot(np.kaiser(sr, 15))


# Help track fftMax's min and max values
dFftMax = { 'raw': {'min': 10, 'max': -10}, 'win': {'min': 10, 'max': -10} }
def fftmax_minmax(dFftMax, fftN, fftWN):
    if dFftMax['raw']['min'] > max(fftN):
        dFftMax['raw']['min'] = max(fftN)
    if dFftMax['raw']['max'] < max(fftN):
        dFftMax['raw']['max'] = max(fftN)
    if dFftMax['win']['min'] > max(fftWN):
        dFftMax['win']['min'] = max(fftWN)
    if dFftMax['win']['max'] < max(fftWN):
        dFftMax['win']['max'] = max(fftWN)


# Check Fft with X.y cycles - ie partial cycle impact
# as the test signal starts at time 0, so also full seconds correspond to full cycles.
# while partial seconds correspond to partial cycle of the signal in the test sample.
startSec = 0
startSec = 0.6
endSec = 1
for i in range(0, 10):
    startIndex = int(sr*startSec)
    endIndex = int(sr*(endSec+(i/20)))
    # Plot Raw signal and fft
    sT = s[startIndex:endIndex]
    print(startIndex, endIndex, len(sT))
    fftN, fftWN = do_it(sT, 21+i*4)
    fftmax_minmax(dFftMax, fftN, fftWN)
print("Raw", dFftMax['raw']['max']/dFftMax['raw']['min'])
print("Win", dFftMax['win']['max']/dFftMax['win']['min'])


# Working with smaller than a second of data
startSec = 0.6
for i in range(0, 10):
    startIndex = int(sr*(startSec+(i/20)))
    sT = s[startIndex:]
    print(startIndex, len(sT))
    fftN, fftWN = do_it_sliding(sT, 0.3, 21+i*4)


# Show the plots
plt.savefig("/tmp/exp4_fft.png")
plt.show()


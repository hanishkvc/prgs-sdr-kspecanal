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


# Generate the signal
s1 = np.sin(2*np.pi*f1*t)
do_it(s1,1, msg="s1")
s2 = np.sin(2*np.pi*f2*t)
do_it(s2,5, msg="s2")
s3 = np.sin(2*np.pi*f3*t)
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
fullSecs = 1
for i in range(0, 10):
    endIndex = int(sr*(fullSecs+(i/20)))
    # Plot Raw signal and fft
    sT = s[0:endIndex]
    print(endIndex, len(sT))
    fftN, fftWN = do_it(sT, 21+i*4)
    fftmax_minmax(dFftMax, fftN, fftWN)
print("Raw", dFftMax['raw']['max']/dFftMax['raw']['min'])
print("Win", dFftMax['win']['max']/dFftMax['win']['min'])


# Show the plots
plt.savefig("/tmp/exp4_fft.png")
plt.show()


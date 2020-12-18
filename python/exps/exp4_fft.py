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
endTime = 8


# Generate time units
t = np.arange(startTime, endTime, 1/sr)
t = np.linspace(startTime, endTime, sr*(endTime-startTime))


# set plot params
pltRows = 35
pltCols = 4
pltWidth = 6*pltCols
pltHeight = 4*pltRows
pltXRes = 0.1/pltWidth
pltYRes = 0.1/pltHeight
plt.figure(figsize=(pltWidth, pltHeight))


def fig_text(text, xD=0, yD=0):
    x0 = plt.gca().get_position().x0
    y0 = plt.gca().get_position().y0
    plt.figtext(x0+xD, y0+yD, text)
    print("******", text)


def pltind(row=0, col=0):
    return row*pltCols + col + 1


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
    #print('rawFftMax', max(fftN), 'kaisFftMax', max(fftWN))
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


def plot_it(s, fftN, sW, fftWN, r=pltRows, c=pltCols, i=1):
    '''
    Plot fft data, till Nyquist freq
    Also match it to its corresponding frequencies
    '''
    print("plot_it", r, c, i)
    plt.subplot(r, c, i)
    plt.plot(s)
    plt.title("NumSamps:{}x".format(len(s)/sr))
    plt.subplot(r, c, i+1)
    plot_fft(fftN)
    plt.subplot(r, c, i+2)
    plt.plot(sW)
    sigSamps = np.count_nonzero(np.abs(sW) < 0.1*np.max(np.abs(sW)))
    plt.title("SignificantSamps:{}x".format(sigSamps/sr))
    plt.subplot(r, c, i+3)
    plot_fft(fftWN)


def do_it(s, i=1, r=pltRows, c=pltCols, msg=None):
    if msg != None:
        print("{}:AreaUnderAmps[{}]".format(msg, np.sum(np.abs(s))))
    s,fftN, sW, fftWN = fft_ex(s)
    plot_it(s, fftN, sW, fftWN, r, c, i)
    return fftN, fftWN


def do_it_sliding(s, winSize=0.5, pi=1, r=pltRows, c=pltCols, msg=None, cType="AVG"):
    winSamples = int(winSize*sr)
    numLoops = int(len(s)/(winSamples*0.1))
    fftNCum = np.zeros(winSamples)
    fftWNCum = np.zeros(winSamples)
    for i in range(numLoops):
        startIndex = int(i*winSamples*0.1)
        endIndex = startIndex + winSamples
        #sT = s[i*winSamples:(i+1)*winSamples]
        sT = s[startIndex:endIndex]
        if len(sT) < winSamples:
            break
        #print("doitslid:", startIndex, endIndex)
        sT,fftN, sW, fftWN = fft_ex(sT)
        if cType == 'AVG':
            fftNCum = (fftNCum + np.abs(fftN))*0.5
        else:
            fftNCum = np.max([fftNCum, np.abs(fftN)], axis=0)
        if cType == 'AVG':
            fftWNCum = (fftWNCum + np.abs(fftWN))*0.5
        else:
            fftWNCum = np.max([fftWNCum, np.abs(fftWN)], axis=0)
    plot_it(s, fftNCum, sW, fftWNCum, r, c, pi)
    return fftNCum, fftWNCum


#### Generate the signal
s1 = 1*np.sin(2*np.pi*f1*t)
do_it(s1, pltind(0), msg="s1")
s2 = 1*np.sin(2*np.pi*f2*t)
do_it(s2, pltind(1), msg="s2")
s3 = 2*np.sin(2*np.pi*f3*t)
do_it(s3, pltind(2), msg="s3")
s = s1 + s2 + s3
do_it(s, pltind(3), msg="All")
print("NumOfSamples:", len(s))


#### Plot the window functions
fig_text("Window Functions", yD=-3*pltYRes)
plt.subplot(pltRows, pltCols, pltind(4))
plt.plot(np.hanning(sr))
for i in range(0,25,2):
    plt.plot(np.kaiser(sr, i))


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



#### Partial Cycles beyond a full sr window
fig_text("Partial Cycles")
# Check Fft with X.y cycles - ie partial cycle impact
# as the test signal starts at time 0, so also full seconds correspond to full cycles.
# while partial seconds correspond to partial cycle of the signal in the test sample.
startSec = 0
startSec = np.random.rand()
endSec = startSec + 1
for i in range(0, 10):
    startIndex = int(sr*startSec)
    endIndex = int(sr*(endSec+(i/20)))
    # Plot Raw signal and fft
    sT = s[startIndex:endIndex]
    print(startIndex, endIndex, len(sT))
    fftN, fftWN = do_it(sT, pltind(5+i))
    fftmax_minmax(dFftMax, fftN, fftWN)
print("Raw", dFftMax['raw']['max']/dFftMax['raw']['min'])
print("Win", dFftMax['win']['max']/dFftMax['win']['min'])



#### For Overlapped sliding
startSec = np.random.rand()


# Overlapped Sliding with partial sr width windows
fig_text("Overlapped sliding - partial window", yD=-3*pltYRes)
for i in range(0, 10):
    startIndex = int(sr*(startSec+(i/10)))
    sT = s[startIndex:]
    print("ForDoItSliding:", startIndex, len(sT))
    fftN, fftWN = do_it_sliding(sT, 0.5, pltind(15+i))


# Overlapped Sliding with full sr width windows
fig_text("Overlapped sliding - full window", yD=-3*pltYRes)
for i in range(0, 10):
    startIndex = int(sr*(startSec+(i/10)))
    sT = s[startIndex:]
    print("ForDoItSliding:", startIndex, len(sT))
    fftN, fftWN = do_it_sliding(sT, 1.0, pltind(25+i))


# Show the plots
#plt.gcf().set_tight_layout(True)
plt.savefig("/tmp/exp4_fft_freq.png")
plt.show()


####
#### Testing signal which is available only for a partial amount of time
####

plt.figure(figsize=(pltWidth, pltHeight))
#plt.gcf().set_tight_layout(True)


def hideandseek(tStart, tEnd, pi):
    iStart = int(tStart*sr)
    iEnd = int(tEnd*sr)
    print("s2", iStart, iEnd)
    s2T = np.copy(s2)
    s2T[iStart:iEnd] = 0
    do_it(s2T, pltind(pi+0), msg="InBtwS2")
    s = s1 + s2T + s3
    do_it(s, pltind(pi+1), msg="All")

    fftN, fftWN = do_it_sliding(s, 0.5, pltind(pi+2))
    fig_text("SlidPartialWindow")
    fftN, fftWN = do_it_sliding(s, 1.0, pltind(pi+3))
    fig_text("SlidFullWindowAvg")
    fftN, fftWN = do_it_sliding(s, 1.0, pltind(pi+4), cType="MAX")
    fig_text("SlidFullWindowMax")


for i in range(4):
    tStart = startSec + i*2
    tEnd = tStart + 3
    hideandseek(tStart, tEnd, i*5)


plt.savefig("/tmp/exp4_fft_time.png")
plt.show()

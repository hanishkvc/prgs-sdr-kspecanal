#!/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

sr = 1024
sr = 10e3
totalTime = 4
totalTime = 3
sampleStartTime = np.random.rand(1)
print("sampleStartTime[{}]".format(sampleStartTime))
t = np.arange(0+sampleStartTime, totalTime+sampleStartTime, 1/sr)

s10k = np.sin(2*np.pi*1e3*t)
s20k = np.sin(2*np.pi*2e3*t)
s17k = np.sin(2*np.pi*1.7e3*t)
s27k = np.sin(2*np.pi*2.7e3*t)
s30k = np.sin(2*np.pi*3e3*t)
s31k = np.sin(2*np.pi*3.1e3*t)
s32k = np.sin(2*np.pi*3.2e3*t)
s33k = np.sin(2*np.pi*3.3e3*t)
s34k = np.sin(2*np.pi*3.4e3*t)
s35k = np.sin(2*np.pi*3.5e3*t)
s36k = np.sin(2*np.pi*3.6e3*t)
s37k = np.sin(2*np.pi*3.7e3*t)
s38k = np.sin(2*np.pi*3.8e3*t)
s39k = np.sin(2*np.pi*3.9e3*t)
s40k = np.sin(2*np.pi*4e3*t)
s47k = np.sin(2*np.pi*4.7e3*t)

s=s10k+s20k+s17k+s27k+s30k+s31k+s32k+s33k+s34k+s35k+s36k+s37k+s38k+s39k+s40k+s47k
s=0
for i in np.arange(0,5e3,0.1e3):
    sT=np.sin(2*np.pi*i*t)
    if (i == 0):
        s = sT
    else:
        s = (s + sT)
        # s = (s + sT)/2 # Averaging in time domain kills many freq components

#DANGER: Note that the comparision of integer with floating point below
# depends on how the language handles such comparison. New versions of
# python seem to handle them better.
def check_multipleOrSubmultiple(mOrS, refNum):
    if (mOrS >= refNum):
        ratio = mOrS/refNum
        ratioInt = int(ratio)
        backFromInt = ratioInt*refNum
        backFromDec = ratio*refNum
        if (ratioInt == ratio):
            return True, "Multiple"
        else:
            return False, "NOT Multiple"
    if (mOrS < refNum):
        ratio = refNum/mOrS
        ratioInt = int(ratio)
        backFromInt = ratioInt*mOrS
        backFromDec = ratio*mOrS
        if (ratioInt == ratio):
            return True, "SubMultiple"
        else:
            return False, "NOT SubMultiple"


def check_fftresult(fftD, dataLen, sr):
    ratio = sr/dataLen
    fftMin = np.min(fftD)
    fftMax = np.max(fftD)
    #fClipped = np.clip(fftD, 0.1, 1.0)
    fRounded = np.round(fftD, 1)
    fftUniqRounded = np.unique(fRounded)
    bMOrS, sType = check_multipleOrSubmultiple(dataLen, sr)
    print("sampleLen[{}], samplingRate[{}], ratio[{}], min[{}], max[{}]\n ratioType[{}], fftUniqRounded[{}]".format(dataLen, sr, ratio, 
        fftMin, fftMax, sType, fftUniqRounded))

def do_fft(s, dataLen, startFreq, endFreq, sr):
    ratio = sr/dataLen
    f = np.zeros(dataLen/2)
    iCnt = max(int(ratio),1)
    for i in range(0, iCnt):
        iStart = i*dataLen
        print("DBG:{}:iStart[{}]".format(i, iStart))
        sT = s[iStart:iStart+dataLen]
        fT=np.fft.fft(sT)
        fT=abs(fT)
        fT=fT/len(fT)
        if (i == 0):
            f = fT[0:len(fT)/2]
        else:
            f=(f+fT[0:len(fT)/2])/2
    fr=np.linspace(startFreq, endFreq,len(f))
    plt.plot(fr,f)
    plt.grid()
    plt.title("sampleLen[{}], ratioWithSR[{}]".format(dataLen, ratio))
    check_fftresult(f, dataLen, sr)

def exp_fft(s, startFreq, endFreq, sr):
    for i in range(1, gLoopCnt):
        plt.subplot(gPlotRows,4,2*i-1)
        plt.plot(s[0:len(s)/i])
        plt.subplot(gPlotRows,4,2*i)
        do_fft(s, int(len(s)/i), startFreq, endFreq, sr)
    # dpi option to savefig didn't seem to help in anyway for ps files
    # Most probably dpi is used for image based formats like svg
    # while papertype seems to help for ps
    plt.savefig("/tmp/exp3_fft.pdf")
    #plt.show()

gLoopCnt = 61
gPlotRows = gLoopCnt/2
plt.figure(num=1, figsize=(24,gPlotRows*4))
exp_fft(s, 0, sr/2, sr)


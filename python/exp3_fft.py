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

s=s10k+s20k+s17k

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

def do_fft(s, startFreq, endFreq, sr):
    f=np.fft.fft(s)
    f=abs(f)
    f=f/len(f)
    f=f[0:len(f)/2]
    fr=np.linspace(startFreq, endFreq,len(f))
    plt.plot(fr,f)
    plt.grid()
    dataLen = len(s)
    ratio = sr/dataLen
    plt.title("sampleLen[{}], ratioWithSR[{}]".format(dataLen, ratio))
    check_fftresult(f, dataLen, sr)

def exp_fft(s, startFreq, endFreq, sr):
    for i in range(1,21):
        plt.subplot(5,4,i)
        do_fft(s[0:len(s)/i], startFreq, endFreq, sr)
    # dpi option to savefig didn't seem to help in anyway for ps files
    # Most probably dpi is used for image based formats like svg
    # while papertype seems to help for ps
    plt.savefig("/tmp/exp3_fft.ps", papertype="a2", orientation="landscape")
    #plt.show()

plt.figure(num=1, figsize=(26,20))
exp_fft(s, 0, sr/2, sr)


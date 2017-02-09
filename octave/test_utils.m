
freq1K = 1000
freq2K = 2000
totalTime = 5
sampRate = 10000
timeDelta=1/sampRate
timeSamples = 0:timeDelta:totalTime-timeDelta;

function data = gen_sin(freq, dBV, sampRate, totalTime)
% data = gen_sin(freq, dBV, sampRate, totalTime)
% freq in Hz
% dBV => 0dBV equals Unity Amplitude, 6dBV equals 2*Unity, -6dBV equals 0.5*Unity
	baseAmp = 10**(dBV/20)
	printf("Genering Sine wave freq:%d dBV:%d sampRate:%d totalTime:%d\n", freq, dBV, sampRate, totalTime)
	timeDelta=1/sampRate
	timeSamples = 0:timeDelta:totalTime-timeDelta;
	data = baseAmp*sin(2*pi*freq*timeSamples);
endfunction

data1K = gen_sin(freq1K, 0, sampRate, totalTime);
data2K = gen_sin(freq2K, 0, sampRate, totalTime);
data1K2K = data1K + data2K;

printf("NOTE: Testing raw and abs of fft result\n")

data1KLen=length(data1K)
d1KfLen=length(fft(data1K))
timeSamplesLen=length(1:sampRate*totalTime)

subplot(2,1,1)
plot(1:sampRate*totalTime, fft(data1K))
subplot(2,1,2)
plot(1:sampRate*totalTime, abs(fft(data1K)))
input("Press Any Key...")

printf("NOTE: Testing need for normalising wrt the amount of data processed\n")
plotRows = 2
plotCols = 4

data1K_1_1000 = data1K(1:1000);
data1K_1_2000 = data1K(1:2000);
data1K_1050_2050 = data1K(1050:2050);
d1Kf_full = abs(fft(data1K));
d1Kf_1_1000 = abs(fft(data1K_1_1000));
d1Kf_1_2000 = abs(fft(data1K_1_2000));
d1Kf_1050_2050 = abs(fft(data1K_1050_2050));

subplot(plotRows, plotCols, 1)
plot(d1Kf_full)
subplot(plotRows, plotCols, 2)
plot(d1Kf_1_1000)
subplot(plotRows, plotCols, 3)
plot(d1Kf_1_2000)
subplot(plotRows, plotCols, 4)
plot(d1Kf_1050_2050)

function r=n1_fft(d)
	r=fft(d)/length(d);
	printf("\nr(1): ")
	r(1)
endfunction

function d_plot(d)
	plot(d)
	printf("\nmin: ")
	min(abs(d))
	printf("\nmax: ")
	max(abs(d))
endfunction

d1Kf_full = abs(n1_fft(data1K));
d1Kf_1_1000 = abs(n1_fft(data1K_1_1000));
d1Kf_1_2000 = abs(n1_fft(data1K_1_2000));
d1Kf_1050_2050 = abs(n1_fft(data1K_1050_2050));

subplot(plotRows, plotCols, 5)
d_plot(d1Kf_full)
subplot(plotRows, plotCols, 6)
d_plot(d1Kf_1_1000)
subplot(plotRows, plotCols, 7)
d_plot(d1Kf_1_2000)
subplot(plotRows, plotCols, 8)
d_plot(d1Kf_1050_2050)

input("Press for next..")

printf("NOTE:Testing the fact of discarding 2nd half of the fft result and multiply by 2")


function r=n2_fft(d)
	printf("fft\n")
	r=fft(d)/length(d);
	r=abs(r);
	mid = length(d)/2;
	%min(r)
	%max(r)
	%r=r(1:mid)+r(end:-1:mid+1);
	r=r(1:mid)*2;
	min(r)
	max(r)
endfunction

plotRows = 3
plotCols = 2

d1Kf = abs(n2_fft(data1K));
d2Kf = abs(n2_fft(data2K));
d1K2Kf = abs(n2_fft(data1K2K));


subplot(plotRows, plotCols, 1)
plot(data1K)
subplot(plotRows, plotCols, 2)
plot(d1Kf)
subplot(plotRows, plotCols, 3)
plot(data2K)
subplot(plotRows, plotCols, 4)
plot(d2Kf)
subplot(plotRows, plotCols, 5)
plot(data1K2K)
subplot(plotRows, plotCols, 6)
plot(d1K2Kf)

input("Press any Key...")

printf("NOTE: Testing out few diff frequencies at different dBVs")

function plot_fft(dF, freqStart, freqStop)
	freqDelta = (freqStop-freqStart)/length(dF)
	freqSamples = freqStart:freqDelta:freqStop-freqDelta;
	freqLen = length(freqSamples)
	fftLen = length(dF)
	plot(freqSamples,dF)
endfunction

sampRate = 10000
d1 = gen_sin(1000,  0, sampRate, 2);
d2 = gen_sin(2000, -6, sampRate, 2);
d4 = gen_sin(4000,  0, sampRate, 2);
d6 = gen_sin(6000,  0, sampRate, 2);
dAll = d1 + d2 + d4 + d6;
dSafe = d1 + d2 + d4;
d1f = n2_fft(d1);
d2f = n2_fft(d2);
d4f = n2_fft(d4);
d6f = n2_fft(d6);
dAllf = n2_fft(dAll);
dSafef = n2_fft(dSafe);

pR = 3
pC = 4

subplot(pR, pC, 1)
plot(d1(1:1024))
title("1KHz")
subplot(pR, pC, 2)
plot(d1f)
title("1KHz fft")

subplot(pR, pC, 3)
plot(d2(1:1024))
title("2KHz")
subplot(pR, pC, 4)
plot(d2f)
title("2KHz fft")

subplot(pR, pC, 5)
plot(d4(1:1024))
title("4KHz")
subplot(pR, pC, 6)
title("4KHz fft")
plot(d4f)

subplot(pR, pC, 7)
plot(d6(1:1024))
title("6KHz")
subplot(pR, pC, 8)
plot(d6f)
plot_fft(d6f, 0, sampRate/2)
title("6KHz plot__fft")

subplot(pR, pC, 9)
plot(dAll(1:1024))
title("1+2+4+6")
subplot(pR, pC, 10)
plot(dAllf)
title("1+2+4+6 fft")

subplot(pR, pC, 11)
plot(dSafe(1:1024))
title("1+2+4")
subplot(pR, pC, 12)
plot_fft(dSafef, 0, sampRate/2)
title("1+2+4 plot__fft")

% 20/10
% 0:20/10:20



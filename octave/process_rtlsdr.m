
function process_rtlsdr(iMode, varargin)
% process_rtlsdr(iMode, rtl-sdr-samples1.bin, rtl-sdr-samples2.bin, ...)
% v20170209, HanishKVC
% iMode =  0, Use the full data
% iMode =  1, Use the 1st 2048 samples from the data
% iMode >= 2, Group the data in the samples into 2048 sample groups and
%             Decimate between these groups
% the sample files are saved from rtl-sdr
	numOfFiles = length(varargin)
	numOfPlotCols = 7
	for i = 1:numOfFiles
		printf("%d=%s\n", i, varargin{i})
		dRaw = load_rtlsdr(varargin{i});
		dLen = length(dRaw)
		dDecimatedLen = 2048
		dNumOfGroups = dLen/dDecimatedLen
		if (iMode == 0)
			dComp = dRaw;
		elseif (iMode == 1)
			dComp = dRaw(1:dDecimatedLen);
		else
			dTemp = reshape(dRaw, dDecimatedLen, dNumOfGroups);
			dComp = sum(dTemp, 2);
		endif

		dReal = real(dComp);
		dRealFft = n_fft(dReal);
		dImag = imag(dComp);
		dImagFft = n_fft(dImag);
		dAbs = abs(dComp);
		dAbsFft = n_fft(dAbs);
		dSum = dReal+dImag;
		dSumFft = n_fft(dSum);

		indexStart = numOfPlotCols*(i-1)
		subplot(numOfFiles, numOfPlotCols, indexStart+1)
		plot(dReal)
		subplot(numOfFiles, numOfPlotCols, indexStart+2)
		plot(dImag)
		subplot(numOfFiles, numOfPlotCols, indexStart+3)
		plot(dComp)
		subplot(numOfFiles, numOfPlotCols, indexStart+4)
		plot(abs(dRealFft))
		subplot(numOfFiles, numOfPlotCols, indexStart+5)
		plot(abs(dImagFft))
		subplot(numOfFiles, numOfPlotCols, indexStart+6)
		plot(abs(dAbsFft))
		subplot(numOfFiles, numOfPlotCols, indexStart+7)
		plot(abs(dSumFft))

		title(varargin{i})
	endfor
endfunction

function r=n_fft(d)
% r=n_fft(d)
% normalised fft
	r = fft(d)/length(d);
	% r(1), r(2), r(end)
endfunction


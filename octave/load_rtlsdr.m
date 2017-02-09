function d = load_rtlsdr(filename)
% d = load_rtlsdr(filename)
% v20170209, HanishKVC
% reads uint8 IQ samples from the file saved by rtl_sdr
% and converts them to complex samples
%

	fData = fopen(filename, 'rb');
	d = fread(fData, 'uint8=>double');

	d = d-127;
	d = d(1:2:end) + i*d(2:2:end);
endfunction


function f0_ceps=ceps(x,fs,len)
ms1=fs/500;                 % maximum speech Fx at 1000Hz
ms20=fs/50;                 % minimum speech Fx at 50Hz
% plot waveform
% do fourier transform of windowed signal
Y=fft(x.*hamming(length(x)));
% Y=fft(x);

% plot spectrum of bottom 5000Hz
% hz5000=5000*length(Y)/fs;
% f=(0:hz5000)*fs/length(Y);
% cepstrum is DFT of log spectrum
C=fft(log(abs(Y)+eps));
% plot between 1ms (=1000Hz) and 20ms (=50Hz)
if ms20>len
    ms20=len;
end
q=(ms1:ms20)/fs;
[c,fx]=max(abs(C(ms1:ms20)));
f0_ceps=fs/(ms1+fx-1);
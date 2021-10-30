function f0_acorr=acorr(x,fs)

ms2=fs/500;                % maximum speech Fx at 500Hz
ms20=fs/50;                % minimum speech Fx at 50Hz
r=xcorr(x,ms20,'coeff');   
% just look at region corresponding to positive delays
r=r(ms20+1:2*ms20+1);
[rmax,tx]=max(r(ms2:ms20));
f0_acorr=fs/(ms2+tx-1);
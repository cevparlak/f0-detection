% calculates lpc, cepstrum, yin, yaapt, and new formants
% clc; clear all; close all;
% % [filename, filepath] = uigetfile('*.wav','Select wave file');
% % [upperPath, deepestFolder, ~] = fileparts(filepath(1:length(filepath)-1));
filepath='E:\PhD\datasets\texas\genderx\';
% filepath='E:\PhD\datasets\texas\genderxbp400_3400\';
% filepath='C:\ex\texas\gender\';
currentdir = filepath;
dir1 = currentdir;
params=[];
fnames = dir(fullfile([dir1 '*.wav']));
size1=size(fnames,1);
data_x={};
data_s={};
power_x={};
% return
alpha1=0.97;
startpos=1;
% % for index=1:size1
% % %     fullfilename = [dir1 namesall{index}]; % for timit
% % %     return
% %     fullfilename = [dir1 fnames(index).name];  % for nontimit
% %     [x fs] = audioread(fullfilename);
% % %     x=[x;x];
% %     l1=size(x,1);
% %     p1=nextpow2(l1);
% %     
% %     l2=2.^p1;
% %     if l2<4096
% %         l2=4096;
% %     end
% %     len=l2;
% %     x=padarray(x,l2-l1,'post');
% % %     s = filter( [1 -alpha11], 1, x ); % pre-emphasize    
% % %     x=bandpass(x,[400 3400],fs);
% % %     x=bandpass(x,[400 3400],fs);        
% % %     s=bandpass(s,[400 3400],fs);
% % %     s=bandpass(s,[400 3400],fs);    
% %     data_x{index}=x;
% % %     data_s{index}=s;
% % len=-1;
% % len=1024;
% % %     s=data_s{index};
% %     if len==-1
% %     len=length(x);
% %     else
% %     x = x(1:len);
% % %     s = s(1:len);
% %     end
% % %     return
% %     res = fs/len;
% % %     t=(0:len-1)/fs;       % times of sampling instants
% %     [power1 freqs] = runfft(x,fs,startpos,len);
% %     power1 = power1*1000;
% %     power_x{index}=power1;
% % if rem(index,100)==0
% %     index
% % end 
% % end

% return

params=[];
e=[];
l=[];
e_min=1000000000;
ll_min=1000000000;
best=[];
load power_tex.mat
load data_x_tex.mat
% load power_tex_bp.mat
% load data_x_tex_bp.mat

load fnames_tex.mat
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii9=1:1
for ii8=1:1
for ii7=1:1
for ii6=1:1
for ii5=1:1
for ii4=1:10
for ii3=1:10
for ii2=1:10
for ii1=1:1
for ii=1:1

rx1=12;
rx2=5;
minf0=65;
hi_freq=1200; % start of high frequency region
ratio1=50; % minimum ratio for frequencies below hi_freq

low_freq1=1900; % this is not used
low_freq2=500; % this is not used
ratio2=10; % unused minimum ratio for frequencies above hi_freq 
max_formant_count=13; % unused
bank_bwidth=150;    
    
% max_formant_count=13;
% minf0=65;
% hi_freq=1300; % start of high frequency region
% ratio1=40; % minimum ratio for frequencies below hi_freq
% ratio2=25; % minimum ratio for frequencies above hi_freq
% low_freq1=1900; % this is not used
% low_freq2=500; % this is not used
% bank_bwidth=150;

    
% % len = power(2,ii9+7);
% max_formant_count=4+ii1;
% minf0=45+ii2*5;
% hi_freq=900+ii3*100; % start of high frequency region
% ratio1=5+ii4*5; % minimum ratio for frequencies below hi_freq
% ratio2=5+ii5*5; % minimum ratio for frequencies above hi_freq
% % we delete formants which are less than 1/10th of the highest
% % formant below hi_freq
% % and formants which are 1/20th of the highest formant
% % above hi_freq
% bank_bwidth=150;
% equ_loudness=1; % if 1 apply equal loudness
% low_freq1=1900; % this is not used
% low_freq2=500; % this is not used   
    
% alpha11=0.97;
% tex for texas vowel dataset
% hil for hillenbrand vowel dataset
% to save the file as .mat, must be 3 char long
% allfeatures=[];
allf0=[];
h1=0;
h2=0;
h3=0;
h4=0;
h5=0;
tt=[];
for index=1:size(fnames,1)
len=-1;
len=1024;
    power1=power_x{index};
    fs=16000;
    res = fs/len;        
%     x=data_x{index};
%     x_yaapt=x(1:2048); % problem when equal to 1024
% %     s=data_s{index};
%     if len==-1
%     len=length(x);
%     else    
%     x = x(1:len);
% %     s = s(1:len);    
%     end
% %     res = fs/len;
%     t=(0:length(x)-1)/fs;       % times of sampling instants
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % get lpc formants
% %     [formants f h]=get_lpc_formants(s,fs,len);
%     [power1 freqs] = runfft(x,fs,startpos,len);
%     power1 = power1*1000;
%     [power2 freqs2] = runfft(s,fs,startpos,len);
%     power2 = power2*1000;

%     logpower=20*log10(abs(power1*1000));
%     logpower2=20*log10(abs(power2*1000));
%     powerx=power1;
%     power2x=power2;
%     logpowerx=logpower;
%     logpower2x=logpower2;
%      winsize=4; % window size for smoothing
%         power=smooth_fft(power,winsize,len);
%         power2=smooth_fft(power2,winsize,len);
%     Lpclist=formants;
%     Lpclist=sortrows(Lpclist,-2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [cepstrum cepsfreqs]=get_cepstrum_formants(s,fs,15);
%     if (size(cepstrum,1)<max_formant_count)
%          cepstrum=padarray(cepstrum,max_formant_count-size(cepstrum,1),'post');
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get formants
% equ_loudness=1;
% t1=tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[f0_1 amp0_1 f0_1x amp0_1x]=get_f0_2(power1,res,len,minf0,hi_freq,ratio1,ratio2,max_formant_count,rx1,rx2);
% t2=toc(t1);
% tt=[tt;t2];
% h1=h1+t2;
% return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=data_x{index};
% x_yaapt=x(1:2048); % problem when equal to 1024
x=x(1:len);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% autocorrelation
t1a=tic;
ham = hamming(len); 
xx = x .* ham;  
f0_acorr=acorr(xx,fs);
t2a=toc(t1a);
tt=[tt; t2a];
h2=h2+t2a;
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % cepstrum
t1b=tic;
f0_ceps=ceps(x,fs,len);
t2b=toc(t1b);
tt=[tt; t2b];
h3=h3+t2b;
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % YIN
% t1c=tic;
% f0yin=YIN(x,fs);
% t2c=toc(t1c);
% tt=[tt;t2c];
% h4=h4+t2c;
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % fullfilename
% % [f0yaapt1, nf1] = yaapt(x_yaapt, fs, 0, [], 0, 1);
% % [f0yaapt2, nf2] = yaapt(x_yaapt, fs, 0, [], 0, 2);
% % ham2 = hamming(2048); 
% xxx=padarray(x,1024,'post');
% x_yaapt = xxx;%.*ham2; 
% t1d=tic;
% [f0yaapt3, nf3] = yaapt(x_yaapt, fs, 0, [], 0, 3);
% t2d=toc(t1d);
% tt=[tt;t2d];
% h5=h5+t2d;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     allf0=[allf0; f0yaapt1(1,1) f0yaapt2(1,1) f0yaapt3(1,1)];
%     allf0=[allf0; f0_1 f0_1x f0_2 f0_2x f0_acorr f0_ceps f0yin(1,1) f0yaapt3(1,1)];
%     allf0=[allf0;f0_1 f0_1x f0_2 f0_2x f0_3 f0_3x f0_4 f0_4x f0_acorr];
%     allf0=[allf0;f0_acorr f0_ceps];
%     allf0=[allf0;f0_1 f0_1x];
%     allf0=[allf0; f0_1 f0_1x f0_acorr f0_ceps f0yin(1,1) f0yaapt3(1,1) ];
    allf0=[allf0; f0_1 f0_1x f0_acorr f0_ceps];
%     allf0=[allf0; f0_1 f0_1x];

end

load texas_gt_names3.mat  % ground truth for files in gender folder
gt1= cell2mat(gt(:, 2));
% a=a(:,[1 3 5 7 9 11 13 15 17 18 19 20 21 22 23]);
allf01=allf0(:,[1:size(allf0,2)]);
% a1=allf0(:,[1 3 5 7 9 11 13 15 16 17 18 19 20 21]);
gt2=repmat(gt1,[1, size(allf01,2)]);
d=abs(allf01-gt2);  % error
s=sum(d);     % sum of error
d1=[d;s]; 

list2=[];
s2=[];
ss1=size(allf01,2);
zz1=zeros(1,4); % used to separate cycles f,k,m + sum(f,k,m)
% find error in genders
ff=gt1(1:1110,1);
kk=gt1(1111:2082,1);
mm=gt1(2083:3314,1);

for index=1:ss1
f=allf01(1:1110,index);
k=allf01(1111:2082,index);
m=allf01(2083:3314,index);

f1=f((f<150) | (f>300),:);
k1=k((k<150) | (k>400),:);
m1=m((m< 50) | (m>200),:);
size2=[size(f1,1) size(k1,1) size(m1,1)];

f1a=f((abs(f-ff) > ff/10),:);
k1a=k((abs(k-kk) > kk/10),:);
m1a=m((abs(m-mm) > mm/10),:);

size2a=[size(f1a,1) size(k1a,1) size(m1a,1)] ;
list1=[size2 sum(size2); size2a sum(size2a)];
list2=[list2;list1];
end

l=[l;list2;zz1];
e=[e;s];
params=[params; alpha1 len max_formant_count minf0 hi_freq bank_bwidth low_freq1 low_freq2 ratio1 ratio2];
ee=min(s);
l2=l;
l2( ~any(l2,2), : ) = [];  %rows
% data( :, ~any(data,1) ) = [];  %columns
[ll idx]=min(l2(:,end));
% l2(idx,:)
if ee<=e_min
    e_min=ee;
    best=[best; alpha1 len max_formant_count minf0 hi_freq bank_bwidth low_freq1 low_freq2 ratio1 ratio2 ii ii1 ii2 ii3 ii4 ii5 ii6 ii7 ii8 ii9 l2(idx,:) s];
end

if ll<=ll_min
    ll_min=ll;
    best=[best; alpha1 len max_formant_count minf0 hi_freq bank_bwidth low_freq1 low_freq2 ratio1 ratio2 ii ii1 ii2 ii3 ii4 ii5 ii6 ii7 ii8 ii9 l2(idx,:) s];
end
% return

e_min
ll_min
disp( [num2str(ii) '  ' num2str(ii1) '  ' num2str(ii2) '  ' num2str(ii3) '  ' num2str(ii4) '  ' num2str(ii5) '  ' num2str(ii9)])

% calculates e10
% first column of allf0 is ground truth
e10list=[gt allf0];
e2=[];
for j=2:size(e10list,2)
e1=0;
for i=1:size(e10list,1)
    if abs(e10list(i,j)-e10list(i,1))>e10list(i,1)/10
        e1=e1+1;
    end

end
e2=[e2;e1];    
end
e3=e2/size(e10list,1);
e2=[e2 e3];

return
end
end
end
end


end
end
end
end
end
end
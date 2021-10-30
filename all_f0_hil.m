% clc; clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dir1='E:\PhD\datasets\hillenbrand\allx\'; % directory of wave files
fnames = dir(fullfile([dir1 '*.wav']));
size1=size(fnames,1);
% data_x={};
% data_s={};
% power_x={};
% 
len = 1024;
alpha1=0.97;
startpos=1;
for index=1:size1
    fullfilename = [dir1 fnames(index).name];  
    [x fs] = audioread(fullfilename);
    l1=size(x,1);
    p1=nextpow2(l1);
    l2=2.^p1;
    if l2<4096
        l2=4096;
    end
    len=l2;    
    x=padarray(x,l2-l1,'post');

%     s = filter( [1   -alpha1], 1, x ); % pre-emphasize    
%     x=bandpass(x,[400 3400],fs);
%     x=bandpass(x,[400 3400],fs);        
%     s=bandpass(s,[400 3400],fs);
%     s=bandpass(s,[400 3400],fs);        
    data_x{index}=x;
%     data_s{index}=s;

% len=-1;
len=1024;
x=x(1:len);
% audiowrite([dir3 fnames(index).name],x,fs);
% x = x .* hamming(len);
% audiowrite([dir2 fnames(index).name],x,fs);
% [dir2 '00_' fnames(index).name]
% %     s=data_s{index};
%     if len==-1
%     len=length(x);
%     else
%     x = x(1:len);
% %     s = s(1:len);
%     end
%     return
    res = fs/len;
%     t=(0:len-1)/fs;       % times of sampling instants
    [power1 freqs] = runfft(x,fs,startpos,len);
%     return
    power1 = power1*1000;
    power_x{index}=power1;

if rem(index,100)==0
    index
end    

end
% load power_hil.mat
% load data_x_hil.mat
% load power_hil_bp.mat
% load data_x_hil_bp.mat
 
% load fnames_hil.mat
params=[];
e=[];
l=[];
e_min=1000000000;
ll_min=1000000000;
best=[];
% % jj1 and jj2 may not be needed
% for jj1=1:1
%     rx1=jj1;
%     for jj2=1:1
%         rx2=jj2;

rx1=12;
rx2=5;
minf0=65;
hi_freq=1200; % start of high frequency region
ratio1=50;    % minimum ratio for frequencies below hi_freq
% rx1=8;
% rx2=2;
% minf0=75;
% hi_freq=1200; % start of high frequency region
% ratio1=30;    % minimum ratio for frequencies below hi_freq
low_freq1=1900; % this is not used
low_freq2=500; % this is not used
ratio2=10; % unused minimum ratio for frequencies above hi_freq 
max_formant_count=13; % unused
bank_bwidth=150;         

allf0=[];
h1=0;
h2=0;
h3=0;
h4=0;
h5=0;
t1=tic;
tt=[];
tt2=[];
tt3=[];
for index=1:size(fnames,1)
len=-1;
len=1024;
    power1=power_x{index};
    fs=16000;
    res=fs/len;    
    
% % % % % % %     s=data_s{index};
% % % % %     if len==-1
% % % % %     len=length(x);
% % % % %     else
% % % % %     x = x(1:len);
% % % % % %     s = s(1:len);
% % % % %     end
% % % % %     res = fs/len;
% % % % % %     t=(0:len-1)/fs;       % times of sampling instants
% % % % % t1=tic;
% % % % % 
% % % % %     [power1 freqs] = runfft(x,fs,startpos,len);
% % % % %     power1 = power1*1000;
% % % %     power_x{index}=power1;

% % % %     return
% %     res = fs/len;
% %     t=(0:len-1)/fs;       % times of sampling instants
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get formants
% equ_loudness=1;

x=data_x{index};
% % % x_yaapt=x(1:2048); % problem when equal to 1024
x=x(1:len);

% t1=tic;
% % ham = hamming(len); 
% % xx = x.* ham;  
[f0_1 amp0_1 f0_1x amp0_1x]=hdm(power1,res,len,minf0,hi_freq,ratio1,ratio2,max_formant_count,rx1,rx2);
% t2=toc(t1);
% h1=h1+t2;
% tt=[tt;t2];
% return
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % autocorrelation
t1a=tic;
% ham = hamming(len); 
xx = x.*hamming(len);  
f0_acorr=acorr(xx,fs);
t2a=toc(t1a);
h2=h2+t2a;
tt2=[tt2; t2a];
% % % 
% % % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cepstrum
t1b=tic;
f0_ceps=ceps(x,fs,len);
t2b=toc(t1b);
h3=h3+t2b;
tt3=[tt3; t2b];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % YIN
% t1c=tic;
% f0yin=YIN(xx,fs);
% t2c=toc(t1c);
% tt=[tt;t2c];
% h4=h4+t2c;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% yaapt
% [f0yaapt1, nf1] = yaapt(x_yaapt, fs, 0, [], 0, 1);
% [f0yaapt2, nf2] = yaapt(x_yaapt, fs, 0, [], 0, 2);
% ham2 = hamming(2048); 
% xxx=padarray(x,1024,'post');
% % x_yaapt = xxx.*ham2; 
% x_yaapt = xxx; 
% t1d=tic;
% [f0yaapt3, nf1] = yaapt2016(x_yaapt, fs, 0, [], 0, 3);
% t2d=toc(t1d);
% tt=[tt;t2d];
% h5=h5+t2d;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     allf0=[allf0; f0_1 f0_1x f0_acorr f0_ceps f0yin(1,1) f0yaapt3(1,1) ];
    allf0=[allf0; f0_1 f0_1x f0_acorr f0_ceps ];
    
end

load hillenb_gt.mat  % ground truth of files in hillenb\all folder
allf01=allf0(:,1:size(allf0,2));
gt1=repmat(gt,[1, size(allf01,2)]);
d=abs(allf01-gt1);  % error
s=sum(d);     % sum of error
d1=[d;s]; 

list2=[];
s2=[];
ss1=size(allf01,2);
zz1=zeros(1,5);
% find error in genders
bb=gt1(1:324,1);
gg=gt1(325:552,1);
mm=gt1(553:1092,1);
ww=gt1(1093:1668,1);

for index=1:ss1
b=allf01(1:324,index);
g=allf01(325:552,index);
m=allf01(553:1092,index);
w=allf01(1093:1668,index);

% gender boundaries
b1=b((b<150) | (b>400),:);
g1=g((g<150) | (g>400),:);
m1=m((m<50) | (m>200),:);
w1=w((w<150) | (w>300),:);
size2=[size(b1,1) size(g1,1) size(m1,1) size(w1,1)];

% b1a=b((b<min(b)) | (b>max(b)),:);
% g1a=g((g<min(g)) | (g>max(g)),:);
% m1a=m((m<min(m)) | (m>max(m)),:);
% w1a=w((w<min(w)) | (w>max(w)),:);
% size2a=[size(b1a,1) size(g1a,1) size(m1a,1) size(w1a,1)];

b1a=b((abs(b-bb) > bb/10),:);
g1a=g((abs(g-gg) > gg/10),:);
m1a=m((abs(m-mm) > mm/10),:);
w1a=w((abs(w-ww) > ww/10),:);
size2a=[size(b1a,1) size(g1a,1) size(m1a,1) size(w1a,1)];

list1=[size2 sum(size2);size2a sum(size2a)];
list2=[list2;list1];
end

l=[l;list2;zz1];
e=[e;s];
params=[params; alpha1 len max_formant_count minf0 hi_freq bank_bwidth low_freq1 low_freq2 ratio1 ratio2];

ee=min(s);
l2=l;
l2( ~any(l2,2), : ) = [];  % rows
% data( :, ~any(data,1) ) = [];  % columns
[ll idx]=min(l2(:,end));

l2(idx,:);
if ee<=e_min
    e_min=ee;
%     best=[best; alpha1 len max_formant_count minf0 hi_freq bank_bwidth low_freq1 low_freq2 ratio1 ratio2 ii ii1 ii2 ii3 ii4 ii5 ii6 ii7 ii8 ii9 jj1 jj2 l2(idx,:) s];
end

if ll<=ll_min
    ll_min=ll;
%     best=[best; alpha1 len max_formant_count minf0 hi_freq bank_bwidth low_freq1 low_freq2 ratio1 ratio2 ii ii1 ii2 ii3 ii4 ii5 ii6 ii7 ii8 ii9 jj1 jj2 l2(idx,:) s];
end

e_min
ll_min
% disp( [num2str(ii) '  ' num2str(ii1) '  ' num2str(ii2) '  ' num2str(ii3) '  ' num2str(ii4) '  ' num2str(ii5) '  ' num2str(ii9)])
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

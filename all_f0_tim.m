% clc; clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [filename, filepath] = uigetfile('*.wav','Select wave file');
% [upperPath, deepestFolder, ~] = fileparts(filepath(1:length(filepath)-1));
dir1='E:\PhD\datasets\timitcorpus\timit_LDC93S1\phn\allex\';
% filepath='E:\PhD\datasets\timitcorpus\timit_LDC93S1\phn\allexbp\';
fnames = dir(fullfile([dir1 '*.wav']));
size1=size(fnames,1);
data_x={};
data_s={};
power_x={};
alpha1=0.97;
startpos=1;
% put the names in F, M order for timit
names1={};
names2={};
i1=0;
i2=0;
for i=1:size1
    str1=fnames(i).name;
    if str1(4)=='F'
        i1=i1+1;
        names1{i1}=str1;

    end
    
    if str1(4)=='M'
        i2=i2+1;
        names2{i2}=str1;
    end
end
names1=names1';
names2=names2';
namesall=[names1;names2];
xx=zeros(1024,1);
for index=1:size1
    fullfilename = [dir1 namesall{index}]; % for timit
%     fullfilename = [dir1 fnames(index).name];  %for nontimit    
%     return
    [x fs] = audioread(fullfilename);
%     return
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
len=-1;
len=1024;
x=x(1:len);
% audiowrite([dir3 fnames(index).name],x,fs);
% x = x .* hamming(len);
% audiowrite([dir2 fnames(index).name],x,fs);

%     s=data_s{index};
    if len==-1
    len=length(x);
    else
    x = x(1:len);
%     s = s(1:len);
    end
%     return
    res = fs/len;
%     t=(0:len-1)/fs;       % times of sampling instants
    [power1 freqs] = runfft(x,fs,startpos,len);
    power1 = power1*1000;
    power_x{index}=power1;
if rem(index,5000)==0
    index
end
end
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params=[];
e=[];
l=[];
e_min=1000000000;
ll_min=1000000000;
best=[];
% load power_tim.mat
% load data_x_tim.mat
% load power_tim_bp.mat
% load data_x_tim_bp.mat

% load fnames_tim.mat
% load namesall_tim

% rx1=8;
% rx2=2;
% minf0=75;
% hi_freq=1200; % start of high frequency region
% ratio1=30; % minimum ratio for frequencies below hi_freq
     
rx1=12;
rx2=5;
minf0=65;
hi_freq=1200; % start of high frequency region
ratio1=50;    % minimum ratio for frequencies below hi_freq

% minf0=45+ii1*5;
% hi_freq=900+ii2*100; % start of high frequency region
% ratio1=0+ii3*5; % minimum ratio for frequencies below hi_freq

low_freq1=1900; % this is not used
low_freq2=500; % this is not used
ratio2=10; % unused minimum ratio for frequencies above hi_freq 
max_formant_count=10; % unused
bank_bwidth=150;   

% alpha1=0.97;
% allfeatures=[];
allf0=[];
h1=0;
h2=0;
h3=0;
h4=0;
h5=0;
tt=[];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% return
for index=1:size(fnames,1)
len=-1;
len=1024;    
% power1=power_x{index};
% fs=16000;
% res = fs/len; 
    
%     len=-1;
%     len=1024;
%     x=data_x{index};
%     x_yaapt=x(1:2048); % problem when equal to 1024
% %    
% %     s=data_s{index};
%     if len==-1
%     len=length(x);
%     else    
%     x = x(1:len);
% %     s = s(1:len);    
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[f0_1 amp0_1 f0_1x amp0_1x]=hdm(power1,res,len,minf0,hi_freq,ratio1,ratio2,max_formant_count,rx1,rx2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=data_x{index};
% % % x_yaapt=x(1:2048); % problem when equal to 1024
x=x(1:len);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % autocorrelation
t1a=tic;
ham = hamming(len); 
xx = x .* ham;  
f0_acorr=acorr(xx,fs);
t2a=toc(t1a);
tt=[tt; t2a];
h2=h2+t2a;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cepstrum
t1b=tic;
f0_ceps=ceps(x,fs,len);
if a(index,1)==0
    f0_ceps=0;
end
t2b=toc(t1b);
tt=[tt; t2b];
h3=h3+t2b;
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YIN
% t1c=tic;
% f0yin=YIN(x,fs);
% t2c=toc(t1c);
% tt=[tt;t2c];
% h4=h4+t2c;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fullfilename
% [f0yaapt1, nf1] = yaapt(x_yaapt, fs, 0, [], 0, 1);
% [f0yaapt2, nf2] = yaapt(x_yaapt, fs, 0, [], 0, 2);
% ham2 = hamming(2048); 
% xxx=padarray(x,1024,'post');
% x_yaapt = xxx.*ham2; 
% x_yaapt = xxx; 
% t1d=tic;
% [f0yaapt3, nf3] = yaapt2016(x_yaapt, fs, 0, [], 0, 3);
% t2d=toc(t1d);
% tt=[tt;t2d];
% h5=h5+t2d;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     allf0=[allf0; f0_1 f0_1x f0_acorr f0_ceps f0yin(1,1) f0yaapt3(1,1) ];
    allf0=[allf0; f0_1 f0_1x f0_acorr f0_ceps];
if rem(index, 2500)==0
    index
end

end

% return
load timit_gt7.mat  % 
a1=allf0;
a=repmat(gt,[1, size(a1,2)]);
d=abs(a1-a);  % error
s=sum(d);     % sum of error
d1=[d;s]; 

list2=[];
s2=[];
ss1=size(a1,2);
zz1=zeros(1,3);
% find error in genders
ff=a(1:24017,1);
mm=a(24018:78374,1);
% gt=ones(16066,1)*220;
% gt=[gt;ones(52851-16066,1)*120];
for index=1:ss1
f=a1(1:24017,index);
m=a1(24018:78374,index);
% gender boundaries
m1=m((m<50) | (m>200),:);
f1=f((f<150) | (f>300),:);
size2=[size(m1,1) size(f1,1)];

% b1a=b((b<min(b)) | (b>max(b)),:);
% g1a=g((g<min(g)) | (g>max(g)),:);
% m1a=m((m<min(m)) | (m>max(m)),:);
% w1a=w((w<min(w)) | (w>max(w)),:);
% size2a=[size(b1a,1) size(g1a,1) size(m1a,1) size(w1a,1)];

m1a=m((abs(m-mm) > mm/10),:);
f1a=f((abs(f-ff) > ff/10),:);
size2a=[size(m1a,1) size(f1a,1)];

list1=[size2 sum(size2);size2a sum(size2a)];
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

l2(idx,:);
if ee<=e_min
    e_min=ee;
%     best=[best; alpha1 len max_formant_count minf0 hi_freq bank_bwidth low_freq1 low_freq2 ratio1 ratio2 ii ii1 ii2 ii3 ii4 ii5 ii6 ii7 ii8 ii9 l2(idx,:) s];
end

if ll<=ll_min
    ll_min=ll;
%     best=[best; alpha1 len max_formant_count minf0 hi_freq bank_bwidth low_freq1 low_freq2 ratio1 ratio2 ii ii1 ii2 ii3 ii4 ii5 ii6 ii7 ii8 ii9 l2(idx,:) s];
end

e_min
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


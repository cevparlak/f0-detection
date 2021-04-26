% s2=[s(1:1000);s(30000:31000)];
% clc;clear all;
% close all;
figure(1)
plot(e1(:,2),'red')
hold on;
plot(e1(:,1),'blue')
% gscatter(e1(:,1),e1(:,2));
xlabel('s a m p l e   n u m b e r','FontSize',10,'FontWeight','bold');
ylabel('f_0','FontSize',14,'FontWeight','bold','Rotation',0, 'VerticalAlignment','top', 'HorizontalAlignment','right');
title('Hillenbrand Dataset');
txt = '\uparrow b  o  y  s';
text(10,170,txt)
txt = '\uparrow g  i  r  l  s';
text(310,170,txt)
txt = '\uparrow m  a  n';
text(550,70,txt)
txt = '\uparrow w  o  m  a  n';
text(1100,160,txt)
return
clc; clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [filename, filepath] = uigetfile('*.wav','Select wave file');
% [upperPath, deepestFolder, ~] = fileparts(filepath(1:length(filepath)-1));
filepath='E:\PhD\datasets\hillenbrand\allx\';
% filepath='E:\PhD\datasets\hillenbrand\allxbp400_3400\';
currentdir = filepath;
dir1 = currentdir;
fnames = dir(fullfile([dir1 '*.wav']));
% return
size1=size(fnames,1);
data_x={};
data_s={};
power_x={};
% return
% len = 1024;
alpha1=0.97;
startpos=1;
% 
% for index=1:size1
% %   fullfilename = [dir1 namesall{index}]; % for timit
%     fullfilename = [dir1 fnames(index).name];  %for nontimit
%     [x fs] = audioread(fullfilename);
% %     x=[x;x];
%     l1=size(x,1);
%     p1=nextpow2(l1);
%     l2=2.^p1;
%     if l2<4096
%         l2=4096;
%     end
%     len=l2;    
%     x=padarray(x,l2-l1,'post');
% 
% %     s = filter( [1   -alpha1], 1, x ); % pre-emphasize    
%     x=bandpass(x,[400 3400],fs);
%     x=bandpass(x,[400 3400],fs);        
% %     s=bandpass(s,[400 3400],fs);
% %     s=bandpass(s,[400 3400],fs);        
%     data_x{index}=x;
% %     data_s{index}=s;
% 
% len=-1;
% len=1024;
% 
% % x=x(1:len);
% % x = x .* hamm;
% % audiowrite([dir1 '00_' fnames(index).name],x,fs);
% 
% %     s=data_s{index};
%     if len==-1
%     len=length(x);
%     else
%     x = x(1:len);
% %     s = s(1:len);
%     end
% %     return
%     res = fs/len;
% %     t=(0:len-1)/fs;       % times of sampling instants
%     [power1 freqs] = runfft(x,fs,startpos,len);
%     power1 = power1*1000;
%     power_x{index}=power1;
% 
% if rem(index,100)==0
%     index
% end    
% 
% end
% return

load power_hil.mat
load data_x_hil.mat

% load power_hil_bp.mat
% load data_x_hil_bp.mat
 
load fnames_hil.mat

params=[];
e=[];
l=[];
e_min=1000000000;
ll_min=1000000000;
best=[];
% for ii9=1:1
% for ii8=1:1% % calculates lpc, cepstrum, yin, yaapt, and new formants for ii7=1:1
% for ii7=1:1
% for ii6=1:1
% for ii5=1:1
% for ii4=1:10
% for ii3=1:10
% for ii2=1:10
% for ii1=1:13
% for ii=1:1

% jj1 and jj2 may not be needed
for jj1=1:1
    rx1=jj1;
    for jj2=1:1
        rx2=jj2;

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

% rx1=10;
% rx2=5;
% max_formant_count=13;
% minf0=65;
% hi_freq=1300; % start of high frequency region
% ratio1=40; % minimum ratio for frequencies below hi_freq
% ratio2=25; % minimum ratio for frequencies above hi_freq
% low_freq1=1900; % this is not used
% low_freq2=500; % this is not used
% bank_bwidth=150;        
        
% rx1=10;
% rx2=5;
% max_formant_count=13;
% minf0=65;
% hi_freq=1300; % start of high frequency region
% ratio1=40; % minimum ratio for frequencies below hi_freq
% ratio2=25; % minimum ratio for frequencies above hi_freq
% low_freq1=1900; % this is not used
% low_freq2=500; % this is not used
% bank_bwidth=150;
    
% % len = power(2,ii9+7);
% % max_formant_count=4+ii1;
% % minf0=45+ii2*5;
% % hi_freq=900+ii3*100; % start of high frequency region
% % ratio1=5+ii4*5; % minimum ratio for frequencies below hi_freq
% % ratio2=5+ii5*5; % minimum ratio for frequencies above hi_freq
% % we delete formants which are less than 1/10th of the highest
% % formant below hi_freq
% % and formants which are 1/20th of the highest formant
% % above hi_freq
% % bank_bwidth=150;
% % equ_loudness=1; % if 1 apply equal loudness
% % low_freq1=1900; % this is not used
% % low_freq2=500; % this is not used

% alpha1=0.97;
% allfeatures=[];
allf0=[];
h1=0;
h2=0;
h3=0;
h4=0;
h5=0;
t1=tic;
tt=[];
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
% %     % get lpc formants
% t1=tic;
%     [power1 freqs] = runfft(x,fs,startpos,len);
%     power1 = power1*1000;
%     [formants f h]=get_lpc_formants(s,fs,len);
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
[f0_1 amp0_1 f0_1x amp0_1x]=get_f0_2(power1,res,len,minf0,hi_freq,ratio1,ratio2,max_formant_count,rx1,rx2);
t2=toc(t1);
% tt=[tt;t2];
% h1=h1+t2;
% return
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YIN
t1c=tic;
f0yin=YIN(xx,fs);
t2c=toc(t1c);
tt=[tt;t2c];
h4=h4+t2c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% yaapt
% [f0yaapt1, nf1] = yaapt(x_yaapt, fs, 0, [], 0, 1);
% [f0yaapt2, nf2] = yaapt(x_yaapt, fs, 0, [], 0, 2);
% ham2 = hamming(2048); 
xxx=padarray(x,1024,'post');
x_yaapt = xxx;%.*ham2; 
t1d=tic;
[f0yaapt3, nf3] = yaapt(x_yaapt, fs, 0, [], 0, 3);
t2d=toc(t1d);
tt=[tt;t2d];
h5=h5+t2d;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     allf0=[allf0; f0yaapt1(1,1) f0yaapt2(1,1) f0yaapt3(1,1)];
%     allf0=[allf0; f0_1 f0_1x f0_2 f0_2x f0_acorr f0_ceps f0yin(1,1) f0yaapt3(1,1)];
%     allf0=[allf0;f0_1 f0_1x f0_2 f0_2x f0_3 f0_3x f0_4 f0_4x f0_acorr];
%     allf0=[allf0;f0_acorr f0_ceps];
%     allf0=[allf0;f0_1 f0_1x];
    allf0=[allf0; f0_1 f0_1x f0_acorr f0_ceps f0yin(1,1) f0yaapt3(1,1) ];
%     allf0=[allf0; f0_1 f0_1x f0_acorr f0_ceps];
%     allf0=[allf0; f0_1 f0_1x];

end

% xx=sum(tt)/1668

% return
% use hillenb_gt2 for correct ground truth of fnames for files in allx folder
% replace the rest of the method, it is necessary
% fnames is sorted to labels
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

% return

end
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



return
% end
% end
% end
% end
% 
% 
% end
% end
% end
% end
% end
% end
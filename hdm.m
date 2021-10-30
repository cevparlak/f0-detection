function [f0 amp0 f0x amp0x]=get_f0_1(powerx,res,len,minf0,hi_freq,ratio1,ratio2,max_formant_count,rx1,rx2)
%function [Flist f0 amp0 f0x f0x2 rx1 rx2 list1 listequ stepsequ]=load_formants(powerx,power,res,len,minf0,low_freq1,low_freq2,hi_freq,ratio1,ratio2,max_formant_count,equ_loudness,bank_bwidth,sort)
% get listx from fft power
listequ=[];
f0listx = [];
max1=max(powerx);
maxfreq=hi_freq+minf0;
maxpow=max1/ratio1;
for i=3:len/10 %hi_freq/res
    f1=i*res;
    p1=powerx(i);
    p2=powerx(i-1);
    p3=powerx(i+1);
    
    if f1>=minf0
        if f1<=maxfreq
            if p1>maxpow
                if(p1>p2 && p1>p3)
                    f0listx=[f0listx; f1 p1];
                end
                
            end
        end
    end
end
% xx=f0listx;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove unnecessary peaks which are nearer than minf0 using frequency values
size1=size(f0listx,1);
if size1>0
    if size1>2
        for i=2:size1-1
          if   f0listx(i,1)>0
            if f0listx(i,1)+minf0>f0listx(i+1,1)
                if f0listx(i,2)<f0listx(i+1,2)
                    f0listx(i,1)=0;
                end
            end
            if f0listx(i,1)-minf0<f0listx(i-1,1)
                if f0listx(i,2)<f0listx(i-1,2)
                    f0listx(i,1)=0;
                end
            end
          end
            
        end
    end
    f0listx = f0listx(f0listx(:,1)~=0,:);
    size1=size(f0listx,1);
    % check first entry
    if size1>1
      if f0listx(1,1)>0        
        if f0listx(1,1)+minf0>f0listx(2,1)
            if f0listx(1,2)<f0listx(2,2)
                f0listx(1,1)=0;
            end
        end
      end
    end
    f0listx = f0listx(f0listx(:,1)~=0,:);
    size1=size(f0listx,1);
    % check last entry
    if size1>1
      if f0listx(size1,1)        
        if f0listx(size1,1)-minf0<f0listx(size1-1,1)
            if f0listx(size1,2)<f0listx(size1-1,2)
                f0listx(size1,1)=0;
            end
        end
      end
    end
    
    f0listx = f0listx(f0listx(:,1)~=0,:);
% xxx=f0listx;    
% continue to remove redundant peaks using amplitudes
    size1=size(f0listx,1);
    if size1>2
        for i=2:size1-1
          if f0listx(i,2)>0            
            if f0listx(i,2)<f0listx(i+1,2)/rx1
                f0listx(i,1)=0;
            end
            if f0listx(i,2)<f0listx(i-1,2)/rx1
                f0listx(i,1)=0;
            end
          end
            
        end
    end
    f0listx = f0listx(f0listx(:,1)~=0,:);
size1=size(f0listx,1);
    % check first entry
    if size1>1
      if f0listx(1,2)>0            
        
        if f0listx(1,2)<f0listx(2,2)/rx1
            f0listx(1,1)=0;
        end
      end
    end
    f0listx = f0listx(f0listx(:,1)~=0,:);
size1=size(f0listx,1);
    % check last entry
    if size1>1
     if f0listx(size1,2)>0            
        
        if f0listx(size1,2)<f0listx(size1-1,2)/rx1
            f0listx(size1,1)=0;
        end
     end
    end
    
    % handle some special cases
    f0listx = f0listx(f0listx(:,1)~=0,:);
    
    size1=size(f0listx,1);
    if f0listx(1,1)<110
        
        if size1>1
            if f0listx(1,2)<f0listx(2,2)/rx2
                f0listx(1,1)=0;
            end
        end
    end
    
    f0listx = f0listx(f0listx(:,1)~=0,:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%
    % extract f0 using differences between the frequencies of the peaks
    list1=[];
    size1=size(f0listx,1);
    if size1>1
        
        for i=1:size1-1
            diff1=abs(f0listx(i+1,1)-f0listx(i,1));
            if diff1>minf0
             list1=[list1;f0listx(i,1) f0listx(i+1,1) diff1];
            end
            if diff1<=minf0
                if f0listx(i,2)>f0listx(i+1,2)
                 list1=[list1;f0listx(i,1) f0listx(i+1,1) diff1];
                else
                 list1=[list1;f0listx(i+1,1) f0listx(i+1,1) diff1];
                end
            end
            
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find the unique values of differences
    
% if size(list1,1)==0    
%     list1
%     f0listx
%     xx
%     xxx
% end

    vals = unique(list1(:,3));

    
    % Count the number of instances of each of the unique values
    s1=size(vals,1);
    if s1>1
        valCount = hist(list1(:,3),vals)';
    else
        valCount=1;
    end
    % f0 is the most repeated difference
    list2 = [vals valCount];
%     list2 = list2(list2(:,1)>minf0,:);  % may be unnecessary
    if size(list2,1)==0
        list2=[list1(1,1) 1];
    end
    list2=sortrows(list2,-2);
    f0x=list2(1,1);
    
    diff1=f0listx(:,1)-f0x;
    list3=[f0listx(:,1) f0listx(:,2) diff1];
    
    list3=sortrows(list3,3);
    
    list3 = list3(list3(:,3)>=0,:);
    f0=list3(1,1);
    amp0=list3(1,2);
    
    if abs(f0-f0x)<minf0/2
        f0=f0x;
    end    
    amp0x=amp0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        % we have only 1 proper f0 candidate
%         diff1=f0listx(1,1);
%         list1=[f0listx(1,1) f0listx(1,2) diff1];
     f0=f0listx(1,1);
     amp0=f0listx(1,2);
     f0x=f0;
     amp0x=amp0;
    end
  
else
    f0=0;
    amp0=0;
    f0x=0;
    amp0x=0;
end
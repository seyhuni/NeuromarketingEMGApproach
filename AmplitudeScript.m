clear; clc; close all;
addpath(genpath('C:\Users\asus\Desktop\linacalisma\'))
[SwallowFileName, path] = uigetfile('*.*','.mat');
cd(path)
File=fullfile(path,SwallowFileName);
EMGData=load(File);
EMGData=EMGData.datas;
EMGData_=EMGData(500:end);
EMGData_=EMGData_-mean(EMGData_);
EMGData_=EMGData_/max(EMGData_);
EMGData_=EMGData_(EMGData_>-1);
% EMGData_(163000:165000)=0;
figure
plot(EMGData_)
xlabel('Time')
ylabel('Amplitude')
xlim([0 length(EMGData_)])
% ylim([-1 1.1])
% [c,l]=wavedec(datas3,6,'sym4');
% a8 = wrcoef('a',c,l,'sym4',6);
% block_duration_datacorrected=datas3-a8;
fs=2000;
% faxis=linspace(-fs/2,fs/2,length(datas3));
% plot(faxis,fftshift(abs(fft(datas3))));
block_duration_datacorrected=EMGData_;
Band    = (2 / fs) * [25, 400];
[B, A]  = butter(6, Band, 'Bandpass');   
fSignal = filtfilt(B, A, double(block_duration_datacorrected));
    % figure
 Band    = (2 / fs) * [45, 55];
[B, A]  = butter(4, Band, 'stop');   
block_duration_datacorrected_last = filtfilt(B, A, double(block_duration_datacorrected));
    N=length(block_duration_datacorrected_last);         %number of points
    t=(0:N-1)/fs;   %time vector
    sgf = sgolayfilt(block_duration_datacorrected_last,3,201);
    % figure
    % plot(sgf,'r')
    % title(['Savitzky-Golay Smoothed - block: ' num2str(blocks{block_num})])
    sgf=sgf(101:end-101);
    ynew=sgf/max(sgf); 
    % initialize filtered signal
    eogF = ynew;
%     figure
%     plot(eogF)
%     xlabel('Time')
%     ylabel('Amplitude')
%     xlim([1020.231193413892 260567.0467979074])
%     ylim([-0.0001904646550288437 0.0002626080637102741])
    % TKEO basic  % Teager–Kaiser energy operator to obtain EMG Bursts
    for i=2:length(eogF)-1
        eogF(i) = ynew(i)^2 - ynew(i-1)*ynew(i+1);
    end   
    % eogF=eogF(1:end-1);
    [c,l] = wavedec(eogF,7,'sym4'); % 8 level decomposition  

    % for t=1:8
    %     D(:,t)=wrcoef('d',c,l,'db6',t);
    % end

    % low frequency components to get swallow patterns, so that approximation
    % of wavelets are obtained
     clear A;
    for t=1:7
        A(:,t)=wrcoef('a',c,l,'sym4',t);
    end

    A8=A(:,7);  % A8 is the filtered and swallow pattern obtained signal

    % figure
    % plot(A8)
    A8_=A8(500:end-500);
%     figure
%     plot(A8_)
%     xlabel('Time')
%     ylabel('Amplitude')
            rmsSwallows = sqrt(movmean(A8_.^2, 2500));   % Burst Detection using RMS Value Over ‘WinLen’ Samples
    figure(3)
    plot(rmsSwallows)
    button = 1;

    while sum(button) <=1   % read input with right click only for threshold
       [xx,yy,button] = ginput(1); % you can choose threshold by right clicking on plot 1 time
    end  % you can change
     thresholdValue=(yy);
     close(figure(3)) 
     x=2000
        for k=x:(length(rmsSwallows)-x)
            interval_=rmsSwallows(k-x+1:k+x);
            if(interval_(ceil(length(interval_)/2))==max(interval_) & interval_(ceil(length(interval_)/2))>thresholdValue)
                swallowpeaks(k) = rmsSwallows(k);
            end
        end
    swallowpeaks_pos=find(swallowpeaks>thresholdValue);
    hv=ones(1,length(rmsSwallows));
    swallowoffset=sqrt(-1)*hv;
    swallowonset=swallowoffset;
    
    for k=1:length(swallowpeaks_pos)
        mt=rmsSwallows(swallowpeaks_pos(k):swallowpeaks_pos(k)+4500);
        mn=min(mt);
        swallowoffset(find(mt==mn)+swallowpeaks_pos(k)-1)=mn;
    end
    swallowoffset_last=find(swallowoffset~=sqrt(-1));

    for m=2:length(swallowpeaks_pos)
        mt2=rmsSwallows(swallowpeaks_pos(m)-5000:swallowpeaks_pos(m));
        mn2=min(mt2);
        swallowonset(swallowpeaks_pos(m)-5000-1+find(mt2==mn2))= mn2;
    end
    swallowonset_last=find(swallowonset~=sqrt(-1));
    
    rmsSwallows=rmsSwallows/max(rmsSwallows);
    figure
    plot(rmsSwallows)
    hold on
    plot(swallowpeaks_pos,rmsSwallows(swallowpeaks_pos),'r+','MarkerFaceColor','r','LineWidth',1)
    hold on
    plot(swallowonset_last,rmsSwallows(swallowonset_last),'g*','MarkerFaceColor','g','LineWidth',1)
    hold on
    plot(swallowoffset_last,rmsSwallows(swallowoffset_last),'k+','MarkerFaceColor','k','LineWidth',1)
    xlabel('Time')
    ylabel('Amplitude')
    
swallow_amp=zeros(1,length(swallowpeaks_pos));
% format longG
for i=1:length(swallow_amp)
    swallow_amp(i)=rmsSwallows(swallowpeaks_pos(i))-rmsSwallows(swallowoffset_last(i));
end
format short g

BSAmp=0;
ASAmp=0;
BlackScreenSwallows=find(swallowpeaks_pos<60*fs);
AdvSwallows=find(swallowpeaks_pos>=60*fs);
for i=1:length(BlackScreenSwallows)
    BSAmp=BSAmp+swallow_amp(i);
end
for j=(length(BlackScreenSwallows)+1):length(BlackScreenSwallows)+length(AdvSwallows)
    ASAmp=ASAmp+swallow_amp(j);
end
MeanBlackScreen=BSAmp/length(BlackScreenSwallows)
MeanAdv=ASAmp/length(AdvSwallows)
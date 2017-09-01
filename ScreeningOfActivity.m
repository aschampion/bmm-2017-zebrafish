cd ~/'Dropbox (MIT)'/MyData/[2017]BMMcourse/Jing;
load('6dpf_2026x2061_fish00_osc01_20170731_232605.mat') ;
%% take a subset of data 
CellInd=randsample(size(Cell_baseline1,1),1000)';
Cell_baseline1_Sample=Cell_baseline1(CellInd,:);
Cell_timesers1_Sample=Cell_timesers1(CellInd,:);
Cell_X_Sample=Cell_X(CellInd,:);
Cell_Y_Sample=Cell_Y(CellInd,:);
Cell_Z_Sample=Cell_Z(CellInd,:);
% cut the begining 
Cell_baseline1_Sample=Cell_baseline1_Sample(:,2e3:end);
Cell_timesers1_Sample=Cell_timesers1_Sample(:,2e3:end);
CellActSample=Cell_timesers1_Sample-Cell_baseline1_Sample;
%% 
figure
CellActSampleZS=zscore(CellActSample');
colormap(flipud(colormap(('bone'))));
s=imagesc(CellActSampleZS);colorbar;
mesh(1:100,1:length(CellActSample),CellActSample)
set(gca,'xtick',floor(linspace(1,length(timeInd),10)))
set(gca,'xticklabel',num2str(floor(timeInd(floor(linspace(1,length(timeInd),10))))'))

%%  
close all;
f=figure;
v = VideoWriter('activity1','MPEG-4');
open(v);
A=max(Cell_timesers1_Sample'-Cell_baseline1_Sample');
Xlim=[min(Cell_X(:,1)),max(Cell_X(:,1))];
Ylim=[min(Cell_Y(:,1)),max(Cell_Y(:,1))];
% 
steps=100;

Ktimes=floor(linspace(1,size(Cell_timesers1,2),steps));
P=subplot('position',[.1,.1,.8,.8]);
for i=Ktimes
    CellAct=-Cell_baseline1_Sample(:,i)+Cell_timesers1_Sample(:,i);
    %CellAct(CellAct==nan)=0;
    CellColor=(CellAct./A')*[1,0,0];
    s=scatter3(Cell_X_Sample(:,1),Cell_Y_Sample(:,1),Cell_Z_Sample(:,1),20,'filled','CData', CellColor);
    daspect([1,1,1]);set(gca,'XLim',Xlim,'YLim',Ylim);
    text(max(get(gca,'xlim')),max(get(gca,'ylim')),num2str(i));
    %drawnow
    F = getframe(f);
    writeVideo(v,F.cdata)
end;
close(v); 

%% run FFT on neurons 
tSample= 1/83.5 %(seconds);
freqMat=[];EnergyMat=[];
for i=1:size(CellActSample,1)
    temp=CellActSample(i,:);
    [freq,mag,power]=computeFFT(temp,tSample,'NoShowPlot');
    freqMat=[freqMat;freq];
    EnergyMat=[EnergyMat;(power)];    
end
% take the mean 
figure;
plot(mean(freqMat,1),mean(10*log10(EnergyMat),1),'-k')
xlabel('Freq');ylabel('Power(db)');grid on;
%% calculate autocorrelation: 
tSample= 1/83.5 %(seconds);
timeLagMat=[];autoCorrMat=[];
numLags=500;
for i=1:size(CellActSample,1)
    temp=CellActSample(i,:);
    [acf,lags]=autocorr(temp,numLags);
    timeLagMat=[timeLagMat;lags*tSample];
    autoCorrMat=[autoCorrMat;acf];    
end
% 
% figure
% for i=1:size(CellActSample,1)
%     plot3(timeLagMat(1,:),i*ones(size(timeLagMat,2)),autoCorrMat(i,:),'ok')
%     hold on 
% end
figure;
plot(mean(timeLagMat,1),mean((autoCorrMat),1),'.k');
grid on;
xlabel('Lag (secs)');ylabel('Auto correlation');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 20170823
clear all;close all;home
cd('~/Dropbox (MIT)/MyData/[2017]BMMcourse/yumu_20160531_singleplane_83.5Hz_mika_processing');
cd plane9/
load('Cells0_data.mat');
forberainX=(Cell_X(:,1)<200 & Cell_X(:,1)>100);
forberainY=(Cell_Y(:,1)>175 & Cell_Y(:,1)<375);
forbrainInd=find(forberainX.*forberainY);
tactumX=(Cell_X(:,1)<315 & Cell_X(:,1)>215);
tactumY=(Cell_Y(:,1)<485 & Cell_Y(:,1)>350);
tactumInd=find(tactumX.*tactumY);
figure; 
scatter(Cell_X(:,1),Cell_Y(:,1));
hold on 
scatter(Cell_X(find(forberainX.*forberainY),1),Cell_Y(find(forberainX.*forberainY),1),10,'filled','r');
scatter(Cell_X(find(tactumX.*tactumY),1),Cell_Y(find(tactumX.*tactumY),1),10,'filled','c');

%% 
tSample= 1/83.5; %(seconds);
cellAct=Cell_timesers1-Cell_baseline1;
forbrainAct=cellAct(forbrainInd,6e3:end);
tectumAct=cellAct(tactumInd,6e3:end);
allAct=[forbrainAct;tectumAct];
[ActPercentile]=prctile(allAct',[5,95]);
minAct=min(ActPercentile(1,:));
maxAct=max(ActPercentile(2,:));
timeInd=[0:size(tectumAct,2)]*tSample;
figure
subplot(2,1,1)
m=256;
cm_viridis=viridis(m);
%colormap(flipud(colormap(('bone'))));
s=imagesc(forbrainAct);
h=colorbar;ylabel(h,'\Deltaf/f');colormap(cm_viridis);
set(gca,'xtick',floor(linspace(1,length(timeInd),10)))
set(gca,'xticklabel',num2str(floor(timeInd(floor(linspace(1,length(timeInd),10))))'))

ylabel('Neurons ')
title('forebrain')


subplot(2,1,2)
%colormap(flipud(colormap(('bone'))));
s=imagesc(tectumAct);h=colorbar;ylabel(h,'\Deltaf/f');colormap(cm_viridis);
set(gca,'xtick',floor(linspace(1,length(timeInd),10)))
set(gca,'xticklabel',num2str(floor(timeInd(floor(linspace(1,length(timeInd),10))))'))

title('Tactum');
ylabel('Neuron');xlabel('time (sec)');
%% fft on cells in each area 
tSample= 1/83.5; %(seconds);
freqMat=[];EnergyMat=[];
for i=1:size(forbrainAct,1)
    temp=forbrainAct(i,:);
    [freq,mag,power]=computeFFT(temp,tSample,'NoShowPlot');
    freqMat=[freqMat;freq];
    EnergyMat=[EnergyMat;(power)];    
end

figure;
subplot(2,1,1)
plot(mean(freqMat,1),mean(10*log10(EnergyMat),1),'-k')
xlabel('Freq');ylabel('Power(db)');grid on;title('Forebrain')
% 
freqMat=[];EnergyMat=[];
for i=1:size(tectumAct,1)
    temp=tectumAct(i,:);
    [freq,mag,power]=computeFFT(temp,tSample,'NoShowPlot');
    freqMat=[freqMat;freq];
    EnergyMat=[EnergyMat;(power)];    
end
% take the mean 
subplot(2,1,2)
plot(mean(freqMat,1),mean(10*log10(EnergyMat),1),'-k')
xlabel('Freq');ylabel('Power(db)');grid on;title('Tactum')
%% autocorrelation on cells in each area 
timeLagMat=[];autoCorrMat=[];
numLags=500;
for i=1:size(forbrainAct,1)
    temp=forbrainAct(i,:)-mean(forbrainAct(i,:));
    [acf,lags]=autocorr(temp,numLags);
    timeLagMat=[timeLagMat;lags*tSample];
    autoCorrMat=[autoCorrMat;acf];    
end
figure;
subplot(2,1,1)
plot(mean(timeLagMat,1),mean((autoCorrMat),1),'.k');
grid on;
xlabel('Lag (secs)');ylabel('Auto correlation');title('Forebrain')

% 
timeLagMat=[];autoCorrMat=[];
numLags=500;
for i=1:size(tectumAct,1)
    temp=tectumAct(i,:);
    [acf,lags]=autocorr(temp,numLags);
    timeLagMat=[timeLagMat;lags*tSample];
    autoCorrMat=[autoCorrMat;acf];    
end
hold on
subplot(2,1,2)
plot(mean(timeLagMat,1),mean((autoCorrMat),1),'.k');
grid on;
xlabel('Lag (secs)');ylabel('Auto correlation');title('Tactum')

%%  look at gauss 
tSample= 1/83.5; %(seconds);
cellAct=Cell_timesers1-Cell_baseline1;
forbrainAct=cellAct(forbrainInd,6e3:end);
tectumAct=cellAct(tactumInd,6e3:end);
timeInd=[0:size(tectumAct,2)]*tSample;
cd(['/Users/eghbalhosseiniasl1/Dropbox (MIT)/MyData/[2017]BMMcourse/yumu_20160531_singleplane_83.5Hz_mika_processing/plane9/GaussProcess'])
load('ForebrainGaus_v3.mat')
figure
subplot(2,1,1)
colormap(flipud(colormap(('bone'))));
s=imagesc(forbrainAct-mean(forbrainAct,2)*ones(1,size(forbrainAct,2)));colorbar
set(gca,'xtick',floor(linspace(1,length(timeInd),10)))
set(gca,'xticklabel',num2str(floor(timeInd(floor(linspace(1,length(timeInd),10))))'))

ylabel('Neurons ')
title('forebrain')

subplot(2,1,2)
colormap(flipud(colormap(('bone'))));
s=imagesc(ForebrainGausEstim);colorbar
set(gca,'xtick',floor(linspace(1,length(timeInd),10)))
set(gca,'xticklabel',num2str(floor(timeInd(floor(linspace(1,length(timeInd),10))))'))
%% 
tSample= 1/83.5; %(seconds);
freqMat=[];EnergyMat=[];
forebrainAct=forbrainAct-mean(forbrainAct,2)*ones(1,size(forbrainAct,2));
forebrainActTemp=forbrainAct;
forbrainAct=forbrainAct(:,1:2000);
for i=1:size(forbrainAct,1)
    temp=forbrainAct(i,:);
    [freq,mag,power]=computeFFT(temp,tSample,'NoShowPlot');
    freqMat=[freqMat;freq];
    EnergyMat=[EnergyMat;(power)];    
end

figure;
subplot(2,1,1)
plot(mean(freqMat,1),mean(10*log10(EnergyMat),1),'-k')
xlabel('Freq');ylabel('Power(db)');grid on;title('Forebrain')
ylim([-100,0]);xlim([-40,40])
% 
ForebrainGausEstim=ForebrainGausProcess.GaussEstimate;
freqMat=[];EnergyMat=[];
ForebrainGausEstim=ForebrainGausEstim-mean(ForebrainGausEstim,2)*ones(1,size(ForebrainGausEstim,2));
for i=1:size(ForebrainGausEstim,1)
    temp=ForebrainGausEstim(i,:);
    [freq,mag,power]=computeFFT(temp,tSample,'NoShowPlot');
    freqMat=[freqMat;freq];
    EnergyMat=[EnergyMat;(power)];    
end


subplot(2,1,2)
plot(mean(freqMat,1),mean(10*log10(EnergyMat),1),'-k')
xlabel('Freq');ylabel('Power(db)');grid on;title('Forebrain')
ylim([-100,0]);xlim([-40,40])

%% 2017/08/31
clear all;close all;home
cd('~/Dropbox (MIT)/MyData/[2017]BMMcourse/yumu_20160531_singleplane_83.5Hz_mika_processing');
cd plane9/
load('Cells0_data.mat');
forberainX=(Cell_X(:,1)<200 & Cell_X(:,1)>100);
forberainY=(Cell_Y(:,1)>175 & Cell_Y(:,1)<375);
forbrainInd=find(forberainX.*forberainY);
tactumX=(Cell_X(:,1)<315 & Cell_X(:,1)>215);
tactumY=(Cell_Y(:,1)<485 & Cell_Y(:,1)>350);
tactumInd=find(tactumX.*tactumY);
%
%% 
tSample= 1/83.5; %(seconds);
cellAct=Cell_timesers1-Cell_baseline1;
forbrainAct=cellAct(forbrainInd,6e3:end);
tectumAct=cellAct(tactumInd,6e3:end);
timeInd=[0:size(tectumAct,2)]*tSample;

close all;
f=figure;
set(f,'color',[1,1,1])
v = VideoWriter('activity1','MPEG-4');
open(v);
A=max(forbrainAct');
B=max(tectumAct');
Xlim=[min(Cell_X(:,1)),max(Cell_X(:,1))];
Ylim=[min(Cell_Y(:,1)),max(Cell_Y(:,1))];
% 
SpeedupRatio=1*1/tSample;
Ktimes=floor([1:SpeedupRatio:size(forbrainAct,2)]);
steps=150;
t=[0:(size(forbrainAct,2)-1)]*tSample;
%Ktimes=floor(linspace(1,size(forbrainAct,2),steps));
P=subplot('position',[.1,.1,.8,.8]);


for i=Ktimes
    CellAct=forbrainAct(:,i);
    CellColor=((CellAct-min(CellAct))./A')*[1,0,0];
    tectum=tectumAct(:,i);
    tactumColor=((tectum-min(tectum))./B')*[1,0,0];
    scatter(Cell_X(:,1),Cell_Y(:,1),'cdata',[.5,.5,.5]);
    hold on 
    s=scatter(Cell_X(forbrainInd,1),Cell_Y(forbrainInd,1),20,'filled','CData', CellColor);
    s=scatter(Cell_X(tactumInd,1),Cell_Y(tactumInd,1),20,'filled','CData', tactumColor);
    hold off
    daspect([1,1,1]);set(gca,'XLim',Xlim,'YLim',Ylim); axis off
    text(max(get(gca,'xlim')),max(get(gca,'ylim')),sprintf('%0.1f',t(i)));
    %drawnow
    F = getframe(f);
    writeVideo(v,F.cdata)
end;
close(v); 

%% 
%% 
tSample= 1/83.5; %(seconds);
cellAct=Cell_timesers1-Cell_baseline1;
fullBrainAct=cellAct(:,6e3:end);
timeInd=[0:size(tectumAct,2)]*tSample;

close all;
f=figure;
set(f,'color',[1,1,1])
v = VideoWriter('activityBrain','MPEG-4');
open(v);
fullBrainActGray=fullBrainAct-[min(fullBrainAct')]'*ones(1,size(fullBrainAct,2));
fullBrainActGray=fullBrainActGray.*([max(fullBrainAct')-min(fullBrainAct')]'*ones(1,size(fullBrainAct,2))).^-1;
Xlim=[min(Cell_X(:,1)),max(Cell_X(:,1))];
Ylim=[min(Cell_Y(:,1)),max(Cell_Y(:,1))];
% 
SpeedupRatio=1*1/tSample;
Ktimes=floor([1:SpeedupRatio:size(fullBrainAct,2)]);
steps=150;
t=[0:(size(fullBrainAct,2)-1)]*tSample;
%Ktimes=floor(linspace(1,size(fullBrainAct,2),steps));
P=subplot('position',[.1,.1,.8,.8]);


for i=Ktimes
    CellAct=fullBrainActGray(:,i);
    CellColor=(CellAct)*[1,0,0];
    scatter(Cell_X(:,1),Cell_Y(:,1),'cdata',[.5,.5,.5]);
    hold on 
    s=scatter(Cell_X(:,1),Cell_Y(:,1),20,'filled','CData', CellColor);
    hold off
    daspect([1,1,1]);set(gca,'XLim',Xlim,'YLim',Ylim); axis off
    text(max(get(gca,'xlim')),max(get(gca,'ylim')),sprintf('%0.1f',t(i)));
    %drawnow
    F = getframe(f);
    writeVideo(v,F.cdata)
end;
close(v); 

%% 
%% 
tSample= 1/83.5; %(seconds);
cellAct=Cell_timesers1-Cell_baseline1;
forbrainAct=cellAct(forbrainInd,6e3:end);
tectumAct=cellAct(tactumInd,6e3:end);
timeInd=[0:size(tectumAct,2)]*tSample;

close all;
f=figure;
set(f,'color',[1,1,1])
v = VideoWriter('activity3','MPEG-4');
open(v);
A=max(forbrainAct');
B=max(tectumAct');
Xlim=[min(Cell_X(:,1)),max(Cell_X(:,1))];
Ylim=[min(Cell_Y(:,1)),max(Cell_Y(:,1))];
% 
SpeedupRatio=1*1/tSample;
Ktimes=floor([1:SpeedupRatio:size(forbrainAct,2)]);
steps=150;
t=[0:(size(forbrainAct,2)-1)]*tSample;
%Ktimes=floor(linspace(1,size(forbrainAct,2),steps));
P=subplot('position',[.1,.1,.8,.8]);
forbrainActGray=mat2gray(forbrainAct);
forbrainActGray=forbrainActGray.*([max(forbrainActGray')]'*ones(1,size(forbrainActGray,2))).^-1;
tectumActGray=mat2gray(tectumAct);
tectumActGray=(tectumActGray.*([max(tectumActGray')]'*ones(1,size(tectumActGray,2))).^-1);
for i=Ktimes
    CellAct=forbrainActGray(:,i);
    CellColor=(CellAct)*[1,0,0];
    tectum=tectumActGray(:,i);
    tactumColor=(tectum)*[1,0,0];
    scatter(Cell_X(:,1),Cell_Y(:,1),'cdata',[.5,.5,.5]);
    hold on 
    s=scatter(Cell_X(forbrainInd,1),Cell_Y(forbrainInd,1),20,'filled','CData', CellColor);
    s=scatter(Cell_X(tactumInd,1),Cell_Y(tactumInd,1),20,'filled','CData', tactumColor);
    hold off
    daspect([1,1,1]);set(gca,'XLim',Xlim,'YLim',Ylim); axis off
    text(max(get(gca,'xlim')),max(get(gca,'ylim')),sprintf(' %0.1f secs',t(i)));
    %drawnow
    F = getframe(f);
    writeVideo(v,F.cdata)
end;
close(v); 

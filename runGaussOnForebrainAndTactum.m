load('Cells0_data.mat');
forberainX=(Cell_X(:,1)<200 & Cell_X(:,1)>100);
forberainY=(Cell_Y(:,1)>175 & Cell_Y(:,1)<375);
forbrainInd=find(forberainX.*forberainY);
tactumX=(Cell_X(:,1)<315 & Cell_X(:,1)>215);
tactumY=(Cell_Y(:,1)<485 & Cell_Y(:,1)>350);
tactumInd=find(tactumX.*tactumY);
% 
tSample= 1/83.5; %(seconds);
cellAct=Cell_timesers1-Cell_baseline1;
forbrainAct=cellAct(forbrainInd,6e3:end);
tectumAct=cellAct(tactumInd,6e3:end);
timeInd=[0:size(tectumAct,2)]*tSample;
% 
ForebrainGausProcess=GaussProcessFreq(forbrainAct,1/tSample);
disp('forebrain done!')
TactumGaussProcess=GaussProcessFreq(tectumAct,1/tSample);
disp('Tactum done!')
% 
save('ForebrainGaus','ForebrainGausProcess','-v7.3');
save('TactumGaus','TactumGaussProcess','-v7.3');

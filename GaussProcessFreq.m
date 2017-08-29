% this function estimate a gaussian process of the spatio temporal input up
% to a certain frequency 
% method: constructing a tensor N*N*w where N is the number of temporal
% seqeuences, and w is the set of frequencies. each N*N matrix is then used
% to create a vector z that contain the temporal equivalent of the
% sequence. 
function varargout=GaussProcessFreq(varargin)
% inputs 
TimeSeries=varargin{1};
Freq=varargin{2};

if size(TimeSeries,2)<size(TimeSeries,1),
    TimeSeries=TimeSeries';
end
TimeSeries=TimeSeries-mean(TimeSeries,2)*ones(1,size(TimeSeries,2));
%% 1- construct the tensor containing N*N 
test=CorrFreq(TimeSeries(1,:),TimeSeries(2,:));
CorrMatrixFreq=zeros(size(TimeSeries,1),size(TimeSeries,1),size(test,2));
for n=1:size(TimeSeries,1)
    for m=n:size(TimeSeries,1)
        temp=CorrFreq(TimeSeries(n,:),TimeSeries(m,:));
        CorrMatrixFreq(n,m,1:length(temp))=temp;
    end
    CorrMatrixFreq(n:end,n,:)=CorrMatrixFreq(n,n:end,:);
end

%% 2- contruct z(w) according to the following equation 
% z(w)=sqrt(C(w)).normal(w)
Zw=[];GausEstimate=[];
for i=1:size(CorrMatrixFreq,3)
NormSample=randn([size(TimeSeries,1),2])*[1,j]';
ZwSample=sqrtm(squeeze(CorrMatrixFreq(:,:,i)))*NormSample;
Zw=[Zw,ZwSample];
end
for i=1:size(Zw,1)
GausEstimate=[GausEstimate;ifft(([Zw(i,:),conj(fliplr(Zw(i,:)))]),'symmetric')];
end
Gauss.GaussEstimate=GausEstimate;
Gauss.CorrMatrixFreq=CorrMatrixFreq;
Gauss.Zw=Zw;
varargout{1}=Gauss;
end
function out=CorrFreq(x,y)
xpad=([x,zeros(1,length(x)+length(y)-1-length(x))]);
ypad=([fliplr(y),zeros(1,length(x)+length(y)-1-length(y))]);
out2=(fft(xpad,length(x)).*fft(ypad,length(y)));
out=out2(1:(length(out2)-1)/2+1);
end 
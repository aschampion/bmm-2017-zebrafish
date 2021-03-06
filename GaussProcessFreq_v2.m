% this function estimate a gaussian process of the spatio temporal input up
% to a certain frequency 
% method: constructing a tensor N*N*w where N is the number of temporal
% seqeuences, and w is the set of frequencies. each N*N matrix is then used
% to create a vector z that contain the temporal equivalent of the
% sequence. 
function varargout=GaussProcessFreq_v2(varargin)
% inputs 
TimeSeries=varargin{1};
Freq=varargin{2};
if size(TimeSeries,2)<size(TimeSeries,1),
    TimeSeries=TimeSeries';
end
TimeSeries=TimeSeries-mean(TimeSeries,2)*ones(1,size(TimeSeries,2));
%% 1- construct the tensor containing N*N 
test=CorrFreq(TimeSeries(1,:),TimeSeries(1,:));
CorrMatrixFreq=zeros(size(TimeSeries,1),size(TimeSeries,1),size(test,2));
for n=1:size(TimeSeries,1)
    for m=n:size(TimeSeries,1)
        [temp,pad,N]=CorrFreq(TimeSeries(n,:),TimeSeries(m,:));
        CorrMatrixFreq(n,m,1:length(temp))=temp;
    end
    CorrMatrixFreq(n:end,n,:)=CorrMatrixFreq(n,n:end,:);
end
%% 2- contruct z(w) according to the following equation 
% z(w)=sqrt(C(w)).normal(w)
Zw=[];GausEstimate=[];
smaller=(length(TimeSeries));
wMax=ceil(size(CorrMatrixFreq,3));
for i=1:wMax
NormSample=randn([size(TimeSeries,1),2])*[1,j]';
NormSampleN=NormSample./norm(NormSample);
ZwSample=sqrtm(squeeze(CorrMatrixFreq(:,:,i)))*NormSampleN;
Zw=[Zw,ZwSample];
end

for i=1:size(Zw,1)
B=([Zw(i,1:N/2),conj(fliplr(Zw(i,1:N/2)))]);
C=ifft(B,'symmetric');
D=C(smaller/2+1:pad-smaller/2+1);
GausEstimate=[GausEstimate;D];
end
Gauss.GaussEstimate=GausEstimate;
varargout{1}=Gauss;
end
function [out,pad,N]=CorrFreq(X,Y)
N=2^nextpow2(length(X)+length(Y)-1);
pad=length(X)+length(Y)-1;
XYConvfft=fft(X,N).*fft(fliplr(Y),N);
%out=XYConvfft((1:(N/4+1)));
out=XYConvfft;
end 

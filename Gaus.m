clear all 
close all 
t=0:pi/1000:pi;
tsample=pi/1000;
x=t;
y=-t;

% 
[XYcorr,lag]=xcorr(x,y,1000);
XYconv=conv(x,conj(fliplr(y)),'same');
figure
plot(t,XYconv)
hold on 
plot(t,XYcorr,'r-');shg
% RESULTS: correlation is equal to convolving x to conj of y(-t)
% 
xf=fft(x);
figure;plot(abs(fftshift(xf)));shg;
figure;plot(t,ifft(xf),t,x);shg;
% RESULTS: fft and ifft return x, fftshift is only for realigning
%
XYconvF=fft(XYconv);
XYcorrF=fft(XYcorr);
figure;plot(abs(fftshift(XYconvF)));
hold on; 
plot(abs(fftshift(XYcorrF)),'r--');shg;
% RESULTS: convolution and correlation agree. 
% 
xpad=([x,zeros(1,length(x)+length(y)-1-length(x))]);
ypad=([fliplr(y),zeros(1,length(x)+length(y)-1-length(y))]);
test=fft(xpad).*fft(ypad);
figure;plot(ifft(test));hold on;plot(XYcorr)
% RESULTS: circular convolution and FFT 
% 
xpad=([x,zeros(1,length(x)+length(y)-1-length(x))]);
ypad=([y,zeros(1,length(x)+length(y)-1-length(y))]);
test=fft(xpad).*conj(fft(ypad));
figure;plot(ifft(test));hold on;plot(XYcorr)
plot(XYcorr,'r-');shg;
% RESULT: this works!!
% 
xpad=([x,zeros(1,length(x)+length(y)-1-length(x))]);
ypad=([fliplr(x),zeros(1,length(x)+length(y)-1-length(y))]);
test=fft(xpad).*fft(ypad);
figure;plot(ifft(test));hold on;plot(xcorr(x,x,1000))

%% coding the whole process 
% 1-testing two signal in different frequencies. 
% xx , xy ,yy


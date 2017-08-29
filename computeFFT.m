function [Freq_vec,Mag,power]=computeFFT(y,tSample,ShowPlot)
N=numel(y);
Time_vec=tSample*(ceil(-N/2):ceil(N/2-1));
Fsample=1/tSample;
Freq_vec=Fsample*(ceil(-N/2):ceil(N/2-1))./N;
yshift=ifftshift(y);
Y=fft(yshift,N)/N;
Yshift=fftshift(Y);
Mag=sqrt(real(Yshift).^2+imag(Yshift).^2);
power=abs(Yshift).^2;
if strcmp('ShowPlot',ShowPlot)
figure; 
plot(Freq_vec,Mag);
xlabel('Freq(Hz)')
end
end 

%figure;
%plot(x,y);hold on;
%plot(x,yshift);
% 
%figure; 
%subplot(3,1,1);
%plot(Time_vec,y);
%xlabel('time');
%subplot(3,1,2);
%plot(Freq_vec,real(Yshift));
%subplot(3,1,3)
%plot(Freq_vec,imag(Yshift));
% 

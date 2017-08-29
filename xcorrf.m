% Matlab function for frequency domain cross correlation
function [Lag,C]=xcorrf(X,Y,L)
% X, Y ---> Input vectors 
% L --->  maximum lag (must be less than minimum of (length of X, Y)
% C ---> correlation vector
% Lag ---> lag times  
X=X(:);
Y=Y(:);
s1=size(X);
s2=size(Y);
D=min(s1(1,1),s2(1,1));
for i=1:L
    X1=ifft(fft(X(1:D-i,:)).*conj(fft(Y(i+1:D,1))));
    C(i,1)=X1(1,1);
end

C=flipud(C);
X1=ifft(fft(X(1:D,:)).*conj(fft(Y(1:D,1))));
C(L+1,1)=X1(1,1);
for i=1:L
    X1=ifft(fft(Y(1:D-i,:)).*conj(fft(X(i+1:D,1))));
    C(i+L+1,1)=X1(1,1);
end
Lag=-L:1:L;
end
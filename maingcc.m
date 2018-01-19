nclear all
close all
ad = 'C:\YiSHI\AD1974Driver\Matlab\branches\wav\test_left.wav';
[x, fs] = audioread(ad);

s1 = x(1:300,1);
s3 = x(1:300,3);
s2 = x(1:300,2);
s4 = x(1:300,4);


figure(4);
subplot(2,2,1);
plot(s1);
subplot(2,2,2);
plot(s2);
subplot(2,2,3);
plot(s3);
subplot(2,2,4);
plot(s4);
pause;
close all;


%{
Nsample = 400;

y1 c= x(1:Nsample,3);
y2 = x(1:Nsample,2);

figure(5);
subplot(211);
plot(y1);
subplot(212);
plot(y2);

FFTsize = 800; 
Y1 = fft(y1,FFTsize);
Y2 = fft(y2,FFTsize);
maxshift = Nsample;
G12 = Y1.*conj(Y2); 
R = real(ifft(G12./(abs(G12))));
%R = fftshift(R); 
[maxv,shift] = max(R);
%shift = shift - maxshift;
figure(1);
plot(R);
shift

Rxy = xcorr(y1,y2);
[maxv1,shift1] = max(Rxy);
maxshift = Nsample;
%shift1 = shift1 - maxshift;
figure(2);
plot(Rxy)
shift1
%}


Nsample = 100;
n=1:Nsample;
y1=sin(2*pi*n/50);
y2=sin(2*pi*(n+5)/50);
y1 = x(1:Nsample,3)';
y2 = x(1:Nsample,1)';

FFTLength=400;
Rxx = xcorr(y1);
Ryy = xcorr(y2);
Rxy = xcorr(y1,y2);
Sxx = fft(Rxx,FFTLength);
Syy = fft(Ryy,FFTLength);
Sxy = fft(Rxy,FFTLength);
% Unfiltered Correlation (plain correlation) 
W = ones(1,FFTLength);
% Apply the filter

G = Sxy.*W;
% Obtain the GCC
R = real(ifft(G));
G = fftshift(R,1);
[maxvalue,shift] = max(G);
figure(1);
plot(G);
shift = shift/48000
title('Plain Time Cross Correlation');


% PHAT Filter
W = 1./(1.0e-4 + abs(Sxy));
% Apply the filter
R = Sxy.*W;
% Obtain the GCC
G = fftshift(real(ifft(R)),1);
figure(2);
plot(G);
[maxvalue1,shift1] = max(G);
shift1
title('GCC PHAT ');

% ROTH Filter
W = 1./(Sxx);
R = Sxy.*W;
% Obtain the GCC
G = real(ifft(R));
G = fftshift(G,1);
figure(3);
plot(G);
[maxvalue2,shift2] = max(G);
shift2
title('GCC ROTH ');

%SCOT Filter
% ROTH Filter
W = 1./(10e-5+(Syy.*Sxx).^0.5);
R = Sxy.*W;
% Obtain the GCC
G = real(ifft(R));
G = fftshift(G,1);
figure(4);
plot(G);
[maxvalue3,shift3] = max(G);
shift3
title('GCC SCOT ');




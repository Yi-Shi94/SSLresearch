clear all
close all
ad = 'C:\YiSHI\AD1974Driver\Matlab\branches\wav\test_left.wav';
[x, fs] = audioread(ad);

micNum = size(x,2);
cpNum = (micNum-1)*micNum/2;

stsample = 1;
nsample = 512;
L = nsample-stsample+1;
T = 1/fs;
x = x(stsample:nsample,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%show spectrum of input sound channel and find prominent freq band%

s1 = x(stsample:nsample,2);
%s2 = x(stsample:nsample,2);
%s3 = x(stsample:nsample,3);
%s4 = x(stsample:nsample,4);

S1 = fft(s1,L);
t = (0:L-1)*T;   
P2 = abs(S1/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;

plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

[val,ind] = max(P1);
MajorFreq = fs/L*ind
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%show gcc-phat of channel1 and channel2 
G1 = gccphat(x(:,4), x(:,2), fs);
plot(G1)
[val,ind]= max(G1);
delay = (ind)/fs

tempC=24.0;
Vsound=331.4*sqrt(1.0+(tempC/273))
endFireDisMin = 0.015; %distance between two mics
endFireDisMinEst = Vsound * delay
ErrorDis = endFireDisMinEst - endFireDisMin 

gccMatrix = zeros(micNum,micNum,L*2);

for m = 1 : micNum
    for d = m+1 : micNum
        gccMatrix(m,d,:)=gccphat(x(:,m), x(:,d), fs);
    end
end

[val,ind] = max(gccMatrix);


 




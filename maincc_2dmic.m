clear all
close all
ad = 'C:\YiSHI\AD1974Driver\Matlab\branches\wav\test_right.wav';
[x, fs] = audioread(ad);

micNum = size(x,2);
cpNum = (micNum-1)*micNum/2;

stsample = 1;
nsample = 2^5;
L = nsample-stsample+1;
T = 1/fs;
x = x(stsample:2*nsample,:);

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
MajorFreq = fs/L*ind;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%show gcc-phat of channel1 and channel2 

tempC = 18.0;
Vsound = 331.4*sqrt(1.0+(tempC/273))
VsoundN = 300;
endFireDisMin = 0.015; %distance between two mics
%{
G1 = gccphat(x(:,1), x(:,2),fs);
[val,ind]= max(G1);
delay = (ind-L)/fs

endFireDisMinEst = Vsound * delay
ErrorDis = abs(abs(endFireDisMinEst) - endFireDisMin)

gccMatrix = zeros(micNum,micNum,L/2+1);
gccMaxMatrix = zeros(micNum,micNum);

for m = 1 : micNum
    for d = m+1 : micNum
        gccMatrix(m,d,:) = gccphat(x(:,m), x(:,d), fs);
        gccMaxMatrix(m,d) = max(gccMatrix(m,d,:));
    end
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%delay and sum...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freqDes = 16000;
freqSpan = 2000;
stride = 1;
disSearch = 1;
azNum = 180/stride + 1;
elNum = azNum;
mic_pos = [1 1 1; 1 1.015 1; 1 1.03 1; 1 1.045 1];
thfL = freqDes - freqSpan; %15khz-17khz
thfH = freqDes + freqSpan;
fdirSeq = 0.00001; 


s1 = x(stsample:nsample,1);
s2 = x(stsample:nsample,2);
s3 = x(stsample:nsample,3);
s4 = x(stsample:nsample,4);

X1 = myfft(s1);
X2 = myfft(s2);
X3 = myfft(s3);
X4 = myfft(s4);

dir  = zeros(4,1);
srpMatrix = zeros(fix(180/stride+1),fix(180/stride+1));
sa = size(X1)


for azAngle = -90:stride:90
    for elAngle = 0:stride:180
        
        
        azRadian = azAngle/180*pi;
        elRadian = elAngle/180*pi;
        x = disSearch*sin(elRadian)*cos(azRadian); 
        y = disSearch*sin(elRadian)*sin(azRadian); 
        z = disSearch*cos(elRadian); 
 
        dir(1,:) = (x*mic_pos(1,1)+y*mic_pos(1,2)+z*mic_pos(1,3))./Vsound; 
        dir(2,:) = (x*mic_pos(2,1)+y*mic_pos(2,2)+z*mic_pos(2,3))./Vsound+fdirSeq; 
        dir(3,:) = (x*mic_pos(3,1)+y*mic_pos(3,2)+z*mic_pos(3,3))./Vsound+fdirSeq*2; 
     	dir(4,:) = (x*mic_pos(4,1)+y*mic_pos(4,2)+z*mic_pos(4,3))./Vsound+fdirSeq*3; 
        p = 0; 
        s = 0; 
        for  frq=fix((thfL/fs)*L):fix((thfH/fs)*L)
            tmp = 2*pi*frq/(L/fs)*dir(1,:)*1j; 
     		s = s+X1(frq)*exp(tmp); 
 			tmp = 2*pi*frq/(L/fs)*dir(2,:)*1j; 
     		s = s+X2(frq)*exp(tmp); 
            tmp = 2*pi*frq/(L/fs)*dir(3,:)*1j; 
     		s = s+X3(frq)*exp(tmp); 
            tmp = 2*pi*frq/(L/fs)*dir(4,:)*1j; 
     		s = s+X4(frq)*exp(tmp); 
            p = p + abs(s); 
        end
        p = p^2;
        srpMatrix((azAngle+90)/stride+1,(elAngle)/stride+1)= srpMatrix((azAngle+90)/stride+1,(elAngle)/stride+1) + p; 
     end
end

surf(srpMatrix)

[rowsOfMaxes,colsOfMaxes] = find(srpMatrix == max(srpMatrix(:)));
azAngleMax = (rowsOfMaxes-1)*stride-90
elAngleMax = (colsOfMaxes-1)*stride
y



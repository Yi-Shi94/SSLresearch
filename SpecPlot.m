% anyway , we shall obtain the time difference here

%REQUIREMENT ON DISTANCE BETWEEN PAIR OF MICROPHONES
%sound speed in air
%vs = 343; 
%when T is 20'C

%vs = 331;
%when T is 0'C

%if max frequency equals fmax 
%fmax = 16000;
% d = 0.01;

%we can calculate max distance in m between pair of microphones
%dmax = vs/(2*fmax);
%fmax = vs/(2*d); 

function [Amp,Pha,f] = SpecPlot(y,Fs)
    [col,row] = size(y)
    T=1/Fs;
    %L = 10000;
    L = col;
    t =(1:L-1)*T;
    f = Fs*(0:(L/2));
    figure(3);
    for i = 1:row
        data = y(:,i);
        subplot(row/2,2,i);
        plot(f(1:(L/10)),data(1:(L/10)));
    end
    figure(1);
    
    for i = 1:row
        subplot(row/2,2,i);    
        ylabel('amplitude');
        xlabel('t/ms');
        F = fft(y(:,i));
        Amp = abs(F/L);
        Amp = Amp(1:L/2+1);
        Amp(2:end-1) = 2* Amp(2:end-1);
        plot(f(1:(L/10)),Amp(1:(L/10)));
    end
    figure(2)
    for i = 1:row
        subplot(row/2,2,i);
        F = fft(y(:,i));
        Pha = angle(F/L);
        Pha = Pha(1:L/2+1);
        Pha(2:end-1) = 2* Pha(2:end-1);
        plot(f(1:(L/10)),Pha(1:(L/10)));
    end
    
   pause
   
%}

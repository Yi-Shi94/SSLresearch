function G = gccphat(y1, y2, fs) 
    FFTLength= size(y1,1)+size(y2,1);
    Rxy = xcorr(y1,y2);
    Sxy = fft(Rxy,FFTLength);
    W = 1./(1.0e-4 + abs(Sxy));
    % Apply the filter
    R = Sxy.*W;
    % Obtain the GCC
    G = fftshift(real(ifft(R)),1);

    
    
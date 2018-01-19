clear all
close all
ad = 'C:\YiSHI\AD1974Driver\Matlab\branches\wav\test_left.wav';
[x, fs] = audioread(ad);

% s is a frame of audio signal, length of data from each channel is 2^12 = 4096 samples. 
% mic_pos is 3D position of the microphones.
% fs is sampling f

mic_pos = [0 0 0; 0.015 0 0; 0.03 0 0; 0.045 0 0];

bitlength = 12;
s = x(1:2^bitlength,:);
size(s)
usb = [16,16,16];
lsb = [0,0,0];

s1 = x(1:2^bitlength,1);
s3 = x(1:2^bitlength,3);
s2 = x(1:2^bitlength,2);
s4 = x(1:2^bitlength,4);

[finalpos,finalsrp] = srpgrid(s, mic_pos, fs, lsb, usb);

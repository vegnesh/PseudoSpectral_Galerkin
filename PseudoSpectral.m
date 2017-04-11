clc
clear all
close all

%% Test of Pseudo Spectral Methods
% Initial Data

TwoN = 100;
L = 2*pi;
xc = L/TwoN * (0:TwoN-1);

funct = @(x) exp(-(x-pi).^2);
dfunct= @(x) exp(-(x-pi).^2).*(-2*(x-pi));
%funct = @(x) sin(x);
%dfunct = @(x) cos(x);
yc = funct(xc);

%%

YN = fft(yc);

ycnew = ifft(YN);
% 
% IKfac = 0 + 1i*[0:TwoN/2-1 -TwoN/2:-1];
% DYN = IKfac.*YN;
% dyc = ifft(DYN);

dyc = Dfun(yc',L);

plot(xc,yc,'Linewidth',2);
hold on 
plot(xc,ycnew,'r')


figure(2)
plot(xc,dfunct(xc),'Linewidth',2)
hold on 
plot(xc,dyc,'r')

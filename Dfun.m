function [Dy] = Dfun(y,L)
%Dfun Obtains the Numerical derivative of a function 
%   The Numerical Derivative is obtained by using the Pseudo Spectral
%   Method

%% Sharanya Code
%Dy = derivative_PS(y,L,length(y)/2);

%%
y = y';
YN = fft(y);

N = length(y);
kvec = 1i*2*pi/L*[0:N/2-1 -N/2:-1];
DYN = kvec.*YN;
Dy = real(ifft(DYN));
Dy = Dy';
end


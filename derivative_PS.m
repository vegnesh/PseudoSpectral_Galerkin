% N- half the number of points
% L- length of spatial domain
% u - function whose PS derivative is required

function [deriv_PS]= derivative_PS (m, L, N)

%k=[0:N (1-N):-1]'*2*pi/L;       %wake number calculation
k=[0:N-1 -N:-1]'*2*pi/L;       %wake number calculation
%compute derivative
deriv_PS=real(ifft(1i*k.*fft(m)));


end
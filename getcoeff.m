function [Ak] = getcoeff(Fun,Nf, Nq)
%GETCOEFF Returns the fourier Coeff of the Fun in the domain
%  Takes in the function, number of modes and number of quadrature points
%  to use and returns the fourier coefficients
%   Fun - Function Handle
%   Nf -  Number of Modes
%   Nq - Number of Quadrature Points

%  Assumptions - Domain goes from 0 to 2pi and the mode coefficients go
%  from -Nf to Nf ie Ak(1) = -Nf_mode and Ak(end) = Nf mode and size of Ak
%  is 2*Nf+1

Ak = zeros(2*Nf+1,1);

xc = linspace(0,2*pi,Nq);

modefn = @(k,x) exp(-1i*k*x)/2/pi;
for ival = -Nf:Nf
    Ak(ival + Nf + 1) = trapz(xc,modefn(ival,xc).*Fun(xc));
end



end


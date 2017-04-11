function [ yval ] = getfunval(Ak,Nf,xc)
%GETFUNVAL Returns the value of the function at given x points
%   Takes in the fourier coefficients and returns the function value at the
%   points given by xc
%   Length of Ak should be 2Nf + 1

modefn = @(k,x) exp(1i*k*x);
yval = 0*xc;
for ival = -Nf:Nf
   yval = yval + Ak(ival+Nf+1)*modefn(ival,xc);
end

yval = real(yval);

end


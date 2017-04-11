function [Rf,Mf,Ef] = GetFlux(Qvec,N,Gam)
%GetFlux Returns the flux for conserved Euler Equations
%  Extracts the conserved or primitive variables from conserved variables
%  vector and returns the Flux Vectors
 
Rho = Qvec(1:N);
RU = Qvec(N+1:2*N);
RE = Qvec(2*N+1:end);

U = RU./Rho;
P = (Gam-1)*(RE - 0.5*RU.*U);

Rf = RU;
Mf = RU.*U + P;
Ef = (RE + P).*U;


end


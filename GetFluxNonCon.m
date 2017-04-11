function [Rf,Mf,Ef] = GetFluxNonCon(Qvec,N,Gam)
%GetFlux Returns the flux for conserved Euler Equations
%  Extracts the conserved or primitive variables from conserved variables
%  vector and returns the Flux Vectors
 
Rho = Qvec(1:N);
RU = Qvec(N+1:2*N);
RP = Qvec(2*N+1:end);

Rf = RU;
Mf = RU.*U + P;
Ef = (RE + P).*U;


end


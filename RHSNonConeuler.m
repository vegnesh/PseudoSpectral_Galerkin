function [Rhs] = RHSNonConeuler(Qvec,N,Gam,L)
%RHSEULER Returns the RHS of the Euler Equation 
%  This function returns the RHS of the Euler equation for the conservative
%  case based on the pseudo spectral method
%  Qvec -  Vector of Conserved Variables
%  N - Number of points in physical space
%  Gam  - Gamma value = 1.4
%  L - Length of the domain
%[A,B,C] = GetFluxNonCon(Qvec,N,Gam);


Rhs = 0*Qvec;
 
Rho = Qvec(1:N);
U = Qvec(N+1:2*N);
P = Qvec(2*N+1:end);

F1 = -Dfun(Rho.*U,L);
F2 = -Dfun(Rho.*U.^2 + P,L)./Rho - U.*F1./Rho;
F3 = -Dfun((P/(Gam-1) + 0.5*Rho.*U.^2 + P).*U,L)*(Gam-1) - 0.5*(Gam-1)*...
    (2*U.*(-Dfun(Rho.*U.^2 + P,L)) - U.^2.*F1);

Rhs(1:N) = F1;
Rhs(N+1:2*N) = F2;
Rhs(2*N+1:end) = F3;

end


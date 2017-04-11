function [Rhs] = RHSeuler(Qvec,N,Gam,L)
%RHSEULER Returns the RHS of the Euler Equation 
%  This function returns the RHS of the Euler equation for the conservative
%  case based on the pseudo spectral method
%  Qvec -  Vector of Conserved Variables
%  N - Number of points in physical space
%  Gam  - Gamma value = 1.4
%  L - Length of the domain
[A,B,C] = GetFlux(Qvec,N,Gam);

Rhs = zeros(length(Qvec),1);
Rhs(1:N) = Dfun(A,L);
Rhs(N+1:2*N) = Dfun(B,L);
Rhs(2*N+1:end) = Dfun(C,L);

Rhs = -Rhs;

end


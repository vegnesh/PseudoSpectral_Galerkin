function [Rhsvec] = RHSgalerkin(Qvec,Nf,Gam)
%RHSGALERKIN Returns the right of Euler Equations for Galerkin Method. 
%   Takes in Fourier Coefficients of 1/rho , u and P and returns the right
%   side of ODE. 
%   Note : Multiply by pi for account for change in x variable
%   Density  -u*d/dx (1/rho) + 1/rho*d/dx (u)
%   Velocity  -u*du/dx - 1/rho*dP/dx
%   Pressure  -Gam*P*du/dx - u*dP/dx

tempN = 2*Nf+1;

Rk = Qvec(1:tempN);
Uk = Qvec(tempN+1:2*tempN);
Pk = Qvec(2*tempN+1:end);

%Rhsvec = 0*Qvec;

Rk_Rhs = 0*Rk;
Uk_Rhs = 0*Uk;
Pk_Rhs = 0*Pk;

kvec = [-Nf:Nf]';

Rk_Rhs = conv(Uk,1i*kvec.*Rk,'same') - conv(Rk,1i*kvec.*Uk,'same');
Uk_Rhs = conv(Uk,1i*kvec.*Uk,'same') + conv(Rk,1i*kvec.*Pk,'same');
Pk_Rhs = conv(Gam*Pk,1i*kvec.*Uk,'same') + conv(Uk,1i*kvec.*Pk,'same');

Rhsvec = [Rk_Rhs;Uk_Rhs;Pk_Rhs];

Rhsvec = -pi*Rhsvec;


end


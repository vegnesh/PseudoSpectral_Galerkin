function [Rhsvec] = RHSgalerkin_old(Qvec,Nf,Gam)
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

Rhsvec = 0*Qvec;

Rk_Rhs = 0*Rk;
Uk_Rhs = 0*Uk;
Pk_Rhs = 0*Pk;

for ival = -Nf:Nf
    for kval = -Nf:Nf
        L_K = ((ival-kval)>=-Nf & (ival-kval)<=Nf);
        if(L_K)
        Rk_Rhs(ival + Nf + 1) = Rk_Rhs(ival + Nf + 1) + ...
            1i*kval*Rk(kval + Nf + 1)*Uk(ival - kval + Nf + 1) - ...
            1i*kval*Uk(kval + Nf + 1)*Rk(ival - kval + Nf + 1);
        
        Uk_Rhs(ival + Nf + 1) = Uk_Rhs(ival + Nf + 1) + ...
            1i*kval*Uk(kval + Nf + 1)*Uk(ival - kval + Nf + 1) + ...
            1i*kval*Pk(kval + Nf + 1)*Rk(ival - kval + Nf + 1);
        
        Pk_Rhs(ival + Nf + 1) = Pk_Rhs(ival + Nf + 1) + ...
            1i*kval*Pk(kval + Nf + 1)*Uk(ival - kval + Nf + 1) + ...
            1i*kval*Uk(kval + Nf + 1)*Pk(ival - kval + Nf + 1)*Gam;
        end
    end
end

Rhsvec = [Rk_Rhs;Uk_Rhs;Pk_Rhs];

Rhsvec = -pi*Rhsvec;


end


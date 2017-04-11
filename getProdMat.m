function [ Rhsvec ] = getProdMat(Ak,Nf)
%GETPRODMAT Returns the Matrix which is in the RHS
%   Detailed explanation goes here
Rhsvec = 0*Ak;
for ival = -Nf:Nf
    for kval = -Nf:Nf
        L_K = ((ival-kval)>=-Nf & (ival-kval)<=Nf);
        if(L_K)
        Rhsvec(ival + Nf + 1) = Rhsvec(ival + Nf + 1) + ...
            1i*kval*Ak(kval + Nf + 1)*Ak(ival - kval + Nf + 1);
        end
    end
end

Rhsvec = -Rhsvec;


end


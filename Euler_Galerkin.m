clc
%clear all 
close all
%clearvars -except Soln xc2
%% Euler Galerkin

u0 = 0.5;
Gam = 1.4;
c0 = 1;

ui = @(x) u0*sin(x);
Ti = @(x) (1+ 0.5*(Gam-1)*ui(x)/c0).^2;

Pi =@(x) Ti(x).^(Gam/(Gam-1));
R = 1/Gam;

rhoi_r =@(x) (R*Ti(x))./Pi(x);

%%
Nf = 90;
Nq = 3000;
Uk = getcoeff(ui,Nf,Nq);
Pk = getcoeff(Pi,Nf,Nq);
Rk = getcoeff(rhoi_r,Nf,Nq);

Qi = [Rk;Uk;Pk];

%% Time Integration

%Qsol = ode45(@(t,Q) RHSgalerkin(Q,Nf,Gam),[0,0.6],Qi);

% RK4 implementation

timetot = 0;
tf = 0.54;
dt = 1e-3;
Qfsoln = zeros(length(Qi),100);
index = 1;
Qc = Qi;
Qfsoln(:,index) = Qc;
while (timetot < tf)
    k1 = RHSgalerkin(Qc,Nf,Gam);
    k2 = RHSgalerkin(Qc + dt*k1*0.5,Nf,Gam);
    k3 = RHSgalerkin(Qc + dt*k2*0.5,Nf,Gam);
    k4 = RHSgalerkin(Qc + dt*k3,Nf,Gam);
    
    Qnew = Qc + dt*(k1+2*k2+k3*2+k4)/6;
    Qc = Qnew;
    index = index+1;
    Qfsoln(:,index) = Qc;
    timetot = timetot + dt;
end



%% Plotting
%Qfsoln = Qsol.y;

xc = linspace(0,2*pi,1000);
tempN = 2*Nf+1;
TwoN = 1000;
for ival = 1:10:length(Qfsoln(1,:))
    rk = Qfsoln(1:tempN,ival);
    uk = Qfsoln(tempN+1:2*tempN,ival);
    pk = Qfsoln(2*tempN+1:end,ival);
    
    %LRho = Soln(1:TwoN,20*(ival-1)+1);
    %LRU = Soln(TwoN+1:2*TwoN,20*(ival-1)+1);
    %LRE = Soln(2*TwoN+1:end,20*(ival-1)+1);

    %LU = LRU./LRho;
    %LP = (Gam-1)*(LRE - 0.5*LRU.*LU);
    %LT = LP./LRho/R;
    figure(1);
    plot(xc/pi,1./(getfunval(rk,Nf,xc)))
    hold on 
    %plot(xc2,LRho,'r--')
    hold off
    figure(2);
    plot(xc/pi,getfunval(uk,Nf,xc),'r')
    hold on 
    %plot(xc2,LU,'r--')
    hold off
    figure(3)
    plot(xc/pi,getfunval(pk,Nf,xc),'k')
    hold on
    %plot(xc2,LP,'r--')
    hold off
    
    pause(0.1)
    
    
end


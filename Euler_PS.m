clc
%clear all
close all

%% Euler Pseudo Spectral

TwoN = 100;
L = 2;
xc = L/TwoN * (0:TwoN-1);
xc = xc';
%% IC
u0 = 0.5;
Gam = 1.4;
c0 = 1;



ui = u0*sin(pi*xc);
Ti = (1+ 0.5*(Gam-1)*ui/c0).^2;

% Isentrope, when T = 1 P = 1

Pi = Ti.^(Gam/(Gam-1));
R = 1/Gam;

rhoi = Pi./(R*Ti);

%% Initial Conserved Variables Vector
Qveci = zeros(length(rhoi)*3,1);

Qveci(1:TwoN,1) = rhoi;
Qveci(TwoN+1:2*TwoN,1) = rhoi.*ui;
Qveci(2*TwoN+1:end,1) = Pi/(Gam-1) + 0.5*rhoi.*ui.^2;

%% Solution

%Qsol = ode45(@(t,Q) RHSeuler(Q,TwoN,Gam,L),[0,0.6],Qveci);
%Qsol = ode45(@(t,Q) RHSeuler(Q,TwoN,Gam,L),[0,0.53],Qveci);
%Qsol = ode45(@(t,Q) RHSeuler(Q,TwoN,Gam,L),[0,0.4],Qveci);

% RK4 implementation

timetot = 0;
tf = 0.54;
dt = 1e-4;
Qfsoln = zeros(length(Qveci),100);
index = 1;
Qc = Qveci;
Qfsoln(:,index) = Qc;
while (timetot < tf)
    k1 = RHSeuler(Qc,TwoN,Gam,L);
    k2 = RHSeuler(Qc + dt*k1*0.5,TwoN,Gam,L);
    k3 = RHSeuler(Qc + dt*k2*0.5,TwoN,Gam,L);
    k4 = RHSeuler(Qc + dt*k3,TwoN,Gam,L);
    
    Qnew = Qc + dt*(k1+2*k2+k3*2+k4)/6;
    Qc = Qnew;
    index = index+1;
    Qfsoln(:,index) = Qc;
    timetot = timetot + dt;
end



%% Plots
%Soln = Qsol.y;
Soln = Qfsoln;
Nsamples = length(Soln(1,:));
figure(1);
for ival = 1:5:Nsamples
    LRho = Soln(1:TwoN,ival);
    LRU = Soln(TwoN+1:2*TwoN,ival);
    LRE = Soln(2*TwoN+1:end,ival);

    LU = LRU./LRho;
    LP = (Gam-1)*(LRE - 0.5*LRU.*LU);
    LT = LP./LRho/R;
    
    %subplot(2,2,1)
    plot(xc,LRho,'Linewidth',2)
    hold on 
    %title('Density')
    %subplot(2,2,2)
    plot(xc,LP,'Linewidth',2)
    %title('Pressure')
    %subplot(2,2,3)
    plot(xc,LT,'Linewidth',2)
    %title('Temp')
    %subplot(2,2,4)
    plot(xc,LU,'Linewidth',2)
    hold off
    %title('Velocity')
    pause(0.001)
end

xc2 = xc;

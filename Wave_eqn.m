clc
clear all 
close all

%% Wave Equation Spectral


TwoN = 100;
L = 2*pi;
xc = L/TwoN * (0:TwoN-1);

funct = @(x) exp(-(x-pi).^2);
%dfunct= @(x) exp(-(x-pi).^2).*(-2*(x-pi));
%funct = @(x) sin(x);
%dfunct = @(x) cos(x);

u0 = funct(xc);
u0 = u0';

%% 1D Wave

rhsfn = @(t,u) (-Dfun(u,L));

u = ode45(rhsfn,[0:0.1:10],u0);

soln = u.y;

figure(1);

% for ival=1:length(soln(1,:))
%     
% 
%     plot(xc,soln(:,ival))
%     xlim([0,max(xc)]);
% ylim([min(u0),max(u0)]);
%     pause(0.1);
% end

%% 1D Burgers

rhsfn = @(t,u) (-u.*Dfun(u,L));

u = ode45(rhsfn,[0:0.1:1],u0);

soln = u.y;

figure(1);

for ival=1:length(soln(1,:))
    

    plot(xc,soln(:,ival))
    xlim([0,max(xc)]);
ylim([min(u0),max(u0)]);
    pause(0.1);
end

xc2 = xc;
